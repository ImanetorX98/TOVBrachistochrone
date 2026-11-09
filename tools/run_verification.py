#!/usr/bin/env python3
"""Run the TOVBrachistochrone verification suite and validate its results.

Usage:
    python3 tools/run_verification.py --profile quick
    python3 tools/run_verification.py --profile full --output-dir results/runs/my-run
    python3 tools/run_verification.py --profile release

Each job is a Wolfram Language script that writes a versioned JSON result to the path given in TOV_RESULT_JSON and
copies the run identifier TOV_RUN_ID into it. The runner decides on those JSON files, never on printed text.

Outcomes of a job:
    passed      exit status 0 and a valid result in which every check passed
    failed      exit status 1 and a valid result in which at least one check failed
    invalid     any other combination: no or malformed result, wrong run or job identifier, missing, extra or
                duplicate checks, wrong evidence kinds, inconsistent counts, abort (exit 2), crash or signal
    incomplete  timeout, or the job could not be launched

Exit status of the runner: 0 if every job passed; 1 if any job failed or was invalid; 2 for a usage or manifest
error (including an output directory that already holds files); 3 if the profile could not be completed (missing
engine, timeout, launch failure, unclean or unknown Git state for the release profile).
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from progress import SuiteProgress  # noqa: E402

SUMMARY_SCHEMA = "tov-summary/2"
RESULT_SCHEMA = "tov-result/1"
MANIFEST_VERSION = "tov-manifest/2"
KINDS = ("exact", "numerical_check", "diagnostic", "negative_control")
CATEGORIES = ("proof_support", "diagnostic")
SHARED_FILES = ["verification/common.wl", "verification/third_kind_lemma.wl"]
ID_PATTERN = re.compile(r"[a-z0-9_]+(\.[A-Za-z0-9_]+)*")
JOB_ID_PATTERN = re.compile(r"[a-z][a-z0-9_]{0,63}")  # used in file names (fullmatch): no separators, dots, ..
RESERVED_JOB_IDS = {"summary", "report", "jobs"}


class ManifestError(Exception):
    pass


# ---------------------------------------------------------------------------------------------- strict JSON

def _reject_duplicates(pairs: list) -> dict:
    keys = [k for k, _ in pairs]
    duplicated = sorted({k for k in keys if keys.count(k) > 1})
    if duplicated:
        raise ValueError(f"duplicate JSON keys: {duplicated}")
    return dict(pairs)


def _reject_constant(name: str):
    raise ValueError(f"non-standard JSON number {name}")


def parse_json_strict(text: str):
    return json.loads(text, object_pairs_hook=_reject_duplicates, parse_constant=_reject_constant)


def load_json_strict(path: Path):
    return parse_json_strict(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


# ---------------------------------------------------------------------------------------------- manifests

MANIFEST_FILES = ("verification/manifests/jobs.json", "verification/manifests/profiles.json")


def load_manifests(root: Path = ROOT, digests: dict | None = None) -> tuple[dict, dict]:
    """Load and validate the job and profile manifests; raise ManifestError on any defect.

    Each manifest is read once as bytes; if `digests` is given, the SHA-256 of exactly those bytes is stored in it,
    so the recorded digest identifies the manifest that was parsed and used."""
    docs = []
    for name in MANIFEST_FILES:
        try:
            raw = (root / name).read_bytes()
            docs.append(parse_json_strict(raw.decode("utf-8")))
        except (OSError, UnicodeDecodeError, ValueError) as error:
            raise ManifestError(f"cannot read {name}: {error}") from error
        if digests is not None:
            digests[name] = hashlib.sha256(raw).hexdigest()
    jobs_doc, profiles_doc = docs
    for name, doc in (("jobs.json", jobs_doc), ("profiles.json", profiles_doc)):
        if not isinstance(doc, dict) or doc.get("manifest_version") != MANIFEST_VERSION:
            raise ManifestError(f"{name}: manifest_version must be {MANIFEST_VERSION!r}")
    if not isinstance(jobs_doc.get("jobs"), list) or not jobs_doc["jobs"]:
        raise ManifestError("jobs.json: 'jobs' must be a non-empty list")
    jobs: dict[str, dict] = {}
    all_ids: set[str] = set()
    for job in jobs_doc["jobs"]:
        if not isinstance(job, dict):
            raise ManifestError("jobs.json: every job must be an object")
        required = {"job_id": str, "script": str, "engine": str, "category": str, "claims": list,
                    "expected_seconds": int, "timeout_seconds": int, "checks": list}
        for key, typ in required.items():
            if not isinstance(job.get(key), typ) or isinstance(job.get(key), bool):
                raise ManifestError(f"jobs.json: field {key!r} of {job.get('job_id')!r} must be {typ.__name__}")
        extra = set(job) - set(required)
        if extra:
            raise ManifestError(f"jobs.json: unexpected fields {sorted(extra)} in {job['job_id']!r}")
        jid = job["job_id"]
        if not JOB_ID_PATTERN.fullmatch(jid) or jid in RESERVED_JOB_IDS:
            raise ManifestError(f"jobs.json: job id {jid!r} is not a valid file-name stem ([a-z][a-z0-9_]*)")
        if not all(isinstance(c, str) for c in job["claims"]) or not job["claims"]:
            raise ManifestError(f"jobs.json: claims of {jid!r} must be a non-empty list of strings")
        if jid in jobs:
            raise ManifestError(f"jobs.json: duplicate job id {jid!r}")
        if job["engine"] != "wolframscript" or job["category"] not in CATEGORIES:
            raise ManifestError(f"jobs.json: unsupported engine or category in {jid!r}")
        if job["expected_seconds"] <= 0 or job["timeout_seconds"] <= 0:
            raise ManifestError(f"jobs.json: durations of {jid!r} must be positive")
        try:
            script = (root / job["script"]).resolve(strict=True)
        except (OSError, ValueError, RuntimeError) as error:
            raise ManifestError(f"jobs.json: script of {jid!r} cannot be resolved: {error}") from error
        if root.resolve() not in script.parents or not script.is_file():
            raise ManifestError(f"jobs.json: script of {jid!r} is missing or outside the repository")
        if not job["checks"]:
            raise ManifestError(f"jobs.json: {jid!r} has no expected checks")
        for check in job["checks"]:
            if not isinstance(check, dict) or set(check) != {"id", "kind", "claim"} or \
                    not all(isinstance(check[k], str) for k in check):
                raise ManifestError(f"jobs.json: malformed check entry in {jid!r}: {check!r}")
            if not ID_PATTERN.fullmatch(check["id"]) or check["kind"] not in KINDS:
                raise ManifestError(f"jobs.json: bad id or kind in {jid!r}: {check!r}")
            if check["claim"] not in job["claims"]:
                raise ManifestError(f"jobs.json: check {check['id']!r} refers to a claim not listed by {jid!r}")
            if check["id"] in all_ids:
                raise ManifestError(f"jobs.json: duplicate check id {check['id']!r}")
            all_ids.add(check["id"])
        jobs[jid] = job
    profiles = profiles_doc.get("profiles")
    if not isinstance(profiles, dict) or set(profiles) != {"quick", "full", "release"}:
        raise ManifestError("profiles.json: exactly the profiles quick, full and release are required")
    for name, profile in profiles.items():
        if not isinstance(profile, dict) or not isinstance(profile.get("jobs"), list) or not profile["jobs"]:
            raise ManifestError(f"profiles.json: profile {name!r} must list at least one job")
        if not all(isinstance(j, str) for j in profile["jobs"]):
            raise ManifestError(f"profiles.json: profile {name!r} must list job ids as strings")
        if len(set(profile["jobs"])) != len(profile["jobs"]):
            raise ManifestError(f"profiles.json: profile {name!r} lists a job twice")
        unknown = [j for j in profile["jobs"] if j not in jobs]
        if unknown:
            raise ManifestError(f"profiles.json: profile {name!r} lists unknown jobs {unknown}")
        if not isinstance(profile.get("description"), str):
            raise ManifestError(f"profiles.json: profile {name!r} needs a description")
    if set(profiles["full"]["jobs"]) != set(jobs) or set(profiles["release"]["jobs"]) != set(jobs):
        raise ManifestError("profiles.json: the full and release profiles must run every job")
    if profiles["release"].get("requires_clean_tree") is not True:
        raise ManifestError("profiles.json: the release profile must require a clean tree")
    return jobs, profiles


# ---------------------------------------------------------------------------------------------- results

def validate_result(job: dict, result_path: Path, exit_code: int | None, run_id: str | None = None) \
        -> tuple[str, list[str], dict]:
    """Return (outcome, problems, counts) with outcome 'passed', 'failed' or 'invalid'."""
    problems: list[str] = []
    if not result_path.is_file():
        return "invalid", ["no result JSON was written"], {}
    try:
        data = load_json_strict(result_path)
    except (OSError, UnicodeDecodeError, ValueError) as error:
        return "invalid", [f"result JSON is not valid: {error}"], {}
    if not isinstance(data, dict):
        return "invalid", ["the result JSON is not an object"], {}
    if data.get("schema_version") != RESULT_SCHEMA:
        problems.append(f"schema_version is {data.get('schema_version')!r}, expected {RESULT_SCHEMA!r}")
    if data.get("job_id") != job["job_id"]:
        problems.append(f"job_id is {data.get('job_id')!r}, expected {job['job_id']!r}")
    if run_id is not None and data.get("run_id") != run_id:
        problems.append(f"run_id is {data.get('run_id')!r}, not the identifier of this run")
    engine = data.get("engine")
    if not (isinstance(engine, dict) and isinstance(engine.get("name"), str) and isinstance(engine.get("version"), str)):
        problems.append("engine must be an object with string fields name and version")
    values = data.get("values")
    if not isinstance(values, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in values.items()):
        problems.append("values must be an object of strings")
    checks = data.get("checks")
    if not isinstance(checks, list) or not all(isinstance(c, dict) for c in checks):
        return "invalid", problems + ["checks must be a list of objects"], {}
    for check in checks:
        for field in ("id", "kind", "description", "status"):
            if not isinstance(check.get(field), str):
                problems.append(f"check {check.get('id')!r}: field {field!r} must be a string")
    if problems:
        return "invalid", problems, {}
    ids = [c["id"] for c in checks]
    expected = {c["id"]: c["kind"] for c in job["checks"]}
    duplicates = sorted({i for i in ids if ids.count(i) > 1})
    missing = sorted(set(expected) - set(ids))
    extra = sorted(set(ids) - set(expected))
    if duplicates:
        problems.append(f"duplicate check ids: {duplicates}")
    if missing:
        problems.append(f"missing check ids: {missing}")
    if extra:
        problems.append(f"unexpected check ids: {extra}")
    for check in checks:
        if check["id"] in expected and check["kind"] != expected[check["id"]]:
            problems.append(f"check {check['id']} has kind {check['kind']!r}, expected {expected[check['id']]!r}")
        if check["status"] not in ("pass", "fail"):
            problems.append(f"check {check['id']} has status {check['status']!r}")
    passed = sum(1 for c in checks if c["status"] == "pass")
    counts = data.get("counts")
    good_counts = {"total": len(checks), "passed": passed, "failed": len(checks) - passed}
    if not (isinstance(counts, dict) and set(counts) == set(good_counts) and
            all(type(counts[k]) is int and counts[k] == good_counts[k] for k in good_counts)):
        problems.append(f"counts {counts!r} disagree with the check list")
    all_pass = passed == len(checks) and len(checks) > 0
    if exit_code == 0 and not all_pass:
        problems.append("exit status 0 although some checks failed")
    elif exit_code == 1 and all_pass:
        problems.append("exit status 1 although every check passed")
    elif exit_code not in (0, 1, None):
        problems.append(f"unexpected exit status {exit_code}")
    by_kind: dict[str, dict[str, int]] = {}
    for check in checks:
        entry = by_kind.setdefault(check["kind"], {"pass": 0, "fail": 0})
        entry["pass" if check["status"] == "pass" else "fail"] += 1
    summary_counts = {"total": len(checks), "passed": passed, "by_kind": by_kind,
                      "engine_version": engine["version"]}
    if problems:
        return "invalid", problems, summary_counts
    return ("passed" if all_pass else "failed"), [], summary_counts


def _tail(path: Path, lines: int = 20) -> list[str]:
    try:
        return path.read_text(encoding="utf-8", errors="replace").splitlines()[-lines:]
    except OSError:
        return []


def run_job(job: dict, out_dir: Path, engine: str, run_id: str, progress: SuiteProgress) -> dict:
    script = ROOT / job["script"]
    result_path = out_dir / "jobs" / f"{job['job_id']}.json"
    log_path = out_dir / "jobs" / f"{job['job_id']}.stdout.log"
    err_path = out_dir / "jobs" / f"{job['job_id']}.stderr.log"
    env = dict(os.environ, TOV_RESULT_JSON=str(result_path), TOV_RUN_ID=run_id)
    progress.job_started(job["job_id"])
    started = time.monotonic()
    exit_code: int | None = None
    outcome, problems, counts = "invalid", [], {}
    try:
        with log_path.open("w", encoding="utf-8") as out, err_path.open("w", encoding="utf-8") as err:
            process = subprocess.Popen([engine, "-file", str(script)], stdin=subprocess.DEVNULL, stdout=out,
                                       stderr=err, env=env, cwd=ROOT, start_new_session=True)
            try:
                exit_code = process.wait(timeout=job["timeout_seconds"])
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
                outcome, problems = "incomplete", [f"timeout after {job['timeout_seconds']} s; process group killed"]
    except OSError as error:
        outcome, problems = "incomplete", [f"the job could not be launched: {error}"]
    if outcome != "incomplete":
        if exit_code == 2:
            outcome, problems = "invalid", ["the script aborted on a violated assumption (exit status 2)"]
        elif exit_code is not None and exit_code < 0:
            outcome, problems = "invalid", [f"the engine was terminated by signal {-exit_code}"]
        elif exit_code not in (0, 1):
            outcome, problems = "invalid", [f"unexpected exit status {exit_code}"]
        else:
            outcome, problems, counts = validate_result(job, result_path, exit_code, run_id)
    progress.job_finished(job["job_id"], outcome)
    return {"job_id": job["job_id"], "script": job["script"], "category": job["category"],
            "claims": job["claims"], "script_sha256": sha256(script), "exit_code": exit_code,
            "elapsed_seconds": round(time.monotonic() - started, 2), "outcome": outcome, "problems": problems,
            "counts": counts, "result_json": f"jobs/{result_path.name}",
            "stdout_log": f"jobs/{log_path.name}", "stderr_log": f"jobs/{err_path.name}", "stderr_tail": _tail(err_path)}


# ---------------------------------------------------------------------------------------------- Git state

def git_state(root: Path = ROOT) -> dict:
    """Fail closed: 'ok' is True only if every Git command succeeded and the repository root is `root`."""
    def git(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
    try:
        top = git("rev-parse", "--show-toplevel")
        head = git("rev-parse", "--verify", "HEAD")
        branch = git("rev-parse", "--abbrev-ref", "HEAD")
        status = git("status", "--porcelain", "--untracked-files=normal")
    except OSError as error:
        return {"ok": False, "clean": False, "problem": f"git could not be run: {error}"}
    failed = [c for c in (top, head, branch, status) if c.returncode != 0]
    if failed:
        return {"ok": False, "clean": False, "problem": failed[0].stderr.strip() or "a git command failed"}
    if Path(top.stdout.strip()).resolve() != root.resolve():
        return {"ok": False, "clean": False, "problem": f"the Git root is {top.stdout.strip()}, not {root}"}
    changes = status.stdout.splitlines()
    return {"ok": True, "commit": head.stdout.strip(), "branch": branch.stdout.strip(),
            "clean": not changes, "changes": changes}


# ---------------------------------------------------------------------------------------------- main

def positive_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError("must be a positive finite number")
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run and validate the TOVBrachistochrone verification suite.")
    parser.add_argument("--profile", choices=["quick", "full", "release"], default="quick")
    parser.add_argument("--output-dir", type=Path,
                        help="new or empty directory for the results (default: results/runs/<time>-<profile>-<id>)")
    parser.add_argument("--jobs", type=int, default=2, help="number of jobs run in parallel (default 2)")
    parser.add_argument("--progress", choices=["auto", "always", "never"], default="auto")
    parser.add_argument("--heartbeat", type=positive_float, default=15.0,
                        help="seconds between plain-text snapshots")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")

    try:
        manifest_digests: dict = {}
        jobs, profiles = load_manifests(digests=manifest_digests)
    except ManifestError as error:
        print(f"MANIFEST ERROR: {error}", file=sys.stderr)
        return 2
    profile = profiles[args.profile]
    selected = [jobs[j] for j in profile["jobs"]]

    run_id = uuid.uuid4().hex
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    try:
        out_dir = (args.output_dir if args.output_dir else
                   ROOT / "results" / "runs" / f"{stamp}-{args.profile}-{run_id[:8]}").expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        # atomic reservation: exactly one process can create the marker; only then is emptiness checked
        marker = out_dir / ".run_id"
        try:
            fd = os.open(marker, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            print(f"USAGE ERROR: the output directory {out_dir} is already used by another run", file=sys.stderr)
            return 2
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(run_id + "\n")
        if any(p.name != ".run_id" for p in out_dir.iterdir()):
            marker.unlink()
            print(f"USAGE ERROR: the output directory {out_dir} is not empty; use a new directory", file=sys.stderr)
            return 2
        (out_dir / "jobs").mkdir()
    except (OSError, ValueError, RuntimeError) as error:
        print(f"USAGE ERROR: cannot prepare the output directory: {error}", file=sys.stderr)
        return 2

    summary: dict = {"schema_version": SUMMARY_SCHEMA, "run_id": run_id, "profile": args.profile,
                     "profile_description": profile["description"], "started_utc": stamp,
                     "python": platform.python_version(), "platform": platform.platform(),
                     "git": git_state(), "jobs": [], "status": "incomplete", "problems": []}
    write_json_atomic(out_dir / "summary.json", summary)

    engine = shutil.which("wolframscript")
    if engine is None:
        summary["problems"].append("wolframscript was not found; the profile cannot be completed")
    elif profile.get("requires_clean_tree"):
        git = summary["git"]
        if not git.get("ok"):
            summary["problems"].append(f"the release profile needs a valid Git checkout: {git.get('problem')}")
        elif not git["clean"]:
            summary["problems"].append("the release profile requires a clean working tree (including untracked files)")
    if summary["problems"]:
        write_json_atomic(out_dir / "summary.json", summary)
        for problem in summary["problems"]:
            print(f"INCOMPLETE: {problem}", file=sys.stderr)
        return 3
    summary["engine"] = {"path": engine}
    files = sorted({j["script"] for j in selected} | set(SHARED_FILES) |
                   {"tools/run_verification.py", "tools/progress.py"})
    summary["manifest_sha256_parsed"] = manifest_digests
    summary["sha256_before"] = {f: sha256(ROOT / f) for f in files}
    summary["sha256_before"].update(manifest_digests)

    progress = SuiteProgress({j["job_id"]: float(j["expected_seconds"]) for j in selected},
                             mode=args.progress, interval=args.heartbeat)
    progress.start()
    results = []
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run_job, job, out_dir, engine, run_id, progress): job for job in selected}
            for future, job in futures.items():
                try:
                    results.append(future.result())
                except Exception as error:  # a crash of the runner itself is reported, never hidden
                    results.append({"job_id": job["job_id"], "script": job["script"], "category": job["category"],
                                    "claims": job["claims"], "outcome": "invalid",
                                    "problems": [f"runner error: {type(error).__name__}: {error}"], "counts": {},
                                    "elapsed_seconds": 0.0, "exit_code": None})
    finally:
        progress.close()

    summary["jobs"] = results
    summary["sha256_after"] = {f: sha256(ROOT / f) for f in list(files) + list(MANIFEST_FILES)}
    if summary["sha256_after"] != summary["sha256_before"]:
        summary["problems"].append("verification files changed while the suite was running")
    versions = sorted({r["counts"].get("engine_version") for r in results if r.get("counts")} - {None})
    summary["engine"]["versions_reported"] = versions
    if len(versions) > 1:
        summary["problems"].append(f"jobs reported different engine versions: {versions}")
    outcomes = [r["outcome"] for r in results]
    if any(o == "incomplete" for o in outcomes):
        summary["status"] = "incomplete"
    elif summary["problems"] or any(o != "passed" for o in outcomes):
        summary["status"] = "failed"
    else:
        summary["status"] = "passed"
    summary["finished_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    git = summary["git"]
    head = (f"commit {git['commit'][:12]} on {git['branch']}" + ("" if git["clean"] else
            " (working tree has changes)")) if git.get("ok") else f"Git state unavailable: {git.get('problem')}"
    lines = [f"TOVBrachistochrone verification, profile '{args.profile}': {summary['status'].upper()}",
             f"run {run_id}; {head}", ""]
    total_checks, diagnostic_checks = 0, 0
    for r in results:
        by_kind = r["counts"].get("by_kind", {}) if r.get("counts") else {}
        total_checks += sum(v["pass"] + v["fail"] for v in by_kind.values())
        diagnostic_checks += sum(v["pass"] + v["fail"] for k, v in by_kind.items() if k == "diagnostic")
        kinds = ", ".join(f"{k} {v['pass']}/{v['pass'] + v['fail']}" for k, v in sorted(by_kind.items()))
        lines.append(f"  {r['outcome']:<10} {r['job_id']:<28} {r['elapsed_seconds']:>7.1f} s   {kinds}")
        lines.extend(f"             problem: {p}" for p in r["problems"])
    lines.extend(f"  problem: {p}" for p in summary["problems"])
    lines += ["", f"{total_checks} checks in {len(results)} jobs, of which {diagnostic_checks} are diagnostics "
              "(sampled statements, not proofs). A count of checks is not a count of theorems.",
              "Not executed by this suite: the analytic proof steps listed in RESULTS_MAP.md.", f"Results: {out_dir}"]
    report = "\n".join(lines) + "\n"
    (out_dir / "report.txt").write_text(report, encoding="utf-8")
    # the terminal summary is the last persistent artefact: a green summary.json attests that the jobs ran and that
    # report.txt and every job artefact were written before it
    write_json_atomic(out_dir / "summary.json", summary)
    print(report, end="")
    return {"passed": 0, "failed": 1, "incomplete": 3}[summary["status"]]


if __name__ == "__main__":
    sys.exit(main())
