#!/usr/bin/env python3
"""Tests of run_verification.py (standard library only; no Wolfram needed).

Unit tests corrupt a valid result or manifest in one way and check that the validator reports exactly that defect.
Integration tests run the real runner on a temporary copy of the repository whose jobs are executed by a fake
`wolframscript`, to exercise the run lifecycle: exclusive output directories, run identifiers, exit statuses,
signals, timeouts with child processes, relative paths, Git state, invalid manifests and a missing engine.

    python3 tools/test_run_verification.py
"""

from __future__ import annotations

import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
from run_verification import ManifestError, load_manifests, validate_result  # noqa: E402

JOB = {"job_id": "example_job", "checks": [{"id": "example.a", "kind": "exact", "claim": "c.one"},
                                           {"id": "example.b", "kind": "numerical_check", "claim": "c.one"}]}
RESULT = {"schema_version": "tov-result/1", "job_id": "example_job", "run_id": "r1",
          "engine": {"name": "fake", "version": "1"},
          "checks": [{"id": "example.a", "kind": "exact", "description": "", "status": "pass"},
                     {"id": "example.b", "kind": "numerical_check", "description": "", "status": "pass"}],
          "values": {}, "counts": {"total": 2, "passed": 2, "failed": 0}}


class ValidatorTests(unittest.TestCase):
    def validate(self, result, job=JOB, exit_code=0, run_id="r1"):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "result.json"
            if result is not None:
                path.write_text(result if isinstance(result, str) else json.dumps(result), encoding="utf-8")
            return validate_result(job, path, exit_code, run_id)

    def assertProblem(self, outcome_problems, fragment):
        outcome, problems, _ = outcome_problems
        self.assertEqual(outcome, "invalid")
        self.assertTrue(any(fragment in p for p in problems), problems)

    def test_valid_result_passes(self):
        self.assertEqual(self.validate(RESULT)[0], "passed")

    def test_failed_check_is_failed_not_invalid(self):
        bad = copy.deepcopy(RESULT)
        bad["checks"][1]["status"] = "fail"
        bad["counts"] = {"total": 2, "passed": 1, "failed": 1}
        self.assertEqual(self.validate(bad, exit_code=1)[0], "failed")

    def test_missing_file(self):
        self.assertProblem(self.validate(None), "no result JSON")

    def test_unparsable_file(self):
        self.assertProblem(self.validate("{not json"), "not valid")

    def test_duplicate_key(self):
        text = json.dumps(RESULT)[:-1] + ', "job_id": "example_job"}'
        self.assertProblem(self.validate(text), "duplicate JSON keys")

    def test_nan(self):
        text = json.dumps(RESULT).replace('"values": {}', '"values": {"x": NaN}')
        self.assertProblem(self.validate(text), "non-standard JSON number")

    def test_non_object_root(self):
        self.assertProblem(self.validate("[1, 2]"), "not an object")

    def test_checks_not_a_list(self):
        self.assertProblem(self.validate(dict(RESULT, checks=None)), "checks must be a list")

    def test_missing_engine_values_description(self):
        bad = copy.deepcopy(RESULT)
        del bad["engine"], bad["values"], bad["checks"][0]["description"]
        outcome, problems, _ = self.validate(bad)
        self.assertEqual(outcome, "invalid")
        self.assertTrue(any("engine" in p for p in problems) and any("values" in p for p in problems), problems)

    def test_mixed_type_ids(self):
        bad = copy.deepcopy(RESULT)
        bad["checks"].append({"id": 3, "kind": "exact", "description": "", "status": "pass"})
        self.assertProblem(self.validate(bad), "field 'id' must be a string")

    def test_boolean_counts(self):
        bad = dict(RESULT, counts={"total": 2, "passed": 2, "failed": False})
        self.assertProblem(self.validate(bad), "counts")

    def test_wrong_schema(self):
        self.assertProblem(self.validate(dict(RESULT, schema_version="tov-result/0")), "schema_version")

    def test_wrong_job_id(self):
        self.assertProblem(self.validate(dict(RESULT, job_id="other")), "job_id")

    def test_wrong_run_id(self):
        self.assertProblem(self.validate(dict(RESULT, run_id="another-run")), "run_id")

    def test_missing_check(self):
        bad = copy.deepcopy(RESULT)
        bad["checks"].pop()
        bad["counts"] = {"total": 1, "passed": 1, "failed": 0}
        self.assertProblem(self.validate(bad), "missing check ids")

    def test_extra_check(self):
        job = copy.deepcopy(JOB)
        job["checks"].pop()
        self.assertProblem(self.validate(RESULT, job=job), "unexpected check ids")

    def test_duplicate_check(self):
        bad = copy.deepcopy(RESULT)
        bad["checks"].append(copy.deepcopy(bad["checks"][0]))
        bad["counts"] = {"total": 3, "passed": 3, "failed": 0}
        self.assertProblem(self.validate(bad), "duplicate check ids")

    def test_wrong_kind(self):
        bad = copy.deepcopy(RESULT)
        bad["checks"][0]["kind"] = "diagnostic"
        self.assertProblem(self.validate(bad), "has kind")

    def test_inconsistent_counts(self):
        self.assertProblem(self.validate(dict(RESULT, counts={"total": 2, "passed": 1, "failed": 1})), "counts")

    def test_exit_status_one_with_all_pass(self):
        self.assertProblem(self.validate(RESULT, exit_code=1), "exit status 1")

    def test_exit_status_zero_with_failure(self):
        bad = copy.deepcopy(RESULT)
        bad["checks"][1]["status"] = "fail"
        bad["counts"] = {"total": 2, "passed": 1, "failed": 1}
        self.assertProblem(self.validate(bad, exit_code=0), "exit status 0")


# ------------------------------------------------------------------------------------------------ integration

FAKE_ENGINE = textwrap.dedent('''\
    #!/usr/bin/env python3
    """Fake wolframscript: the 'script' is a JSON file describing what to do."""
    import json, os, signal, subprocess, sys, time
    spec = json.load(open(sys.argv[sys.argv.index("-file") + 1]))
    out, run_id = os.environ["TOV_RESULT_JSON"], spec.get("run_id", os.environ.get("TOV_RUN_ID"))
    if spec.get("child_sleep"):
        subprocess.Popen([sys.executable, "-c", "import time; time.sleep(%d)" % spec["child_sleep"]])
    time.sleep(spec.get("sleep", 0))
    status = spec.get("status", "pass")
    if spec.get("write", True):
        checks = [{"id": c, "kind": "exact", "description": "", "status": status} for c in spec["ids"]]
        passed = sum(c["status"] == "pass" for c in checks)
        json.dump({"schema_version": "tov-result/1", "job_id": spec["job_id"], "run_id": run_id,
                   "engine": {"name": "fake", "version": "1"}, "checks": checks, "values": {},
                   "counts": {"total": len(checks), "passed": passed, "failed": len(checks) - passed}}, open(out, "w"))
    if spec.get("touch_manifest"):
        path = os.path.join(os.path.dirname(os.path.dirname(spec["touch_manifest"])), "manifests", "profiles.json")
        with open(path, "a") as handle:
            handle.write(" ")
    if spec.get("signal"):
        os.kill(os.getpid(), signal.SIGKILL)
    sys.exit(spec.get("exit", 0 if status == "pass" else 1))
''')


class RunnerIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.repo = self.tmp / "repo"
        (self.repo / "tools").mkdir(parents=True)
        for name in ("run_verification.py", "progress.py"):
            shutil.copy(TOOLS / name, self.repo / "tools" / name)
        for name in ("verification/common.wl", "verification/third_kind_lemma.wl"):
            (self.repo / name).parent.mkdir(parents=True, exist_ok=True)
            (self.repo / name).write_text("(* placeholder *)\n")
        (self.repo / "verification/manifests").mkdir(parents=True)
        self.bin = self.tmp / "bin"
        self.bin.mkdir()
        engine = self.bin / "wolframscript"
        engine.write_text(FAKE_ENGINE)
        engine.chmod(0o755)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write_jobs(self, specs: dict, timeout: int = 30, profiles: dict | None = None):
        jobs = []
        for jid, spec in specs.items():
            script = self.repo / "verification" / f"{jid}.json"
            script.write_text(json.dumps(dict(spec, job_id=jid, ids=[f"{jid}.check"])))
            jobs.append({"job_id": jid, "script": f"verification/{jid}.json", "engine": "wolframscript",
                         "category": "proof_support", "claims": ["c.one"], "expected_seconds": 1,
                         "timeout_seconds": timeout,
                         "checks": [{"id": f"{jid}.check", "kind": "exact", "claim": "c.one"}]})
        (self.repo / "verification/manifests/jobs.json").write_text(
            json.dumps({"manifest_version": "tov-manifest/2", "jobs": jobs}))
        names = list(specs)
        profiles = profiles or {"quick": {"description": "q", "jobs": names},
                                "full": {"description": "f", "jobs": names},
                                "release": {"description": "r", "jobs": names, "requires_clean_tree": True}}
        (self.repo / "verification/manifests/profiles.json").write_text(
            json.dumps({"manifest_version": "tov-manifest/2", "profiles": profiles}))

    def run_runner(self, *args, path_with_engine=True, cwd=None):
        env = dict(os.environ)
        env["PATH"] = (str(self.bin) + os.pathsep if path_with_engine else "") + "/usr/bin:/bin"
        return subprocess.run([sys.executable, str(self.repo / "tools/run_verification.py"), "--progress", "never",
                               *args], capture_output=True, text=True, env=env, cwd=cwd or self.tmp, timeout=120)

    def test_passing_run(self):
        self.write_jobs({"a": {}, "b": {}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 0)
        summary = json.loads((self.tmp / "out" / "summary.json").read_text())
        for job in summary["jobs"]:
            result = json.loads((self.tmp / "out" / job["result_json"]).read_text())
            self.assertEqual((result["schema_version"], result["job_id"], result["run_id"]),
                             ("tov-result/1", job["job_id"], summary["run_id"]))

    def test_manifest_changed_during_the_run_fails(self):
        self.write_jobs({"a": {"touch_manifest": str(self.repo / "verification/a.json")}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 1)

    def test_output_directory_must_be_new_or_empty(self):
        self.write_jobs({"a": {}})
        (self.tmp / "out").mkdir()
        (self.tmp / "out" / "old.json").write_text("{}")
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 2)

    def test_relative_output_directory_is_resolved_from_the_caller(self):
        self.write_jobs({"a": {}})
        caller = self.tmp / "caller"
        caller.mkdir()
        result = self.run_runner("--output-dir", "rel", cwd=caller)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((caller / "rel" / "jobs" / "a.json").is_file())

    def test_result_of_another_run_is_rejected(self):
        self.write_jobs({"a": {"run_id": "someone-else"}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 1)

    def test_signal_with_a_failure_json_is_invalid(self):
        self.write_jobs({"a": {"status": "fail", "signal": True}})
        result = self.run_runner("--output-dir", str(self.tmp / "out"))
        summary = json.loads((self.tmp / "out" / "summary.json").read_text())
        self.assertEqual(result.returncode, 1)
        self.assertEqual(summary["jobs"][0]["outcome"], "invalid")

    def test_abort_exit_two_is_invalid(self):
        self.write_jobs({"a": {"exit": 2, "write": False}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 1)

    def test_silent_exit_zero_without_result_is_invalid(self):
        self.write_jobs({"a": {"write": False}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 1)

    def test_timeout_kills_the_process_group(self):
        self.write_jobs({"a": {"sleep": 30, "child_sleep": 30}}, timeout=2)
        started = time.monotonic()
        result = self.run_runner("--output-dir", str(self.tmp / "out"))
        self.assertEqual(result.returncode, 3)
        self.assertLess(time.monotonic() - started, 20)

    def test_release_outside_git_is_incomplete(self):
        self.write_jobs({"a": {}})
        self.assertEqual(self.run_runner("--profile", "release", "--output-dir", str(self.tmp / "out")).returncode, 3)

    def test_empty_profile_is_a_manifest_error(self):
        self.write_jobs({"a": {}}, profiles={"quick": {"description": "q", "jobs": []},
                                             "full": {"description": "f", "jobs": ["a"]},
                                             "release": {"description": "r", "jobs": ["a"], "requires_clean_tree": True}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 2)

    def test_missing_engine_is_incomplete(self):
        self.write_jobs({"a": {}})
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out"), path_with_engine=False).returncode, 3)

    def test_duplicate_check_id_in_manifest(self):
        self.write_jobs({"a": {}})
        doc = json.loads((self.repo / "verification/manifests/jobs.json").read_text())
        doc["jobs"][0]["checks"].append(dict(doc["jobs"][0]["checks"][0]))
        (self.repo / "verification/manifests/jobs.json").write_text(json.dumps(doc))
        with self.assertRaises(ManifestError):
            load_manifests(self.repo)


    def test_directory_reserved_by_another_run_is_refused(self):
        self.write_jobs({"a": {}})
        (self.tmp / "out").mkdir()
        (self.tmp / "out" / ".run_id").write_text("another-run\n")
        self.assertEqual(self.run_runner("--output-dir", str(self.tmp / "out")).returncode, 2)

    def test_concurrent_runs_on_one_directory(self):
        self.write_jobs({"a": {"sleep": 2}})
        env = dict(os.environ, PATH=str(self.bin) + os.pathsep + "/usr/bin:/bin")
        cmd = [sys.executable, str(self.repo / "tools/run_verification.py"), "--progress", "never",
               "--output-dir", str(self.tmp / "shared")]
        procs = [subprocess.Popen(cmd, env=env, cwd=self.tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 for _ in range(4)]
        for proc in procs:
            proc.communicate(timeout=120)
        codes = sorted(p.returncode for p in procs)
        self.assertEqual(codes, [0, 2, 2, 2])

    def _manifest_error(self, mutate):
        self.write_jobs({"a": {}})
        for name in ("jobs.json", "profiles.json"):
            path = self.repo / "verification/manifests" / name
            doc = json.loads(path.read_text())
            mutate(name, doc)
            path.write_text(json.dumps(doc))
        with self.assertRaises(ManifestError):
            load_manifests(self.repo)

    def test_job_id_with_path_separators(self):
        self._manifest_error(lambda n, d: d["jobs"][0].update(job_id="../../escape") if n == "jobs.json" else None)

    def test_nested_list_in_profile(self):
        self._manifest_error(lambda n, d: d["profiles"]["quick"].update(jobs=[["a"]]) if n == "profiles.json" else None)

    def test_reserved_job_id(self):
        self._manifest_error(lambda n, d: d["jobs"][0].update(job_id="summary") if n == "jobs.json" else None)

    def test_job_id_with_trailing_newline(self):
        self._manifest_error(lambda n, d: d["jobs"][0].update(job_id="a\n") if n == "jobs.json" else None)

    def test_check_id_with_trailing_newline(self):
        self._manifest_error(lambda n, d: d["jobs"][0]["checks"][0].update(id="a.check\n") if n == "jobs.json" else None)

    def test_script_path_with_nul(self):
        self._manifest_error(lambda n, d: d["jobs"][0].update(script="verification/a\u0000.json") if n == "jobs.json" else None)

    def test_non_string_claim(self):
        self._manifest_error(lambda n, d: d["jobs"][0].update(claims=[["c.one"]]) if n == "jobs.json" else None)


class RepositoryManifestTests(unittest.TestCase):
    def test_repository_manifests_are_valid_and_claims_are_mapped(self):
        root = TOOLS.parent
        jobs, profiles = load_manifests(root)
        result_map = (root / "RESULTS_MAP.md").read_text(encoding="utf-8")
        claims = {c for job in jobs.values() for c in job["claims"]}
        missing = sorted(c for c in claims if f"`{c}`" not in result_map)
        self.assertEqual(missing, [], "claims used by the manifest but absent from RESULTS_MAP.md")


if __name__ == "__main__":
    unittest.main(verbosity=2)
