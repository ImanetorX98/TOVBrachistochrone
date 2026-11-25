# Results

Run outputs of `tools/run_verification.py` go to `results/runs/<UTC time>-<profile>-<run id>/` (not tracked by git):
`.run_id` (the reservation marker), `jobs/` (one JSON result, one stdout log and one stderr log per job),
`report.txt`, and `summary.json`, which is written last.

## Release reports

Each `results/release/vX.Y.Z/` (`v2.1.0/` current, `v2.0.0/` previous) contains `summary.json` and `report.txt` of the
release-profile run made on a fresh clone of the release-candidate commit recorded in `summary.json` (field `git`). The only later change before the tag is the
addition of this report; neither the per-job logs nor the per-job JSON results of that run are
tracked, so the report does not keep every digit of the individual jobs. In `report.txt` the local output path of the run
is redacted; nothing else is changed.
