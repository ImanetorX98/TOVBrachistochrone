"""Progress reporting for slow verification jobs.

Progress goes to stderr only and never into the result files. In an interactive terminal a single status line is
redrawn; otherwise plain-text snapshots are printed at a fixed interval. The estimated time remaining is derived
from the expected durations in the job manifest and is labelled as an estimate.
"""

from __future__ import annotations

import sys
import threading
import time


def _format_seconds(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    minutes, sec = divmod(seconds, 60)
    return f"{minutes}m{sec:02d}s" if minutes else f"{sec}s"


class SuiteProgress:
    """Aggregated progress of a verification suite, safe to update from worker threads."""

    def __init__(self, weights: dict[str, float], mode: str = "auto", interval: float = 15.0) -> None:
        self._weights = weights
        self._total_weight = sum(weights.values()) or 1.0
        self._running: dict[str, float] = {}
        self._finished: dict[str, str] = {}
        self._lock = threading.Lock()
        self._start = time.monotonic()
        self._interval = interval
        self._interactive = sys.stderr.isatty()
        self._enabled = mode != "never"
        self._redraw = self._enabled and self._interactive
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def start(self) -> None:
        if not self._enabled:
            return
        self._thread = threading.Thread(target=self._heartbeat, daemon=True)
        self._thread.start()

    def job_started(self, job_id: str) -> None:
        with self._lock:
            self._running[job_id] = time.monotonic()
        self._emit(f"started {job_id}")

    def job_finished(self, job_id: str, outcome: str) -> None:
        with self._lock:
            started = self._running.pop(job_id, time.monotonic())
            self._finished[job_id] = outcome
        self._emit(f"finished {job_id}: {outcome} in {_format_seconds(time.monotonic() - started)}")

    def close(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)
        if self._redraw:
            sys.stderr.write("\n")
            sys.stderr.flush()

    def _status(self) -> str:
        with self._lock:
            done_weight = sum(self._weights.get(j, 0.0) for j in self._finished)
            running = sorted(self._running)
            n_done, n_total = len(self._finished), len(self._weights)
        elapsed = time.monotonic() - self._start
        fraction = done_weight / self._total_weight
        if 0 < fraction < 1:
            eta = f"ETA ~{_format_seconds(elapsed * (1 - fraction) / fraction)} (estimate)"
        else:
            eta = "ETA unavailable"
        current = ", ".join(running) if running else "-"
        return f"[{n_done}/{n_total} jobs] elapsed {_format_seconds(elapsed)}, {eta}; running: {current}"

    def _emit(self, event: str) -> None:
        if not self._enabled:
            return
        if self._redraw:
            sys.stderr.write("\r\x1b[K" + event + "\n" + self._status())
        else:
            sys.stderr.write(event + "\n")
        sys.stderr.flush()

    def _heartbeat(self) -> None:
        period = 1.0 if self._redraw else self._interval
        while not self._stop.wait(period):
            if self._redraw:
                sys.stderr.write("\r\x1b[K" + self._status())
            else:
                sys.stderr.write("heartbeat: " + self._status() + "\n")
            sys.stderr.flush()
