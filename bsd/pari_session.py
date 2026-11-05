"""A single PARI/GP session shared by the package, plus a timeout helper.

PARI calls are C calls; a Python-level signal cannot interrupt them.  We use
cysignals' alarm, which cypari2 honours inside every sig_on()-protected call.
"""
from __future__ import annotations

import contextlib
import time

import cypari2
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt

_PARI: cypari2.Pari | None = None
DEFAULT_STACK = 2**30          # 1 GiB PARI stack
DEFAULT_DIGITS = 38            # ~128 bits of real precision


def get_pari(stack: int = DEFAULT_STACK, digits: int = DEFAULT_DIGITS) -> cypari2.Pari:
    """Return the process-wide PARI instance (created on first use)."""
    global _PARI
    if _PARI is None:
        _PARI = cypari2.Pari()
        try:
            _PARI.allocatemem(stack, silent=True)
        except TypeError:  # older cypari2
            _PARI.allocatemem(stack)
        _PARI.set_real_precision(digits)
    return _PARI


class PariTimeout(Exception):
    """Raised when a PARI computation exceeds its time budget."""


@contextlib.contextmanager
def time_limit(seconds: float | None):
    """Context manager: abort the enclosed PARI computation after `seconds`."""
    if seconds is None or seconds <= 0:
        yield
        return
    t0 = time.time()
    try:
        alarm(seconds)
        yield
    except AlarmInterrupt as exc:  # pragma: no cover - timing dependent
        raise PariTimeout(f"PARI computation exceeded {seconds}s (ran {time.time()-t0:.1f}s)") from exc
    finally:
        cancel_alarm()


def gp(expr: str):
    """Evaluate a GP expression in the shared session."""
    return get_pari()(expr)
