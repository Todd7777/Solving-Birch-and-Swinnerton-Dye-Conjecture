"""
bsd — computational verification toolkit for the Birch and Swinnerton-Dyer
conjecture over Q, built on PARI/GP (via cypari2).

Nothing in this package proves the conjecture.  It computes every quantity
that appears in the BSD formula for explicit elliptic curves, checks the
formula numerically, records what is *rigorously* established for each
curve by the theorems in the literature, and keeps the two clearly apart.
"""

__version__ = "0.1.0"

from .verify import verify_curve, BSDRecord  # noqa: F401
