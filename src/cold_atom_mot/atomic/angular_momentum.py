"""Small-integer angular-momentum coefficients used by the 87Rb model.

The implementation follows the Racah factorial expression and is intentionally
limited to integer F,m values needed for 87Rb hyperfine manifolds.  It avoids an
opaque table of Zeeman strengths while keeping the dependency surface small.
"""

from math import factorial, sqrt


def _fact(value: int) -> int:
    if value < 0:
        return 0
    return factorial(value)


def wigner_3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """Return a Wigner 3-j symbol for integer angular momenta."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    triangle = _fact(j1 + j2 - j3) * _fact(j1 - j2 + j3) * _fact(-j1 + j2 + j3)
    triangle /= _fact(j1 + j2 + j3 + 1)
    magnetic = (
        _fact(j1 + m1) * _fact(j1 - m1)
        * _fact(j2 + m2) * _fact(j2 - m2)
        * _fact(j3 + m3) * _fact(j3 - m3)
    )
    prefactor = (-1) ** (j1 - j2 - m3) * sqrt(triangle * magnetic)
    lower = max(0, j2 - j3 - m1, j1 - j3 + m2)
    upper = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    total = 0.0
    for z in range(lower, upper + 1):
        denominator = (
            _fact(z) * _fact(j1 + j2 - j3 - z)
            * _fact(j1 - m1 - z) * _fact(j2 + m2 - z)
            * _fact(j3 - j2 + m1 + z) * _fact(j3 - j1 - m2 + z)
        )
        total += (-1) ** z / denominator
    return prefactor * total


def clebsch_gordan(j1: int, m1: int, j2: int, m2: int, total_j: int, total_m: int) -> float:
    """Return <j1,m1;j2,m2|J,M> for integer angular momenta."""
    if m1 + m2 != total_m:
        return 0.0
    return (
        (-1) ** (j1 - j2 + total_m)
        * sqrt(2 * total_j + 1)
        * wigner_3j(j1, j2, total_j, m1, m2, -total_m)
    )
