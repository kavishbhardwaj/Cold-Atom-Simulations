import pytest
from sympy import S
from sympy.physics.wigner import clebsch_gordan, wigner_3j


def test_known_clebsch_gordan_values():
    assert float(clebsch_gordan(2, 1, 3, 2, 1, 3)) == pytest.approx(1.0)
    assert float(clebsch_gordan(1, 1, 1, 0, 0, 0)) == pytest.approx(0.0, abs=1e-15)
    assert abs(float(clebsch_gordan(1, 1, 0, 1, -1, 0))) == pytest.approx(1 / 3**0.5)


def test_wigner_selection_rules():
    assert wigner_3j(1, 1, 1, 1, 1, -2) == 0
    assert wigner_3j(1, 1, 3, 0, 0, 0) == 0
