import pytest
from cold_atom_mot.atomic.angular_momentum import clebsch_gordan, wigner_3j


def test_known_clebsch_gordan_values():
    assert clebsch_gordan(2, 2, 1, 1, 3, 3) == pytest.approx(1.0)
    assert clebsch_gordan(1, 0, 1, 0, 1, 0) == pytest.approx(0.0, abs=1e-15)
    assert clebsch_gordan(1, 1, 1, -1, 0, 0) == pytest.approx(1 / 3**0.5)


def test_wigner_selection_rules():
    assert wigner_3j(1, 1, 1, 1, 1, -2) == 0.0
    assert wigner_3j(1, 1, 3, 0, 0, 0) == 0.0
