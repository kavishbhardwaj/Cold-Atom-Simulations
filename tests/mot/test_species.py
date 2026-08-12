import numpy as np
import pytest
from cold_atom_mot.atomic.species import ATOMIC_LINES, get_atomic_line


def test_registry_distinguishes_both_isotopes_and_lines():
    assert set(ATOMIC_LINES) == {('87Rb', 'D1'), ('87Rb', 'D2'), ('85Rb', 'D1'), ('85Rb', 'D2')}
    assert get_atomic_line('87Rb', 'D2').wavelength_m < get_atomic_line('87Rb', 'D1').wavelength_m
    assert 'multilevel rate-equation MOT' in get_atomic_line('87Rb', 'D2').model_support
    assert 'multilevel rate-equation MOT' in get_atomic_line('85Rb', 'D2').model_support
    rb87_d2 = get_atomic_line('87Rb', 'D2')
    rb85_d2 = get_atomic_line('85Rb', 'D2')
    assert rb87_d2.recoil_velocity_m_s > 0
    assert rb87_d2.recoil_temperature_k > 0
    assert rb87_d2.doppler_temperature_k > rb87_d2.recoil_temperature_k
    assert rb85_d2.natural_abundance > rb87_d2.natural_abundance


def test_registry_rejects_unknown_line():
    with pytest.raises(ValueError, match='choose one of'):
        get_atomic_line('86Rb', 'D2')

from cold_atom_mot.atomic.species import build_atomic_basis

@pytest.mark.parametrize("isotope,line,size", [("87Rb","D2",24),("87Rb","D1",16),("85Rb","D2",36),("85Rb","D1",24)])
def test_generated_basis_dimensions_and_branching(isotope,line,size):
    basis=build_atomic_basis(isotope,line)
    assert basis.state_count == size
    np.testing.assert_allclose(basis.spontaneous_branching.sum(axis=1),1.0,atol=1e-14)
    assert all(t.q in (-1,0,1) for t in basis.transitions)
    assert all(abs(basis.excited[t.excited_index].F-basis.ground[t.ground_index].F)<=1 for t in basis.transitions)


def test_d2_cycling_transitions_and_d1_capability_boundary():
    assert get_atomic_line("87Rb","D2").cooling_transition == (2,3)
    assert get_atomic_line("85Rb","D2").cooling_transition == (3,4)
    assert get_atomic_line("87Rb","D1").rate_equation_mot is False
    for isotope,fg,fe in (("87Rb",2,3),("85Rb",3,4)):
        basis=build_atomic_basis(isotope,"D2")
        stretched=[t for t in basis.transitions if basis.ground[t.ground_index].F==fg and basis.ground[t.ground_index].m==fg and basis.excited[t.excited_index].F==fe and basis.excited[t.excited_index].m==fe]
        assert len(stretched)==1 and stretched[0].strength == pytest.approx(1.0)
