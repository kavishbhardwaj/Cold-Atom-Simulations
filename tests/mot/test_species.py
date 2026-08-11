import pytest
from cold_atom_mot.atomic.species import ATOMIC_LINES, get_atomic_line


def test_registry_distinguishes_both_isotopes_and_lines():
    assert set(ATOMIC_LINES) == {('87Rb', 'D1'), ('87Rb', 'D2'), ('85Rb', 'D1'), ('85Rb', 'D2')}
    assert get_atomic_line('87Rb', 'D2').wavelength_m < get_atomic_line('87Rb', 'D1').wavelength_m
    assert get_atomic_line('87Rb', 'D2').model_support.startswith('Level A')
    assert 'not yet implemented' in get_atomic_line('85Rb', 'D2').model_support


def test_registry_rejects_unknown_line():
    with pytest.raises(ValueError, match='choose one of'):
        get_atomic_line('86Rb', 'D2')
