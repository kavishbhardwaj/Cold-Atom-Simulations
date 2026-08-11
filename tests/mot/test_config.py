import copy
import pytest
from cold_atom_mot.io.config import load_config, validate_config


def test_config_rejects_negative_power():
    config = load_config('configs/rb87_standard_mot.yaml')
    config = copy.deepcopy(config); config['laser']['power_per_beam_w'] = -1
    with pytest.raises(ValueError, match="power"):
        validate_config(config)


def test_standard_config_loads():
    config = load_config('configs/rb87_standard_mot.yaml')
    assert config['monte_carlo']['seed'] == 20260811


def test_phase2_config_rejects_negative_repump_power():
    config = load_config('configs/rb87_phase2_rate_equation.yaml')
    config['repump']['power_per_beam_w'] = -1e-3
    with pytest.raises(ValueError, match="repump"):
        validate_config(config)
