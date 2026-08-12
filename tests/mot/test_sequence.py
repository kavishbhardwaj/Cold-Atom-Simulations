from dataclasses import replace
import numpy as np
import pytest

from cold_atom_mot.io.config import build_effective_model, build_multilevel_model, load_config
from cold_atom_mot.simulation.sequence import (
    EddyCurrentResponse, ExperimentalSequence, Ramp, SequencedForce, Stage,
)


def stage(name="hold", duration=1.0, kind="linear"):
    return Stage(name, duration, Ramp(-1, -3, kind), Ramp(1, .2, kind),
                 Ramp(.1, .05, kind), Ramp(.1, 0, kind))


def test_ramp_endpoints_and_smooth_continuity():
    for kind in ("step", "linear", "smooth"):
        ramp = Ramp(2, 8, kind)
        assert ramp.value(0) == 2 and ramp.value(1) == 8
    smooth = Ramp(0, 1, "smooth")
    assert smooth.value(1e-6) < 1e-10
    assert 1-smooth.value(1-1e-6) < 1e-10


def test_stage_boundaries_and_zero_duration_are_deterministic():
    sequence = ExperimentalSequence((stage("zero", 0), stage("run", 2)))
    assert sequence.snapshot(0).stage == "run"
    assert sequence.snapshot(0).fraction == 0
    assert sequence.snapshot(2).fraction == 1
    with pytest.raises(ValueError, match="outside"):
        sequence.snapshot(2.1)


def test_eddy_current_and_coil_decay_are_analytical():
    response = EddyCurrentResponse(1, .2, .01, [3e-5, 0, 0], .02, [1e-6, 0, 0])
    assert response.gradient(1+.01) == pytest.approx(.2/np.e)
    assert response.field([0, 0, 0], 1+.02)[0] == pytest.approx(1e-6+3e-5/np.e)
    assert not response.is_time_independent


def test_sequence_controls_beams_and_static_sequence_reproduces_model():
    base = build_effective_model(load_config("configs/rb87_d2_mot.yaml"))
    beam = base.beams[0]
    hold = Stage("static", 1, Ramp(beam.detuning), Ramp(beam.power), Ramp(0),
                 Ramp(.1), cooling_frequency_offset=Ramp(beam.frequency_offset))
    sequence = ExperimentalSequence((hold,))
    wrapped = SequencedForce(base, sequence)
    point, velocity = np.array([1e-4, 0, 0]), np.array([.02, 0, 0])
    np.testing.assert_allclose(wrapped.force(point, velocity, .4),
                               base.force(point, velocity, .4))


def test_sequenced_trajectory_force_is_deterministically_reproducible():
    base = build_effective_model(load_config("configs/rb87_d2_mot.yaml"))
    sequence = ExperimentalSequence((stage(duration=.002),))
    force = SequencedForce(base, sequence)
    args = (np.zeros(3), np.array([.01, 0, 0]), .0007)
    np.testing.assert_array_equal(force.force(*args), force.force(*args))


def test_beam_step_linear_and_smooth_controls_reach_exact_endpoints():
    base = build_effective_model(load_config("configs/rb87_d2_mot.yaml"))
    controlled = Stage("pgc", .01, Ramp(-2, -6, "smooth"), Ramp(.01, .001, "linear"),
                       Ramp(.001, 0, "step"), Ramp(.1, 0, "smooth"),
                       cooling_coherence_groups=("x","x","y","y","z","z"),
                       cooling_phases=(0,1,2,3,4,5))
    sequence = ExperimentalSequence((controlled,))
    start, stop = sequence.apply_beams(base.beams, 0), sequence.apply_beams(base.beams, .01)
    assert start[0].power == .01 and stop[0].power == .001
    assert stop[0].detuning == -6 and stop[0].coherence_group == "x"
    assert stop[0].phase == 0


def test_internal_state_beam_families_consume_sequence_snapshot():
    rate = build_multilevel_model(load_config("configs/rb87_d2_multilevel.yaml"))
    sequence = ExperimentalSequence((stage(duration=1),))
    families = sequence.apply_beam_families(rate.beam_families, .5)
    assert len(families) == 12
    assert families[0].beam.power == pytest.approx(.6)
    assert families[6].beam.power == pytest.approx(.075)
