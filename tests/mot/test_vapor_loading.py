import numpy as np
import pytest

from cold_atom_mot.atomic.species import get_atomic_line
from cold_atom_mot.io.config import build_vapor_state, load_config, validate_config
from cold_atom_mot.simulation.capture import (
    CaptureCriterion,
    estimate_stratified_vapor_capture_rate,
    estimate_adaptive_vapor_capture_rate,
)
from cold_atom_mot.vacuum import (
    VaporState,
    flux_speed_cdf,
    loading_curve,
    one_sided_thermal_flux_m2_s,
    sample_flux_speeds,
    sample_spherical_inward_flux,
    flux_speed_pdf,
    gaussian_two_body_effective_volume,
    steady_state_population,
    wilson_interval,
)


class FreeParticle:
    atom = get_atomic_line("87Rb", "D2")
    gravity = np.zeros(3)

    def force(self, position, velocity, time=0):
        return np.zeros(3)


def test_flux_speed_sampler_matches_analytic_mean_and_cdf():
    atom = get_atomic_line("87Rb", "D2")
    rng = np.random.default_rng(12)
    speeds = sample_flux_speeds(300, atom.mass_kg, 100_000, rng)
    expected = 3 * np.sqrt(np.pi) / 4 * np.sqrt(2 * 1.380649e-23 * 300 / atom.mass_kg)
    assert speeds.mean() == pytest.approx(expected, rel=0.006)
    assert flux_speed_cdf(0, 300, atom.mass_kg) == 0
    assert flux_speed_cdf(2000, 300, atom.mass_kg) > 0.999


def test_spherical_flux_is_inward_and_reproducible():
    atom = get_atomic_line("87Rb", "D2")
    first = sample_spherical_inward_flux(0.01, 300, atom.mass_kg, 2000, seed=4)
    second = sample_spherical_inward_flux(0.01, 300, atom.mass_kg, 2000, seed=4)
    for a, b in zip(first, second):
        np.testing.assert_array_equal(a, b)
    positions, velocities, _ = first
    assert np.all(np.sum(positions * velocities, axis=1) < 0)
    np.testing.assert_allclose(np.linalg.norm(positions, axis=1), 0.01)


def test_stratified_capture_weights_and_loading_rate_are_consistent():
    vapor = VaporState(299, 301, 305, 1e-7, 2e-7, {"85Rb": 0.0, "87Rb": 1.0})
    criterion = CaptureCriterion(0.02, 1e6, 2e-5, 5e-6)
    estimate = estimate_stratified_vapor_capture_rate(
        FreeParticle(),
        vapor,
        criterion,
        capture_surface_radius_m=0.001,
        speed_bin_edges_m_s=[0, 20, 40],
        atoms_per_bin=4,
        max_step_s=1e-6,
        seed=8,
    )
    assert estimate.capture_probability == pytest.approx(estimate.sample_weights.sum())
    assert estimate.loading_rate_s == pytest.approx(
        estimate.incident_flux_s * estimate.capture_probability
    )
    assert estimate.omitted_high_speed_probability > 0.99


def test_loading_config_keeps_rb_and_background_pressures_independent():
    config = load_config("configs/rb_vapor_loading.yaml")
    vapor = build_vapor_state(config)
    assert vapor.rb_partial_pressure_pa != vapor.background_gas_pressure_pa
    modified = load_config("configs/rb_vapor_loading.yaml", validate=False)
    modified["vacuum"]["rb_partial_pressure_pa"] = 3e-8
    modified["vacuum"]["background_gas_pressure_pa"] = 7e-7
    validate_config(modified)
    direct = build_vapor_state(modified)
    assert direct.rb_partial_pressure_pa == 3e-8
    assert direct.background_gas_pressure_pa == 7e-7
    assert direct.rb_reservoir_temperature_k == 300
    assert direct.vapor_temperature_k == 300
    assert direct.background_temperature_k == 300


def test_loading_ode_covers_zero_loss_and_two_body_validation():
    time = np.linspace(0, 2, 5)
    np.testing.assert_allclose(loading_curve(time, 12, 0), 12 * time)
    with pytest.raises(ValueError, match="effective volume"):
        loading_curve(time, 12, 0.1, two_body_coefficient=1e-18)


def test_one_sided_flux_matches_n_vbar_over_four():
    atom = get_atomic_line("87Rb", "D2")
    density, temperature = 2e14, 310
    mean_speed = np.sqrt(8 * 1.380649e-23 * temperature / (np.pi * atom.mass_kg))
    assert one_sided_thermal_flux_m2_s(density, temperature, atom.mass_kg) == pytest.approx(
        density * mean_speed / 4
    )


def test_capture_dwell_uses_actual_adaptive_times():
    criterion = CaptureCriterion(0.02, 1e6, 2e-5, 5e-6)
    from cold_atom_mot.simulation.capture import evaluate_capture
    caught, times = evaluate_capture(
        FreeParticle(), [[-1e-3, 0, 0]], [[20, 0, 0]], criterion, max_step=1e-6
    )
    assert caught[0] and np.isfinite(times[0])

from cold_atom_mot.vacuum import background_collision_loss_rate_s


def test_background_collision_loss_requires_explicit_cross_section():
    atom = get_atomic_line("87Rb", "D2")
    zero = background_collision_loss_rate_s(1e-7, 300, atom.mass_kg, 4.65e-26, 0.0)
    finite = background_collision_loss_rate_s(1e-7, 300, atom.mass_kg, 4.65e-26, 1e-18)
    assert zero == 0 and finite > 0


def test_vapour_pressure_range_and_regression_values():
    from cold_atom_mot.vacuum import rubidium_vapor_pressure_pa
    assert rubidium_vapor_pressure_pa(300) == pytest.approx(10**(7.738-4215/300))
    assert rubidium_vapor_pressure_pa(400) == pytest.approx(10**(7.193-4040/400))
    with pytest.raises(ValueError, match="298.15"):
        rubidium_vapor_pressure_pa(280)
    with pytest.warns(RuntimeWarning):
        assert rubidium_vapor_pressure_pa(280,allow_extrapolation=True)>0


def test_zero_capture_wilson_interval_has_nonzero_upper_bound():
    low,high=wilson_interval(0,24)
    assert low == 0 and high > 0


def test_loading_arbitrary_initial_population_and_steady_state():
    t=np.linspace(0,3,20); n0=17.; rate=12.; gamma=.4
    expected=rate/gamma+(n0-rate/gamma)*np.exp(-gamma*t)
    np.testing.assert_allclose(loading_curve(t,rate,gamma,initial_population=n0),expected)
    np.testing.assert_allclose(loading_curve(t,0,gamma,initial_population=n0),n0*np.exp(-gamma*t))
    steady=steady_state_population(rate,gamma,two_body_coefficient=2e-18,effective_volume_m3=1e-9)
    coefficient=2e-18/1e-9
    assert rate-gamma*steady-coefficient*steady**2 == pytest.approx(0,abs=1e-10)


def test_gaussian_two_body_effective_volume():
    sx,sy,sz=1e-3,2e-3,3e-3
    assert gaussian_two_body_effective_volume(sx,sy,sz)==pytest.approx(8*np.pi**1.5*sx*sy*sz)


def test_capture_classification_stable_under_solver_refinement():
    from cold_atom_mot.simulation.capture import evaluate_capture
    criterion=CaptureCriterion(.02,1e6,2e-5,5e-6)
    coarse=evaluate_capture(FreeParticle(),[[-1e-3,0,0]],[[20,0,0]],criterion,
                            max_step=2e-6,rtol=1e-6,atol=1e-9)[0]
    fine=evaluate_capture(FreeParticle(),[[-1e-3,0,0]],[[20,0,0]],criterion,
                          max_step=5e-7,rtol=1e-9,atol=1e-12)[0]
    np.testing.assert_array_equal(coarse,fine)


def test_background_component_additivity_and_double_count_guard():
    atom=get_atomic_line("87Rb","D2")
    first=background_collision_loss_rate_s(1e-7,300,atom.mass_kg,4.65e-26,1e-18)
    second=background_collision_loss_rate_s(2e-7,300,atom.mass_kg,6.64e-27,2e-18)
    assert first+second > first
    config=load_config("configs/rb_vapor_loading.yaml",validate=False)
    config["loading"]["background_collision_model"]={
        "particle_mass_kg":4.65e-26,"effective_loss_cross_section_m2":1e-18}
    config["loading"]["background_gas_components"]=[{
        "species":"N2","partial_pressure_pa":1e-7,"particle_mass_kg":4.65e-26,
        "effective_loss_cross_section_m2":1e-18}]
    with pytest.raises(ValueError,match="double count"):
        validate_config(config)
    config["loading"]["background_collision_model"] = None
    with pytest.raises(ValueError, match="must be zero"):
        validate_config(config)


def test_adaptive_tail_reports_nonzero_zero_capture_bound():
    vapor=VaporState(300,300,300,1e-7,2e-7,{"85Rb":0.,"87Rb":1.})
    impossible=CaptureCriterion(1e-9,1e-9,2e-5,5e-6)
    estimate=estimate_adaptive_vapor_capture_rate(
        FreeParticle(),vapor,impossible,capture_surface_radius_m=.001,
        initial_speed_edges_m_s=[0,20],atoms_per_bin=4,maximum_speed_m_s=40,
        tail_relative_loading_tolerance=1e-3,max_step_s=1e-6,seed=4)
    assert estimate.capture_probability == 0
    assert estimate.omitted_capture_probability_upper > 0
    assert not estimate.tail_converged


def test_adaptive_tail_can_meet_statistical_stopping_rule():
    vapor=VaporState(300,300,300,1e-7,2e-7,{"85Rb":0.,"87Rb":1.})
    criterion=CaptureCriterion(.02,25.,2e-5,5e-6)
    estimate=estimate_adaptive_vapor_capture_rate(
        FreeParticle(),vapor,criterion,capture_surface_radius_m=.001,
        initial_speed_edges_m_s=[0,20,40],atoms_per_bin=4,maximum_speed_m_s=1600,
        tail_relative_loading_tolerance=1e-3,max_step_s=2e-6,seed=5)
    assert estimate.capture_probability > 0
    assert estimate.tail_converged
    assert estimate.last_simulated_speed_m_s <= 1600


def test_temperature_dependent_flux_weighting_changes_capture_integral():
    speed=np.linspace(0,80,2000); response=(speed<20).astype(float)
    atom=get_atomic_line("87Rb","D2")
    cold=np.trapezoid(response*flux_speed_pdf(speed,280,atom.mass_kg),speed)
    hot=np.trapezoid(response*flux_speed_pdf(speed,350,atom.mass_kg),speed)
    assert cold > hot > 0
