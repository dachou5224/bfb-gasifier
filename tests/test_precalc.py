# -*- coding: utf-8 -*-
"""T25b：预计算、初值与两状态释放。"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from hamel.hydrodynamics.bed import geldart_group
from hamel.properties.data import GAS_SPECIES, atoms, molar_mass_kg_mol
from hamel.properties.thermo import mixture_enthalpy_j
from hamel.pyrolysis.release import FuelTransport, release_distribution
from hamel.pyrolysis.yields import lignite_yields
from hamel.freeboard import account_freeboard, terminal_settling_velocity
from hamel.cases.config import SolverConfig
from hamel.solver.precalc import (
    PrecalcError, _adaptive_bed, bed_downflow, downward_mass_flow, freeboard_return_kg_s,
    initial_temperature_k, require_geldart_b, wall_coefficient,
)
from hamel.solver.precalc import _equilibrium, _gas_enthalpy, _load, _primary_temperature, _primary_vector
from hamel.solver.precalc import _transport as pack_transport
from hamel.solver.residual import wall_area_m2

_CASE = Path("data/cases/vtt_pressurized.json")
_DOC = json.loads(_CASE.read_text(encoding="utf-8"))
_FUEL = _DOC["fuels"]["rhein_braunkohle"]
_ASH = float(_FUEL["ash_mass_fraction_wf"]["value"])
_MOISTURE = float(_DOC["cases"]["1"]["moisture_mass_fraction"]["value"])
_BED_K = 1092.15
_PRESSURE = 5.0e5


def test_geldart_b_and_d():
    """ADR 0024 第 4 项。d' 用 μm，密度差用 g/cm³。"""
    assert geldart_group(300.0e-6, 1500.0, 0.0) == "B"
    assert require_geldart_b(300.0e-6, 1500.0, 0.0) == "B"
    assert geldart_group(2.0e-3, 1500.0, 0.0) == "D"
    with pytest.raises(PrecalcError, match="不是 B"):
        require_geldart_b(2.0e-3, 1500.0, 0.0)


def test_downward_flow_reports_the_bottom_residual():
    """ADR 0009 第 2 项。闭合时残差为 0；最底向下流保持 0。"""
    up = np.array([1.0, 1.0, 1.0])
    feed = np.array([0.0, 0.0, 0.2])
    down, residual = downward_mass_flow(up, feed, np.zeros(3), np.array([0.0, 0.0, 0.2]), np.zeros(3), 1.0)
    assert down == pytest.approx([1.0, 1.0, 0.0])
    assert residual == pytest.approx(0.0)
    up = np.array([1.0, 1.2, 1.5])
    down, residual = downward_mass_flow(
        up, np.array([0.4, 0.0, 0.0]), np.zeros(3), np.array([0.0, 0.0, 0.3]),
        np.array([0.1, 0.1, 0.05]), 0.8,
    )
    assert down[-1] == 0.0
    assert residual == pytest.approx(0.15)


def test_wall_coefficient_returns_the_measured_loss():
    """ADR 0029。面积用残差模块的 wall_area_m2。"""
    bed_t = np.full(8, 1100.0)
    free_t = np.full(8, 1000.0)
    loss = 21700.0
    coefficient = wall_coefficient(bed_t, free_t, 0.15, 0.2875, 0.15, 0.25, loss, adiabatic=False)
    total = 0.0
    for index, temperature in enumerate(bed_t):
        total += wall_area_m2("bed", index, 8, 0.15, 0.15, 0.25) * (temperature - 298.15)
    for index, temperature in enumerate(free_t):
        total += wall_area_m2("freeboard", index, 8, 0.2875, 0.15, 0.25) * (temperature - 298.15)
    assert coefficient * total == pytest.approx(loss)
    assert wall_coefficient(bed_t, free_t, 0.15, 0.2875, 0.15, 0.25, loss, adiabatic=True) == 0.0
    with pytest.raises(PrecalcError):
        wall_coefficient(np.full(2, 298.15), np.full(2, 200.0), 0.15, 0.2875, 0.15, 0.25, loss, adiabatic=False)


def test_feed_equilibrium_conserves_elements_and_enthalpy():
    """ADR 0017 第 9 项。产物元素等于进料元素，产物焓等于进料焓减 21.7 kW。"""
    case = _load(_CASE, "1")
    result = _equilibrium(case, adiabatic=False)
    expected = {}
    for name, amount in result.amounts.items():
        for element, count in atoms(name).items():
            expected[element] = expected.get(element, 0.0) + count * amount
    feed = _feed_elements(case)
    for element, amount in feed.items():
        assert expected.get(element, 0.0) == pytest.approx(amount, rel=1.0e-6, abs=1.0e-8)
    assert result.product_enthalpy == pytest.approx(result.feed_enthalpy - case.heat_loss_w, rel=1.0e-6, abs=1.0)
    assert mixture_enthalpy_j(result.temperature, result.amounts) == pytest.approx(result.product_enthalpy)
    hot = _equilibrium(case, adiabatic=True)
    assert hot.product_enthalpy == pytest.approx(hot.feed_enthalpy, rel=1.0e-6, abs=1.0)
    assert hot.temperature > result.temperature


def test_release_with_capture_partial_conversion_and_absent_free_state():
    """ADR 0032。非零捕获、部分未转化、目标格没有自由状态。"""
    yields = _yields()
    captured = _transport(False)
    partial = _run(np.array([2.0]), np.array([3.0]), captured)
    assert partial.undried_kg.sum() > 0.0
    assert partial.dried_kg.sum() > 0.0
    assert partial.char_hover_kg_s is not None and partial.char_hover_kg_s.sum() > 0.0
    assert partial.fuel_up_kg_s[0].sum() > 0.0
    absent = _run(np.array([1.0]), np.array([1.0]), _transport(True))
    assert absent.undried_free_kg[1] == pytest.approx(0.0)
    assert absent.undried_hover_kg[1] + absent.dried_hover_kg[1] > 0.0
    assert absent.char_free_kg_s[1] == pytest.approx(0.0)
    assert absent.char_hover_kg_s[1] > 0.0


def test_hover_drying_stays_in_the_hover_state():
    """ADR 0032。悬停中干燥得到的 P 不回到自由状态。"""
    yields = _yields()
    transport = FuelTransport(
        1, np.array([[True]]), np.array([[False]]), np.array([[0.5]]), np.zeros((1, 1)),
        np.zeros((1, 1)), np.zeros((1, 1)), np.array([0.0]), np.array([1.0]),
    )
    result = _run_cells(np.array([1.0]), np.array([1.0]), transport, feed=1.0, up=1.0)
    assert result.undried_free_kg[0, 0] == pytest.approx(0.0)
    assert result.undried_hover_kg[0, 0] == pytest.approx(0.5)
    assert result.dried_free_kg[0, 0] == pytest.approx(0.0)
    assert result.dried_hover_kg[0, 0] == pytest.approx(0.75)
    assert result.char_free_kg_s[0, 0] == pytest.approx(0.0)
    assert result.fuel_up_kg_s[0, 0] == pytest.approx(0.0)
    waf = (1.0 - _MOISTURE) * (1.0 - _ASH)
    assert result.char_hover_kg_s[0, 0] == pytest.approx(waf * yields.x_char * 0.75)
    instant = _run_cells(np.array([1.0]), np.array([0.0]), transport, feed=1.0, up=1.0)
    assert instant.char_free_kg_s[0, 0] == pytest.approx(0.0)
    assert instant.char_hover_kg_s[0, 0] == pytest.approx(waf * yields.x_char * 0.5)


def test_hover_exists_when_theta_is_zero():
    """捕获区决定悬停。θ = 0 只表示当地轨迹不是起步即被捕获。"""
    edges = np.array([0.0, 1.0, 1.25, 2.5])
    rho = np.full(3, 0.4)
    mu = np.full(3, 1.5e-5)
    diameter, density, phi = 1.0e-3, 1500.0, 0.8
    settling = terminal_settling_velocity(0.4, 1.5e-5, diameter, density, phi)

    def gas(_cell, height, settling=settling):
        return settling + 8.0 * (1.01 - float(height))

    book = account_freeboard(
        edges, gas, rho, mu, diameter, density, phi, np.array([0.1]), np.array([1.0]), 0.0,
        SolverConfig(freeboard_delta_h_m=1.0e-3, freeboard_max_step_s=1.0e-3, freeboard_time_limit_s=20.0),
    )
    assert book.capture_cells == (1,)
    assert all(cell.theta == 0 and cell.free_state == "存在" for cell in book.cells)
    assert all(cell.mean_slip_m_s is not None and cell.mean_slip_m_s > 0.0 for cell in book.cells)
    _packed, theta, hover, absent = pack_transport((book,))
    assert absent.sum() == 0
    assert theta.sum() == 0.0
    assert hover[1, 0]
    assert not hover[0, 0] and not hover[2, 0]


def test_return_uses_the_freeboard_inventory_and_the_drain_enters_the_residual():
    """ADR 0016 第 10 条。同一床顶、不同自由存量，回落不同；排灰进入底格残差。"""
    coeff = np.array([0.25, 0.25])
    empty = freeboard_return_kg_s(coeff, np.zeros(2), np.zeros(2), np.zeros(2))
    filled = freeboard_return_kg_s(coeff, np.zeros(2), np.array([1.0, 1.0]), np.array([0.5, 0.5]), 0.1)
    assert empty == 0.0
    assert filled == pytest.approx(0.85)
    up = np.array([1.0, 1.0, 1.0])
    _down, closed = bed_downflow(up, 0.0, np.zeros(3), 1.0, 0.0)
    _down, drained = bed_downflow(up, 0.0, np.zeros(3), 1.0, 0.3)
    _down, reacted = bed_downflow(up, 0.0, np.array([0.0, 0.0, -0.2]), 1.0, 0.0)
    assert closed == pytest.approx(0.0)
    assert drained == pytest.approx(-0.3)
    assert reacted == pytest.approx(-0.2)


def test_primary_enthalpy_keeps_the_purge_temperature():
    """ADR 0014 C2。冻结的一次风温度复现分温度的混合焓。"""
    case = _load(_CASE, "1")
    temperature = _primary_temperature(case)
    vector = _primary_vector(case)
    flows = {name: float(amount) for name, amount in zip(GAS_SPECIES, vector) if amount}
    air = float(case.physical.primary_air_amount_flow)
    purge = float(case.physical.nitrogen_purge_volume_flow)
    steam = float(case.physical.steam_mass_flow) / molar_mass_kg_mol("H2O")
    inlet = float(case.physical.inlet_temperature)
    target = (
        _gas_enthalpy({"O2": 0.21 * air, "N2": 0.79 * air}, inlet)
        + _gas_enthalpy({"H2O": steam}, inlet)
        + _gas_enthalpy({"N2": purge}, case.purge_temperature_k)
    )
    assert _gas_enthalpy(flows, temperature) == pytest.approx(target, rel=1.0e-9)
    assert temperature != pytest.approx(inlet)


def _yields():
    return lignite_yields(
        _PRESSURE, _BED_K, float(_FUEL["carbon_mass_fraction_wf"]["value"]),
        float(_FUEL["hydrogen_mass_fraction_wf"]["value"]), float(_FUEL["oxygen_mass_fraction_wf"]["value"]),
        float(_FUEL["nitrogen_mass_fraction_wf"]["value"]), float(_FUEL["sulfur_mass_fraction_wf"]["value"]), _ASH,
    )


def _transport(bottom_absent: bool) -> FuelTransport:
    hover = np.array([[False], [True]])
    absent = np.array([[False], [bottom_absent]])
    if bottom_absent:
        hover[1, 0] = True
    a_free = np.array([[0.4], [0.0 if bottom_absent else 0.2]])
    a_hover = np.array([[0.0], [0.0 if bottom_absent else 0.3]])
    b_free = np.array([[0.2], [0.0]])
    b_hover = np.array([[0.0], [0.0]])
    to_free = np.array([0.0 if bottom_absent else 0.25])
    to_hover = np.array([1.0 if bottom_absent else 0.75])
    return FuelTransport(1, hover, absent, a_free, a_hover, b_free, b_hover, to_free, to_hover)


def _run_cells(tau_dry, tau_pyr, transport, feed=0.01, up=0.5):
    upward = np.array([[0.0], [up]])
    downward = np.zeros((2, 1))
    leaving = np.zeros((2, 1))
    return release_distribution(
        np.array([feed]), 0.2, [1.0, 0.0], [2.0, 1.0], upward, downward, leaving,
        tau_dry, tau_pyr, _MOISTURE, _ASH, _yields(), _BED_K, _PRESSURE, transport=transport,
    )


def _run(tau_dry, tau_pyr, transport):
    upward = np.array([[0.0], [0.0], [0.5]])
    downward = np.zeros((3, 1))
    leaving = np.array([[0.0], [0.0], [0.2]])
    return release_distribution(
        np.array([0.01]), 0.05, [2.0, 1.0, 0.0], [3.0, 2.0, 1.0], upward, downward, leaving,
        tau_dry, tau_pyr, _MOISTURE, _ASH, _yields(), _BED_K, _PRESSURE, transport=transport,
    )


def _feed_elements(case) -> dict[str, float]:
    from hamel.properties.data import atomic_weight
    from hamel.solver.precalc import _elements
    physical = case.physical
    air = float(physical.primary_air_amount_flow)
    steam = float(physical.steam_mass_flow) / molar_mass_kg_mol("H2O")
    purge = float(physical.nitrogen_purge_volume_flow)
    primary = {name: 0.0 for name in GAS_SPECIES}
    primary["O2"], primary["N2"], primary["H2O"] = 0.21 * air, 0.79 * air + purge, steam
    secondary_n = float(physical.secondary_air_amount_flow)
    secondary = {name: 0.0 for name in GAS_SPECIES}
    secondary["O2"], secondary["N2"] = 0.21 * secondary_n, 0.79 * secondary_n
    return _elements(case, primary, secondary)


def test_yield_range_follows_only_the_substituted_temperature():
    """产率式只约束实际代入的温度。自由空间气温不另成整炉约束。"""
    from hamel.solver.precalc import _load, _yield_at

    case = _load(_CASE, "1")
    notes = []
    legal = _yield_at(case, 1100.0, None, None, notes, 0, True)
    assert notes == []
    assert legal.x_char > 0.0
    recorded = []
    _yield_at(case, 1000.0, None, None, recorded, 0, False)
    assert recorded == [(0, "产率式温度或压力", 1000.0, None, None)]
    with pytest.raises(PrecalcError, match="产率式"):
        _yield_at(case, 1000.0, None, None, [], 0, True)


def test_initial_temperature_is_projected_onto_the_yield_interval():
    """求解手段：ADR 0035。报告里的 T_ad 不改，初值用投影后的温度。"""
    assert initial_temperature_k(971.5) == 1073.0
    assert initial_temperature_k(1095.7) == pytest.approx(1095.7)
    assert initial_temperature_k(1400.0) == 1273.0


def _run_bounded(code: str, seconds: float) -> subprocess.CompletedProcess:
    """墙钟到点就杀掉整个进程组，并记 TIMEOUT。CHARTER §7。"""
    proc = subprocess.Popen(
        [sys.executable, "-c", code],
        cwd=str(Path(__file__).resolve().parents[1]),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        out, err = proc.communicate(timeout=seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        pytest.fail("TIMEOUT")
    return subprocess.CompletedProcess(proc.args, proc.returncode, out, err)


def test_precalculation_timeout_is_recorded():
    with pytest.raises(pytest.fail.Exception, match="TIMEOUT"):
        _run_bounded("import time; time.sleep(30)", 0.2)


_CASE1_SMOKE = r"""
import json, math
from hamel.solver.precalc import initial_precalculation
result = initial_precalculation("data/cases/vtt_pressurized.json", "1")
_previous, _finer, rho_rel, masses = result.substep_history[-1]
print(json.dumps({
    "geldart": result.geldart,
    "t_init_k": result.t_init_k,
    "bed_intervals": result.bed_intervals,
    "rho_rel": rho_rel,
    "max_mass": max(masses),
    "slip_ok": bool(result.slip_ratios) and all(item[4] > 0.0 and math.isfinite(item[5]) for item in result.slip_ratios),
}))
"""


@pytest.mark.slow
def test_case1_first_iteration_finishes():
    """工况 1 第一次外层。默认集不跑；墙钟 900 s 到点记 TIMEOUT。"""
    completed = _run_bounded(_CASE1_SMOKE, 900)
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout.strip().splitlines()[-1])
    assert payload["geldart"] == "B" and payload["t_init_k"] == 1073.0
    assert payload["bed_intervals"] >= 16
    assert payload["rho_rel"] < 1.0e-6 and payload["max_mass"] < 1.0e-6
    assert payload["slip_ok"] is True


def test_bed_substeps_stop_at_the_cap_without_splicing_cells():
    """求解器参数：ADR 0036。到上限仍不过门就失败，并带上各次整床相对变化。"""
    from dataclasses import replace

    from hamel.solver.precalc import _bed_gas, _blank_state, _class_fractions, _diameters, _equilibrium, _flows, _load

    case = _load(_CASE, "1")
    case = replace(case, solver=replace(case.solver, bed_integration_max_intervals_per_cell=16))
    equilibrium = _equilibrium(case, False)
    temperature = initial_temperature_k(equilibrium.temperature)
    fractions = _class_fractions(None, case)
    d_bar, _sv = _diameters(case.diameters_m, fractions)
    flow = _flows(case, _blank_state(case, temperature), equilibrium, True, temperature)
    with pytest.raises(PrecalcError, match="上限") as caught:
        _adaptive_bed(case, flow, _bed_gas(flow), d_bar)
    text = str(caught.value)
    assert "16" in text and text.count("(") > 1


def test_density_search_recovers_the_requested_bed_mass():
    """ADR 0014 A4。先由 bed_profile 得存量，再反求 ρ_P，存量必须回到该值。"""
    from dataclasses import replace

    from hamel.hydrodynamics.bed import BedConditions, bed_profile, minimum_fluidization_velocity, particle_density_from_bed_mass
    from hamel.solver.precalc import PHI_S, SEARCH

    edges = np.linspace(0.0, 0.6, 5)
    conditions = BedConditions(1.0e-3, 0.46, 900.0, 1.2, 1.8e-5, 1.5e-5, 101300.0, 2.0e-5, 0.15, "technical", "B", True)
    minimum = minimum_fluidization_velocity(0.46, conditions.nu_g, PHI_S, conditions.d_p, 900.0, conditions.rho_g)
    velocities = np.full(4, minimum + 0.4)
    target = float(bed_profile(edges, velocities, minimum, 0.02, 40, 0.0, 4, conditions, SEARCH).solids_mass.sum())
    density = particle_density_from_bed_mass(target, edges, velocities, 0.02, 40, 0.0, 4, PHI_S, conditions, SEARCH)
    checked = minimum_fluidization_velocity(0.46, conditions.nu_g, PHI_S, conditions.d_p, density, conditions.rho_g)
    found = bed_profile(edges, velocities, checked, 0.02, 40, 0.0, 4, replace(conditions, rho_p=density), SEARCH)
    assert found.solids_mass.sum() == pytest.approx(target, rel=1.0e-6)
