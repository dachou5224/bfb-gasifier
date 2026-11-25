"""T31 反应路径账：默认关闭、带符号 q、φ／焦油分项与 C.0 三列。"""
from __future__ import annotations

import numpy as np
import pytest

from hamel.properties.data import GAS_SPECIES, molar_mass_kg_mol
from hamel.solver.residual import _apply, evaluate_residual
from hamel.solver.vector import pack_state, unpack_state
from test_assembly import make_frozen, make_state
from tools.ledger_account import account, aggregated, gas_residual_rows, registry

_INDEX = {name: index for index, name in enumerate(GAS_SPECIES)}
_MOLAR = molar_mass_kg_mol("C_gr")


def _one(reaction, extent, **kwargs):
    """一次 _apply，返回它的账目记录（q 已去掉结构零）。"""
    rows = []
    gas = np.zeros(len(GAS_SPECIES))
    _apply(gas, reaction, extent, _MOLAR, sink=lambda r, p, t, q: rows.append({"reaction": r, "phi": p, "tar": t, "q": q}), **kwargs)
    return rows[0]


def _record(zone, cell, phase, record):
    return {"zone": zone, "cell": cell, "phase": phase, **record}


def test_closed_sink_keeps_the_residual_bitwise_identical():
    frozen = make_frozen()
    values = pack_state(make_state(frozen), frozen)
    off = evaluate_residual(values, frozen)
    again = evaluate_residual(values, frozen)
    on = evaluate_residual(values, frozen, ledger=True)
    assert np.array_equal(off.values, again.values)
    assert np.array_equal(off.values, on.values)
    assert off.ledger is None and again.ledger is None
    assert len(on.ledger) > 0


def test_ledger_reproduces_the_assembled_reaction_sources():
    frozen = make_frozen()
    values = pack_state(make_state(frozen), frozen)
    result = evaluate_residual(values, frozen, ledger=True)
    sums = np.zeros_like(result.reaction_bed)
    free = np.zeros_like(result.reaction_free)
    char_bed = np.zeros(frozen.n_bed)
    char_free = np.zeros(frozen.n_freeboard)
    for row in result.ledger:
        cell, species = int(row["cell"]), row["q"]
        for name, q in species.items():
            if name == "C":
                continue
            if row["zone"] == "bed":
                sums[cell, 0 if row["phase"] == "suspension" else 1, _INDEX[name]] += q
            else:
                free[cell, _INDEX[name]] += q
    for row in result.ledger:
        if "C" in row["q"]:
            target = char_bed if row["zone"] == "bed" else char_free
            target[int(row["cell"])] += row["q"]["C"] * _MOLAR
    assert np.allclose(sums, result.reaction_bed, rtol=0.0, atol=1.0e-15)
    assert np.allclose(free, result.reaction_free, rtol=0.0, atol=1.0e-15)
    assert np.allclose(char_bed, -result.char_consumption_bed_kg_s.sum(axis=1), rtol=0.0, atol=1.0e-15)
    assert np.allclose(char_free, -(result.char_consumption_free_kg_s + result.char_consumption_hover_kg_s).sum(axis=1), rtol=0.0, atol=1.0e-15)


def test_r1_records_keep_each_phi_separately():
    first = _one("R1", 2.0, phi=1.0)
    second = _one("R1", 2.0, phi=2.0)
    assert first["q"]["CO"] == 0.0 and first["q"]["CO2"] == 2.0
    assert second["q"]["CO"] == 2.0 and second["q"]["CO2"] == 0.0
    assert second["q"]["O2"] == -1.0  # −1/φ
    rows = aggregated([_record("bed", 0, "suspension", first), _record("bed", 0, "suspension", second)])
    assert len(rows) == 2 and {row["phi"] for row in rows} == {1.0, 2.0}


def test_two_tar_species_are_separate_records():
    benzene = _one("R10", 1.0, tar="C6H6")
    naphthalene = _one("R10", 1.0, tar="C10H8")
    assert benzene["q"]["O2"] == -3.0 and benzene["q"]["CO"] == 6.0
    assert naphthalene["q"]["O2"] == -5.0 and naphthalene["q"]["CO"] == 10.0
    assert benzene["tar"] == "C6H6" and naphthalene["tar"] == "C10H8"


def test_freeboard_with_hover_is_in_the_ledger():
    hover = np.zeros((2, 2), dtype=bool)
    hover[0, 1] = True
    frozen = make_frozen(n_free=2, n_class=2, hover=hover)
    values = pack_state(make_state(frozen), frozen)
    result = evaluate_residual(values, frozen, ledger=True)
    phases = {(row["zone"], row["phase"]) for row in result.ledger}
    assert ("freeboard", "free") in phases and ("freeboard", "hover") in phases
    free_rows = sum(1 for row in result.ledger if row["zone"] == "freeboard" and row["phase"] == "hover")
    assert free_rows > 0


def test_negative_r8_generates_water_and_consumes_hydrogen():
    reverse = _record("bed", 0, "suspension", _one("R8", -0.7))
    book = account([reverse])
    assert book["species"]["H2O"]["generate"] == pytest.approx(0.7)
    assert book["species"]["H2O"]["consume"] == 0.0
    assert book["species"]["H2"]["consume"] == pytest.approx(0.7)
    assert book["species"]["H2"]["generate"] == 0.0
    forward = _record("bed", 0, "suspension", _one("R8", 0.5))
    other = account([forward])
    assert other["species"]["H2O"]["consume"] == pytest.approx(0.5)


def test_r11_water_uses_the_actual_tar_coefficient():
    benzene = _record("bed", 0, "suspension", _one("R11", 2.0, tar="C6H6"))
    naphthalene = _record("bed", 0, "suspension", _one("R11", 1.5, tar="C10H8"))
    book = account([benzene, naphthalene])
    assert book["species"]["H2O"]["consume"] == pytest.approx(2.0 * 4.5 + 1.5 * 23.0 / 3.0)
    # 只按 −1·R11 会得到 3.5，必须不是这个数。
    assert book["species"]["H2O"]["consume"] != pytest.approx(2.0 + 1.5)


def test_char_q_is_negative_for_positive_progress_and_needs_no_extra_sign():
    for reaction in ("R1", "R2", "R3", "R4"):
        gas = np.zeros(len(GAS_SPECIES))
        rows = []
        returned = _apply(gas, reaction, 2.0, _MOLAR, phi=1.2 if reaction == "R1" else None,
                          sink=lambda r, p, t, q: rows.append(q))
        assert returned < 0.0
        assert rows[0]["C"] == pytest.approx(-2.0)
        assert rows[0]["C"] * _MOLAR == pytest.approx(returned)


def test_opposite_r8_in_two_cells_keeps_the_turnover():
    forward = _record("bed", 0, "suspension", _one("R8", 1.0))
    reverse = _record("bed", 3, "free", _one("R8", -0.8))
    water = account([forward, reverse])["species"]["H2O"]
    assert water["consume"] == pytest.approx(1.0)
    assert water["generate"] == pytest.approx(0.8)
    assert water["net"] == pytest.approx(-0.2)


def test_registry_reports_each_path_contribution_not_the_total():
    """占比表要给出该路径自己的贡献；分母为零时记“分母可忽略”，不退回 0 占比。"""
    r2 = _record("bed", 0, "suspension", _one("R2", 1.5))
    r8 = _record("bed", 1, "suspension", _one("R8", 0.5))
    reg = registry([r2, r8])
    assert reg["water_consume"] == pytest.approx(2.0)
    assert reg["r2_water_contribution"] == pytest.approx(1.5)
    assert reg["r2_water_share"] == pytest.approx(0.75)
    assert reg["tar_water_contribution"] == 0.0 and reg["tar_water_share"] == 0.0
    only_generate = registry([_record("bed", 0, "bubble", _one("R6", 1.0))])
    assert only_generate["water_consume"] == 0.0
    assert only_generate["r2_water_share"] is None
    assert only_generate["tar_water_share"] is None


def test_gas_residual_rows_match_the_assembled_residual_with_exchange_sign():
    """C.3／D.2 共用的逐格逐相残差必须等于装配值，交换项只带一次相符号。"""
    frozen = make_frozen()
    state = make_state(frozen, temperature=900.0)
    state.bed_gas_mol_s[:, 1, _INDEX["N2"]] += 0.01  # 两相组成不同，交换项才非零
    values = pack_state(state, frozen)
    result = evaluate_residual(values, frozen, ledger=True)
    assert float(np.max(np.abs(result.gas_exchange_mol_s))) > 0.0  # 符号错误会被抓到
    rows = gas_residual_rows(frozen, unpack_state(values, frozen), result, result.ledger)
    for item in rows:
        cell, phase = int(item["cell"]), item["phase"]
        if item["zone"] == "bed":
            expected = result.gas_suspension_mol_s[cell] if phase == "suspension" else result.gas_bubble_mol_s[cell]
        else:
            expected = result.gas_free_mol_s[cell]
        assert np.allclose(item["residual"], expected, rtol=0.0, atol=1.0e-18), (item["zone"], cell, phase)
    assert any(float(np.max(np.abs(item["exchange"]))) > 0.0
               for item in rows if item["zone"] == "bed")
