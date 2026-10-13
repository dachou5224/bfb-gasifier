"""Q9：解包投影与下游隐式裁剪审计。

隔离的 ``strict_domain_no_projection`` 只去掉 unpack 的 max/fence。
速率函数内部把负浓度或 T<300 K 写成正值仍是投影，不能声称原式路径已干净。

Source: docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md Q9
"""

from __future__ import annotations

from typing import Literal, TypedDict

ClipKind = Literal[
    "state_projection",
    "maps_negative_concentration_to_zero",
    "maps_T_below_300",
    "overflow_guard",
    "geometric_floor",
]


class ClipRecord(TypedDict):
    path: str
    symbol: str
    kind: ClipKind
    on_strict_unpack_path: bool
    note: str


DOWNSTREAM_CLIPS: tuple[ClipRecord, ...] = (
    {
        "path": "src/solvers/nr_indexing.py",
        "symbol": "unpack_cell",
        "kind": "state_projection",
        "on_strict_unpack_path": False,
        "note": "生产 N=max(N,0) 与 300–2500 K fence；strict 入口不走此函数",
    },
    {
        "path": "src/kinetics/gas_reactions.py",
        "symbol": "rate_R5_suspension / rate_R5_bubble / rate_R6 / rate_R7 / rate_R8 / rate_R9 / rate_R12",
        "kind": "maps_negative_concentration_to_zero",
        "on_strict_unpack_path": True,
        "note": "合法域 C=0 边界可保留；负 C 仍被写成 0",
    },
    {
        "path": "src/kinetics/gas_reactions.py",
        "symbol": "rate_R5_bubble / rate_R6 / rate_R7 / rate_R8 / rate_R9 / rate_R12",
        "kind": "maps_T_below_300",
        "on_strict_unpack_path": True,
        "note": "max(T, 300) 不是 F 的定义域边界",
    },
    {
        "path": "src/kinetics/char_reactions.py",
        "symbol": "rate_R1 / rate_R2 / rate_R3 / extent_R4",
        "kind": "maps_negative_concentration_to_zero",
        "on_strict_unpack_path": True,
        "note": "负浓度/负分压裁零；extent_R4 还裁 m_char",
    },
    {
        "path": "src/kinetics/char_reactions.py",
        "symbol": "rate_R1",
        "kind": "maps_T_below_300",
        "on_strict_unpack_path": True,
        "note": "T = max(T, 300)",
    },
    {
        "path": "src/kinetics/arrhenius.py",
        "symbol": "k_jensen_h2",
        "kind": "maps_T_below_300",
        "on_strict_unpack_path": True,
        "note": "工厂函数内 T fence，R12 会走到",
    },
    {
        "path": "src/kinetics/arrhenius.py",
        "symbol": "k_hobbs / k_standard / k_jensen_r7 / _clamp_exp_arg",
        "kind": "overflow_guard",
        "on_strict_unpack_path": True,
        "note": "只裁指数幅值，不改写状态向量",
    },
    {
        "path": "src/kinetics/char_reactions.py",
        "symbol": "k_dg / rate_R1 denominators",
        "kind": "geometric_floor",
        "on_strict_unpack_path": True,
        "note": "粒径/扩散分母下限，不是物种投影",
    },
)


def unpack_projection_removed_from_strict_path() -> bool:
    """strict unpack 不再对状态做 max/fence。"""
    return all(
        (not rec["on_strict_unpack_path"])
        for rec in DOWNSTREAM_CLIPS
        if rec["symbol"] == "unpack_cell"
    )


def strict_path_still_projects_negative_state() -> bool:
    """下游速率仍会把负浓度变成 0；审计诚实，不是已完成声明。"""
    return any(
        rec["on_strict_unpack_path"] and rec["kind"] == "maps_negative_concentration_to_zero"
        for rec in DOWNSTREAM_CLIPS
    )
