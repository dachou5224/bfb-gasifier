"""具名 F 登记：一次取消一项省略。不切生产默认。

Source: docs/audit_reports/wirsum_af_open_questions_answers_2026-09-14.md
"""

from __future__ import annotations

OMITTED_BASE: tuple[str, ...] = (
    "energy_wall_heat",
    "cyclone_K_update_P",
    "K_eject_in_this_P",
    "DAEM_X_from_T",
    "R2_R12",
    "tar_transport",
    "ash_inventory",
    "freeboard_reaction",
    "return_leg_holdup",
    "production_hamel_outer_refresh",
)

# 阶段 k 取消的省略名；旋风 K 与生产 refresh 永不取消。
_CANCEL: dict[int, str] = {
    1: "R2_R4",
    2: "R5_R8_R12",
    3: "ash_inventory",
    4: "DAEM_X_from_T",
    5: "energy_wall_heat",
    6: "tar_transport",
    7: "freeboard_reaction",
    8: "K_eject_in_this_P",
}

_STAGE_NAMES: tuple[str, ...] = (
    "F_E0",
    "F_E1",
    "F_E2",
    "F_E3",
    "F_E4",
    "F_E5",
    "F_E6",
    "F_E7",
    "F_E8",
    "F_E9",
    "F_E10",
    "F_E11",
)


def f_name(stage: int) -> str:
    if not 0 <= int(stage) <= 11:
        raise ValueError(f"named-F stage {stage} not in 0..8")
    return _STAGE_NAMES[int(stage)]


def n_unknowns(stage: int) -> int:
    """40 + 灰库存 + 焦油。温度始终在最后。

    E6 起两个焦油代表物各一个出口。E8 再给每个代表物一个气泡相未知量，
    悬浮相仍占 E6 的那两个槽。Hamel §2.2：气体物种两相分别衡算。
    """
    n = 40
    if int(stage) >= 3:
        n += 1
    if int(stage) >= 6:
        n += 2
    if int(stage) >= 8:
        n += 2
    if int(stage) >= 9:
        n += 1
    if int(stage) >= 10:
        n += 2
    if int(stage) >= 11:
        n += 2
    return n


def omitted_for_stage(stage: int) -> tuple[str, ...]:
    """本阶段仍省略的项。R2_R12 在 E1 拆成剩余均相。"""
    s = int(stage)
    out: list[str] = []
    for name in OMITTED_BASE:
        if name == "R2_R12":
            if s < 1:
                out.append("R2_R12")
            elif s < 2:
                out.append("R5_R8_R12")
            continue
        if name == "ash_inventory" and s >= 3:
            continue
        if name == "DAEM_X_from_T" and s >= 4:
            continue
        if name == "energy_wall_heat" and s >= 5:
            continue
        if name == "tar_transport" and s >= 6:
            continue
        if name == "freeboard_reaction" and s >= 7:
            continue
        if name == "K_eject_in_this_P" and s >= 8:
            continue
        out.append(name)
    return tuple(out)


def cancelled_this_stage(stage: int) -> str | None:
    return _CANCEL.get(int(stage))


def omitted_named(stage: int) -> dict[str, str]:
    return {name: "omitted, named" for name in omitted_for_stage(stage)}
