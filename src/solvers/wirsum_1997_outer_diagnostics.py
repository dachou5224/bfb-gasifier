"""Wirsum (1997) §2.2 Vorabrechnung outer map 诊断。

记录

    h = (ε_b, K_bd, K_s, u_b, T)

在相邻 outer 步的欧氏增量，并估计

    ρ ≈ ||h_{k+1} − h_k|| / ||h_k − h_{k−1}||

不改水力公式，不加 under-relaxation / Anderson / inner clip。
生产 ``outer_loop.py`` 仍走工程 Abgleich；本模块是 reduced / 对照核。

Source: Wirsum (1997) §2.2 Vorabrechnung / Zellenmodell 交替
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import numpy.typing as npt

FIELD_NAMES: tuple[str, ...] = ("eps_b", "K_bd", "K_s", "u_b", "T")
HYDRO_FIELDS: tuple[str, ...] = ("eps_b", "K_bd", "K_s", "u_b")

# 仅用于拼向量后的单一 ρ；逐场 ||Δ|| 仍是主输出（量纲不同，不混用）
_FIELD_SCALE: dict[str, float] = {
    "eps_b": 1.0,  # [-]
    "K_bd": 1.0,  # [1/s]
    "K_s": 1.0,  # [1/s]
    "u_b": 1.0,  # [m/s]
    "T": 100.0,  # [K]  拼向量时除以 100 K，避免温度主导
}


@dataclass(frozen=True)
class VorabFields:
    """单次 outer 的 Vorabrechnung 场。长度均为 n_cells。"""

    eps_b: npt.NDArray[np.float64]  # [-]
    K_bd: npt.NDArray[np.float64]  # [1/s]
    K_s: npt.NDArray[np.float64]  # [1/s]  per-cell 固相输运系数
    u_b: npt.NDArray[np.float64]  # [m/s]
    T: npt.NDArray[np.float64]  # [K]

    def __post_init__(self) -> None:
        n = int(np.asarray(self.eps_b).size)
        for name in FIELD_NAMES:
            arr = np.asarray(getattr(self, name), dtype=np.float64)
            object.__setattr__(self, name, arr.reshape(-1).astype(np.float64, copy=False))
            if int(arr.size) != n:
                raise ValueError(f"{name} size {arr.size} != eps_b size {n}")

    @property
    def n_cells(self) -> int:
        return int(self.eps_b.size)

    def as_lists(self) -> dict[str, list[float]]:
        return {name: np.asarray(getattr(self, name), dtype=np.float64).tolist() for name in FIELD_NAMES}


def _as_1d(values: Any) -> npt.NDArray[np.float64]:
    return np.asarray(values, dtype=np.float64).reshape(-1)


def make_vorab_fields(
    *,
    eps_b: npt.ArrayLike,
    K_bd: npt.ArrayLike,
    K_s: npt.ArrayLike,
    u_b: npt.ArrayLike,
    T: npt.ArrayLike,
) -> VorabFields:
    """从同长度数组构造场快照。"""
    return VorabFields(
        eps_b=_as_1d(eps_b),
        K_bd=_as_1d(K_bd),
        K_s=_as_1d(K_s),
        u_b=_as_1d(u_b),
        T=_as_1d(T),
    )


def cell_solid_transport_K_s(cell: Any) -> float:
    """单格 ``K_s`` [1/s]。

    优先 ``mean(K_solid_auf) + mean(K_solid_ab)``；否则读标量 ``K_s``。
    这是诊断用的固相输运强度，不是新的本构。
    """
    auf = getattr(cell, "K_solid_auf", None)
    ab = getattr(cell, "K_solid_ab", None)
    parts: list[float] = []
    if auf is not None:
        parts.append(float(np.mean(np.asarray(auf, dtype=np.float64))))
    if ab is not None:
        parts.append(float(np.mean(np.asarray(ab, dtype=np.float64))))
    if parts:
        return float(np.sum(parts))
    return float(getattr(cell, "K_s", 0.0))


def snapshot_vorab_fields(cells: Sequence[Any]) -> VorabFields:
    """从 cell 列表读场。不刷新水力，不写回 cell。"""
    if len(cells) == 0:
        raise ValueError("snapshot_vorab_fields requires at least one cell")
    return make_vorab_fields(
        eps_b=[float(getattr(c, "eps_b")) for c in cells],
        K_bd=[float(getattr(c, "K_bd")) for c in cells],
        K_s=[cell_solid_transport_K_s(c) for c in cells],
        u_b=[float(getattr(c, "u_b")) for c in cells],
        T=[float(getattr(c, "T")) for c in cells],
    )


def field_l2_deltas(prev: VorabFields, curr: VorabFields) -> dict[str, float]:
    """逐场 ``||f^{k+1} − f^k||_2``。"""
    if prev.n_cells != curr.n_cells:
        raise ValueError("VorabFields cell counts differ")
    out: dict[str, float] = {}
    for name in FIELD_NAMES:
        delta = np.asarray(getattr(curr, name)) - np.asarray(getattr(prev, name))
        out[name] = float(np.linalg.norm(delta))
    return out


def scaled_concat_delta(prev: VorabFields, curr: VorabFields) -> float:
    """量纲归一后的拼接欧氏增量。只用于单一 ρ，不替代逐场范数。"""
    chunks: list[npt.NDArray[np.float64]] = []
    for name in FIELD_NAMES:
        scale = _FIELD_SCALE[name]
        delta = (np.asarray(getattr(curr, name)) - np.asarray(getattr(prev, name))) / scale
        chunks.append(delta)
    return float(np.linalg.norm(np.concatenate(chunks)))


def contractivity_ratio(delta_next: float, delta_prev: float) -> float | None:
    """``ρ = ||Δh_{k+1}|| / ||Δh_k||``。上一步增量为 0 时无法定义。"""
    prev = float(delta_prev)
    nxt = float(delta_next)
    if not np.isfinite(prev) or not np.isfinite(nxt) or prev <= 0.0:
        return None
    return nxt / prev


def apply_linear_outer_map(h: VorabFields, h_star: VorabFields, rho: float) -> VorabFields:
    """仿射 outer map ``h ↦ h* + ρ (h − h*)``。ρ<1 收缩，ρ>1 扩张。"""
    kwargs = {
        name: np.asarray(getattr(h_star, name))
        + float(rho) * (np.asarray(getattr(h, name)) - np.asarray(getattr(h_star, name)))
        for name in FIELD_NAMES
    }
    return make_vorab_fields(**kwargs)


@dataclass
class OuterIterationRecorder:
    """累积 outer 场快照与 inner residual，计算增量 / 收缩比。"""

    inner_residual_kind: str = "raw_l2"
    snapshots: list[VorabFields] = field(default_factory=list)
    inner_residuals: list[float] = field(default_factory=list)

    def record(self, fields: VorabFields, inner_residual: float) -> dict[str, Any]:
        self.snapshots.append(fields)
        self.inner_residuals.append(float(inner_residual))
        return self._step_payload(len(self.snapshots) - 1)

    def _step_payload(self, index: int) -> dict[str, Any]:
        curr = self.snapshots[index]
        payload: dict[str, Any] = {
            "outer_iter": index,
            "n_cells": curr.n_cells,
            "inner_residual": float(self.inner_residuals[index]),
            "inner_residual_kind": self.inner_residual_kind,
            "snapshot": curr.as_lists(),
            "delta": None,
            "delta_concat_scaled": None,
            "contractivity": None,
            "contractivity_concat": None,
        }
        if index == 0:
            return payload
        prev = self.snapshots[index - 1]
        delta = field_l2_deltas(prev, curr)
        payload["delta"] = delta
        payload["delta_concat_scaled"] = scaled_concat_delta(prev, curr)
        if index >= 2:
            prev2 = self.snapshots[index - 2]
            delta_prev = field_l2_deltas(prev2, prev)
            payload["contractivity"] = {
                name: contractivity_ratio(delta[name], delta_prev[name]) for name in FIELD_NAMES
            }
            payload["contractivity_concat"] = contractivity_ratio(
                float(payload["delta_concat_scaled"]),
                scaled_concat_delta(prev2, prev),
            )
        return payload

    def steps(self) -> list[dict[str, Any]]:
        return [self._step_payload(i) for i in range(len(self.snapshots))]

    def verdict(self) -> dict[str, Any]:
        """用后半段水力场 ρ 判断 map 是否扩张。不据此改求解器。"""
        ratios: list[float] = []
        hydro_ratios: list[float] = []
        for step in self.steps():
            contr = step.get("contractivity") or {}
            rho_c = step.get("contractivity_concat")
            if rho_c is not None and np.isfinite(float(rho_c)):
                ratios.append(float(rho_c))
            hydro_vals = [
                float(contr[name])
                for name in HYDRO_FIELDS
                if contr.get(name) is not None and np.isfinite(float(contr[name]))
            ]
            if hydro_vals:
                hydro_ratios.append(float(max(hydro_vals)))
        sample = hydro_ratios[-3:] if len(hydro_ratios) >= 3 else hydro_ratios
        mean_rho = float(np.mean(sample)) if sample else None
        max_rho = float(np.max(sample)) if sample else None
        return {
            "n_contractivity_samples": len(hydro_ratios),
            "mean_hydro_rho_tail": mean_rho,
            "max_hydro_rho_tail": max_rho,
            "mean_concat_rho": float(np.mean(ratios)) if ratios else None,
            "non_contractive": bool(mean_rho is not None and mean_rho > 1.0),
        }

    def payload(self, *, label: str = "") -> dict[str, Any]:
        return {
            "source": "Wirsum (1997) §2.2 Vorabrechnung outer map",
            "label": label,
            "fields": list(FIELD_NAMES),
            "hydro_fields": list(HYDRO_FIELDS),
            "norm": "euclidean_l2",
            "inner_residual_kind": self.inner_residual_kind,
            "steps": self.steps(),
            "verdict": self.verdict(),
        }

    def write_json(self, path: str | Path, *, label: str = "") -> Path:
        out = Path(path)
        out.write_text(json.dumps(self.payload(label=label), indent=2) + "\n", encoding="utf-8")
        return out


def run_linear_outer_history(
    *,
    h0: VorabFields,
    h_star: VorabFields,
    rho: float,
    n_outer: int,
    inner_residual0: float,
    inner_residual_decay: float,
    inner_residual_kind: str = "raw_l2",
) -> OuterIterationRecorder:
    """对仿射 map 迭代 ``n_outer`` 步，inner residual 按几何衰减（可与 ρ 脱钩）。"""
    rec = OuterIterationRecorder(inner_residual_kind=inner_residual_kind)
    h = h0
    residual = float(inner_residual0)
    rec.record(h, residual)
    for _ in range(int(n_outer)):
        h = apply_linear_outer_map(h, h_star, rho)
        residual = residual * float(inner_residual_decay)
        rec.record(h, residual)
    return rec


def build_outer_map_ab_payload() -> dict[str, Any]:
    """三格仿射 map：收缩 ρ=0.4 vs 扩张 ρ=1.3；两边 inner residual 都减半。"""
    h_star = make_vorab_fields(
        eps_b=[0.20, 0.18, 0.15],
        K_bd=[4.0, 3.5, 3.0],
        K_s=[0.40, 0.30, 0.20],
        u_b=[0.80, 0.85, 0.90],
        T=[1100.0, 1080.0, 1050.0],
    )
    h0 = make_vorab_fields(
        eps_b=[0.35, 0.30, 0.28],
        K_bd=[1.0, 1.2, 1.5],
        K_s=[0.10, 0.10, 0.10],
        u_b=[0.40, 0.45, 0.50],
        T=[900.0, 920.0, 940.0],
    )
    contractive = run_linear_outer_history(
        h0=h0,
        h_star=h_star,
        rho=0.4,
        n_outer=8,
        inner_residual0=1.0,
        inner_residual_decay=0.5,
    )
    expanding = run_linear_outer_history(
        h0=h0,
        h_star=h_star,
        rho=1.3,
        n_outer=8,
        inner_residual0=1.0,
        inner_residual_decay=0.5,
    )
    return {
        "source": "Wirsum (1997) §2.2 reduced outer map A/B",
        "n_cells": 3,
        "inner_residual_decay": 0.5,
        "note": (
            "inner residual 两边都几何衰减；扩张 map 的 non_contractive=true "
            "说明 inner 收敛不能代替 outer fixed-point 收缩。"
        ),
        "contractive_rho_0p4": contractive.payload(label="rho=0.4"),
        "expanding_rho_1p3": expanding.payload(label="rho=1.3"),
    }
