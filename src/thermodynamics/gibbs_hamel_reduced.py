"""Hamel (1999) 附录 A1 风格的 lambda-only 降维 Gibbs 初始化求解器。"""

from __future__ import annotations

import itertools
from typing import Any, Dict, List, Protocol

import numpy as np
from scipy.optimize import least_squares

from src.core.constants import P0, Rg


class SpeciesDB(Protocol):
    """组分热力学数据库接口。"""

    def get_mu0(self, species: str, T: float) -> float:
        """标准化学势 μ°(T) [J/mol]。"""
        ...

    def get_atom_count(self, species: str, element: str) -> int:
        """组分中某元素的原子数。"""
        ...


class ReducedHamelGibbsSolver:
    """Hamel A1 降维 lambda-only Newton 求解器（major species 初始化用）。"""

    def __init__(self, species_db: SpeciesDB) -> None:
        self.species_db = species_db

    def solve(
        self,
        T: float,
        P: float,
        elements: Dict[str, float],
        candidates: List[str],
        tol: float = 1e-10,
        max_iter: int = 80,
        omega: float = 1.0,
        lambda0: np.ndarray | None = None,
        return_diag: bool = False,
    ) -> Dict[str, float] | tuple[Dict[str, float], Dict[str, Any]]:
        """求解入口（Hamel A1 Eq. A1-8/A1-27/A1-28 对应的降维系统）。"""
        orig_elem_names = list(elements.keys())
        orig_candidates = list(candidates)
        fixed_species: Dict[str, float] = {}
        elem_names = orig_elem_names
        candidates = orig_candidates

        n_elem = len(elem_names)
        n_sp = len(candidates)
        if n_elem < 2:
            raise ValueError("Reduced Hamel solver requires at least 2 tracked elements.")
        if n_sp <= 0:
            raise ValueError("Reduced Hamel solver requires non-empty candidate species.")

        a = self._build_atom_matrix(elem_names, candidates)
        assert a.shape == (n_elem, n_sp)

        b = np.array([float(elements[e]) for e in elem_names], dtype=np.float64)
        b = np.maximum(b, 0.0)
        mu0 = np.array([self.species_db.get_mu0(s, T) for s in candidates], dtype=np.float64)
        c = mu0 / (Rg * float(T)) + np.log(float(P) / P0)

        ref_idx = self._select_reference_element(a, b)
        active_idx = [i for i in range(n_elem) if i != ref_idx]
        active_arr = np.array(active_idx, dtype=np.int64)

        lam = np.zeros(n_elem, dtype=np.float64)
        if lambda0 is not None:
            lam0 = np.asarray(lambda0, dtype=np.float64).reshape(-1)
            if lam0.size == n_elem and np.all(np.isfinite(lam0)):
                lam = np.clip(lam0, -80.0, 80.0)
                lam[ref_idx] = 0.0
        lam = self._coarse_lambda_seed(
            lam=lam,
            a=a,
            c=c,
            b=b,
            ref_idx=ref_idx,
            active_idx=active_arr,
        )
        omega0 = float(np.clip(omega, 1e-6, 1.0))

        diag: Dict[str, Any] = {
            "n_iter": 0,
            "converged": False,
            "n_backtracks": 0,
            "final_residual": np.inf,
            "lambda": lam.copy(),
            "ln_N": float("nan"),
            "reference_element_idx": int(ref_idx),
            "reference_element_name": elem_names[ref_idx],
            "max_abs_element_closure": float("inf"),
        }

        for it in range(max_iter):
            diag["n_iter"] = int(it + 1)
            q = self._calc_q(lam, a, c)
            aq = self._calc_aq(a, q)
            f = self._reduced_residual(aq=aq, b=b, ref_idx=ref_idx, active_idx=active_arr)
            fnorm = self._scaled_residual_norm(f=f, aq=aq, b=b, ref_idx=ref_idx, active_idx=active_arr)
            diag["final_residual"] = fnorm
            if fnorm < tol:
                diag["converged"] = True
                break
            g = self._calc_g(a, q)
            j = self._reduced_jacobian(g=g, b=b, ref_idx=ref_idx, active_idx=active_arr)
            if not np.all(np.isfinite(j)):
                j = self._reduced_jacobian_fd(
                    lam=lam,
                    a=a,
                    c=c,
                    b=b,
                    ref_idx=ref_idx,
                    active_idx=active_arr,
                )
            row_scale = self._residual_row_scale(aq=aq, b=b, ref_idx=ref_idx, active_idx=active_arr)
            j_scaled = j / row_scale[:, None]
            f_scaled = f / row_scale
            try:
                du = np.linalg.solve(j_scaled, -f_scaled)
            except np.linalg.LinAlgError:
                du = np.linalg.lstsq(j_scaled, -f_scaled, rcond=None)[0]

            eta = omega0
            accepted = False
            for _ in range(20):
                lam_try = lam.copy()
                lam_try[active_arr] = lam_try[active_arr] + eta * du
                lam_try[ref_idx] = 0.0
                q_try = self._calc_q(lam_try, a, c)
                aq_try = self._calc_aq(a, q_try)
                f_try = self._reduced_residual(aq=aq_try, b=b, ref_idx=ref_idx, active_idx=active_arr)
                fn_try = self._scaled_residual_norm(
                    f=f_try,
                    aq=aq_try,
                    b=b,
                    ref_idx=ref_idx,
                    active_idx=active_arr,
                )
                if np.isfinite(fn_try) and fn_try <= fnorm * (1.0 - 1e-4 * eta):
                    lam = lam_try
                    accepted = True
                    break
                eta *= 0.5
                diag["n_backtracks"] = int(diag["n_backtracks"]) + 1
                if eta < 1e-8:
                    break
            if not accepted:
                omega0 *= 0.5
                if omega0 < 1e-4:
                    break

        if not bool(diag["converged"]):
            lam_ls, ls_norm = self._least_squares_refine(
                lam=lam,
                a=a,
                c=c,
                b=b,
                ref_idx=ref_idx,
                active_idx=active_arr,
            )
            if np.isfinite(ls_norm) and ls_norm < float(diag["final_residual"]):
                lam = lam_ls
                diag["final_residual"] = float(ls_norm)
            diag["used_least_squares"] = True
            if float(diag["final_residual"]) < tol:
                diag["converged"] = True

        q_final = self._calc_q(lam, a, c)
        aq_final = self._calc_aq(a, q_final)
        n_final, n_total = self._recover_moles_from_ratio(
            q=q_final,
            aq=aq_final,
            b=b,
            ref_idx=ref_idx,
        )
        out = {species: 0.0 for species in orig_candidates}
        for j, species in enumerate(candidates):
            out[species] = float(max(n_final[j], 0.0))
        for species, n_fix in fixed_species.items():
            out[species] = float(max(out.get(species, 0.0) + n_fix, 0.0))

        n_full = np.array([out[sp] for sp in orig_candidates], dtype=np.float64)
        a_full = self._build_atom_matrix(orig_elem_names, orig_candidates)
        b_full = np.array([float(elements[e]) for e in orig_elem_names], dtype=np.float64)
        closure = a_full @ n_full - b_full
        n_total_full = float(np.sum(n_full))

        diag["lambda"] = lam.copy()
        diag["ln_N"] = float(np.log(max(n_total_full, 1e-30)))
        diag["max_abs_element_closure"] = float(np.max(np.abs(closure), initial=0.0))
        diag["element_closure"] = {orig_elem_names[i]: float(closure[i]) for i in range(len(orig_elem_names))}
        diag["N_total"] = n_total_full
        total_elements = max(float(np.sum(np.maximum(b_full, 0.0))), 1e-12)
        closure_rel = float(diag["max_abs_element_closure"]) / total_elements
        diag["closure_rel"] = closure_rel
        if (
            not bool(diag.get("converged", False))
            and closure_rel <= 1e-3
            and float(diag.get("final_residual", np.inf)) <= 1e-2
        ):
            diag["converged"] = True
            diag["converged_by"] = "closure_relaxed"
        if fixed_species:
            diag["fixed_species"] = {k: float(v) for k, v in fixed_species.items()}
        if return_diag:
            return out, diag
        return out

    def _calc_q(self, lam: np.ndarray, a: np.ndarray, c: np.ndarray) -> np.ndarray:
        """按 Hamel Eq. A1-8 的指数项构造 q_j(lambda)。"""
        assert a.shape[0] == lam.size and a.shape[1] == c.size
        expo = np.clip(a.T @ lam - c, -80.0, 80.0)
        q = np.exp(expo)
        return np.maximum(q, 1e-300)

    def _calc_aq(self, a: np.ndarray, q: np.ndarray) -> np.ndarray:
        """元素加权和 A_i(lambda)=Σ_j a_ij q_j（Hamel A1 元素约束项）。"""
        assert a.shape[1] == q.size
        return a @ q

    def _calc_g(self, a: np.ndarray, q: np.ndarray) -> np.ndarray:
        """G_il=Σ_j a_ij a_lj q_j（用于 Hamel A1 降维 Jacobian 组装）。"""
        assert a.shape[1] == q.size
        return a @ (q[:, None] * a.T)

    def _reduced_residual(
        self,
        *,
        aq: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> np.ndarray:
        """降维残差：b_ref*A_i - b_i*A_ref = 0（Hamel A1 消元后方程）。"""
        assert aq.size == b.size
        a_ref = float(aq[ref_idx])
        b_ref = float(b[ref_idx])
        return b_ref * aq[active_idx] - b[active_idx] * a_ref

    def _reduced_jacobian(
        self,
        *,
        g: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> np.ndarray:
        """降维 Jacobian（对应 Hamel A1 Eq. A1-27/A1-28 的数值组装）。"""
        n_act = active_idx.size
        j = np.zeros((n_act, n_act), dtype=np.float64)
        for ri, i in enumerate(active_idx):
            for cj, l in enumerate(active_idx):
                j[ri, cj] = float(b[ref_idx]) * g[i, l] - float(b[i]) * g[ref_idx, l]
        return j

    def _reduced_jacobian_fd(
        self,
        *,
        lam: np.ndarray,
        a: np.ndarray,
        c: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> np.ndarray:
        """Finite-difference Jacobian for reduced equations (robust fallback)."""
        n = active_idx.size
        j = np.zeros((n, n), dtype=np.float64)
        h = 1e-6
        for col, l in enumerate(active_idx):
            lam_p = lam.copy()
            lam_m = lam.copy()
            lam_p[l] += h
            lam_m[l] -= h
            lam_p[ref_idx] = 0.0
            lam_m[ref_idx] = 0.0
            fp = self._reduced_residual(
                aq=self._calc_aq(a, self._calc_q(lam_p, a, c)),
                b=b,
                ref_idx=ref_idx,
                active_idx=active_idx,
            )
            fm = self._reduced_residual(
                aq=self._calc_aq(a, self._calc_q(lam_m, a, c)),
                b=b,
                ref_idx=ref_idx,
                active_idx=active_idx,
            )
            j[:, col] = (fp - fm) / (2.0 * h)
        return j

    def _recover_moles_from_ratio(
        self,
        *,
        q: np.ndarray,
        aq: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
    ) -> tuple[np.ndarray, float]:
        """由 z=A_ref/b_ref 与 n_j=q_j/z 恢复摩尔流（Hamel A1 消元回代）。"""
        a_ref = max(float(aq[ref_idx]), 1e-300)
        b_ref = max(float(b[ref_idx]), 1e-300)
        z = a_ref / b_ref
        n = q / max(z, 1e-300)
        n = np.maximum(n, 0.0)
        n_total = float(np.sum(n))
        return n, n_total

    def _scaled_residual_norm(
        self,
        *,
        f: np.ndarray,
        aq: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> float:
        """Residual scaled norm for robust stopping in large-magnitude atom pools."""
        scale = self._residual_row_scale(aq=aq, b=b, ref_idx=ref_idx, active_idx=active_idx)
        return float(np.linalg.norm(f / scale))

    def _residual_row_scale(
        self,
        *,
        aq: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> np.ndarray:
        """Row scaling for reduced residual/Jacobian equations."""
        a_ref = float(aq[ref_idx])
        b_ref = float(b[ref_idx])
        return np.abs(b_ref * aq[active_idx]) + np.abs(b[active_idx] * a_ref) + 1.0

    def _extract_fixed_singleton_species(
        self,
        *,
        elements: Dict[str, float],
        candidates: List[str],
    ) -> tuple[Dict[str, float], List[str], Dict[str, float]]:
        """Extract pure-singleton species fixed directly by one elemental pool."""
        rem_elements = {k: float(max(v, 0.0)) for k, v in elements.items()}
        rem_candidates = list(candidates)
        fixed: Dict[str, float] = {}

        changed = True
        while changed and rem_candidates and rem_elements:
            changed = False
            elem_names = list(rem_elements.keys())
            for elem in elem_names:
                if rem_elements[elem] <= 1e-18:
                    continue
                participating = [
                    sp for sp in rem_candidates if self.species_db.get_atom_count(sp, elem) > 0
                ]
                if len(participating) != 1:
                    continue
                sp = participating[0]
                if any(
                    self.species_db.get_atom_count(sp, e2) > 0 for e2 in elem_names if e2 != elem
                ):
                    continue
                atom_n = float(self.species_db.get_atom_count(sp, elem))
                if atom_n <= 0.0:
                    continue
                n_fix = rem_elements[elem] / atom_n
                fixed[sp] = fixed.get(sp, 0.0) + n_fix
                rem_candidates = [name for name in rem_candidates if name != sp]
                for e2 in elem_names:
                    rem_elements[e2] = max(
                        rem_elements[e2] - n_fix * float(self.species_db.get_atom_count(sp, e2)),
                        0.0,
                    )
                changed = True
                break
            if changed:
                rem_elements = {k: v for k, v in rem_elements.items() if v > 1e-18}

        return rem_elements, rem_candidates, fixed

    def _coarse_lambda_seed(
        self,
        *,
        lam: np.ndarray,
        a: np.ndarray,
        c: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> np.ndarray:
        """Grid-search a robust initial lambda in highly non-convex cases."""
        if active_idx.size == 0:
            return lam
        if active_idx.size <= 2:
            grid = (
                -60.0,
                -50.0,
                -45.0,
                -40.0,
                -39.0,
                -35.0,
                -30.0,
                -25.0,
                -20.0,
                -15.0,
                -10.0,
                -5.0,
                0.0,
                5.0,
                10.0,
                20.0,
            )
        else:
            grid = (-60.0, -45.0, -35.0, -25.0, -15.0, -5.0, 5.0, 15.0)
        best = lam.copy()
        q0 = self._calc_q(best, a, c)
        aq0 = self._calc_aq(a, q0)
        f0 = self._reduced_residual(aq=aq0, b=b, ref_idx=ref_idx, active_idx=active_idx)
        best_norm = self._scaled_residual_norm(f=f0, aq=aq0, b=b, ref_idx=ref_idx, active_idx=active_idx)
        for vals in itertools.product(grid, repeat=active_idx.size):
            trial = lam.copy()
            trial[active_idx] = np.array(vals, dtype=np.float64)
            trial[ref_idx] = 0.0
            q = self._calc_q(trial, a, c)
            aq = self._calc_aq(a, q)
            f = self._reduced_residual(aq=aq, b=b, ref_idx=ref_idx, active_idx=active_idx)
            fn = self._scaled_residual_norm(f=f, aq=aq, b=b, ref_idx=ref_idx, active_idx=active_idx)
            if np.isfinite(fn) and fn < best_norm:
                best = trial
                best_norm = fn
        return best

    def _least_squares_refine(
        self,
        *,
        lam: np.ndarray,
        a: np.ndarray,
        c: np.ndarray,
        b: np.ndarray,
        ref_idx: int,
        active_idx: np.ndarray,
    ) -> tuple[np.ndarray, float]:
        """Trust-region least-squares refinement for hard reduced-system cells."""
        def _fun_full(x: np.ndarray) -> np.ndarray:
            trial = np.asarray(x, dtype=np.float64)
            q = self._calc_q(trial, a, c)
            aq = self._calc_aq(a, q)
            f_ratio = self._reduced_residual(aq=aq, b=b, ref_idx=ref_idx, active_idx=active_idx)
            s_ratio = self._residual_row_scale(aq=aq, b=b, ref_idx=ref_idx, active_idx=active_idx)
            f_norm = np.array([float(np.sum(q) - 1.0)], dtype=np.float64)
            s_norm = np.array([float(np.sum(q) + 1.0)], dtype=np.float64)
            return np.concatenate([f_ratio / s_ratio, f_norm / s_norm])

        starts: list[np.ndarray] = [np.array(lam, dtype=np.float64), np.zeros(lam.size, dtype=np.float64)]
        if lam.size <= 4:
            for vals in itertools.product((0.0, -20.0), repeat=lam.size):
                starts.append(np.array(vals, dtype=np.float64))
        else:
            starts.append(np.full(lam.size, -20.0, dtype=np.float64))

        # de-duplicate starts while preserving order
        dedup: list[np.ndarray] = []
        seen: set[tuple[float, ...]] = set()
        for s in starts:
            key = tuple(np.round(s.astype(float), 12))
            if key in seen:
                continue
            seen.add(key)
            dedup.append(s)

        best_lam = lam.copy()
        best_norm = float(np.linalg.norm(_fun_full(best_lam)))
        for x0 in dedup:
            try:
                sol = least_squares(
                    _fun_full,
                    x0=x0,
                    method="trf",
                    max_nfev=2000,
                    ftol=1e-12,
                    xtol=1e-12,
                    gtol=1e-12,
                )
            except Exception:
                continue
            trial_lam = np.asarray(sol.x, dtype=np.float64)
            fn = float(np.linalg.norm(_fun_full(trial_lam)))
            if np.isfinite(fn) and fn < best_norm:
                best_norm = fn
                best_lam = trial_lam
        return best_lam, best_norm

    def _select_reference_element(self, a: np.ndarray, b: np.ndarray) -> int:
        """选择消元参考元素（优先选择参与多组分的元素，避免退化参考行）。"""
        positive = np.where(b > 0.0)[0]
        if positive.size == 0:
            raise ValueError("All element pools are zero; cannot build reduced Hamel system.")
        participation = np.sum(a[positive, :] > 1e-18, axis=1)
        mixed = positive[participation >= 2]
        if mixed.size > 0:
            return int(mixed[np.argmax(b[mixed])])
        return int(positive[np.argmax(b[positive])])

    def _build_atom_matrix(self, elem_names: List[str], candidates: List[str]) -> np.ndarray:
        """原子矩阵 a[i,j]（Hamel A1 基础量）。"""
        a = np.zeros((len(elem_names), len(candidates)), dtype=np.float64)
        for i, elem in enumerate(elem_names):
            for j, species in enumerate(candidates):
                a[i, j] = float(self.species_db.get_atom_count(species, elem))
        return a
