"""Gibbs 自由焓最小化求解器。

Hamel 附录 A1（pp. 163–168）/ TechSpec §5.6。
拉格朗日乘子法 + 带阻尼 Newton-Raphson 迭代。
"""

from __future__ import annotations

from typing import Any, Dict, List, Protocol

import warnings

import numpy as np

from src.core.constants import Rg, P0
from src.core.hamel_parity_defaults import (
    HAMEL_GIBBS_INITIAL_GUESS_MODE,
    HAMEL_GIBBS_TOL_MODE,
    LEGACY_GIBBS_INITIAL_GUESS_MODE,
    LEGACY_GIBBS_TOL_MODE,
)

GIBBS_TOL_MODES = frozenset({"l2_residual", "holub_vonka_a34"})
GIBBS_INITIAL_GUESS_MODES = frozenset({"taylor_linearized", "hamel_a40"})


def normalize_gibbs_initial_guess_mode(mode: str | None) -> str:
    """Normalize Gibbs initial-guess construction mode."""
    guess_mode = str(mode or HAMEL_GIBBS_INITIAL_GUESS_MODE).strip().lower()
    if guess_mode not in GIBBS_INITIAL_GUESS_MODES:
        raise ValueError(
            "Unsupported initial_guess_mode="
            f"{mode!r}; expected one of {sorted(GIBBS_INITIAL_GUESS_MODES)}"
        )
    return guess_mode


def normalize_gibbs_tol_mode(tol_mode: str | None) -> str:
    """Normalize Gibbs solver convergence mode."""
    mode = str(tol_mode or HAMEL_GIBBS_TOL_MODE).strip().lower()
    if mode not in GIBBS_TOL_MODES:
        raise ValueError(
            f"Unsupported tol_mode={tol_mode!r}; expected one of {sorted(GIBBS_TOL_MODES)}"
        )
    return mode


def holub_vonka_a34_relative_change(
    new: np.ndarray,
    old: np.ndarray,
    *,
    floor: float = 1e-300,
) -> float:
    """Holub & Vonka (1976) / Hamel (1999) (A.34) 相对变化判据。

    η = Σ_j |x_j^{(k+1)} − x_j^{(k)}| / max(Σ_j x_j^{(k)}, floor)

    增广路径对摩尔量 n_j 使用；降阶路径对 (A.23) 的 q_j 使用。
    """
    old_arr = np.maximum(np.asarray(old, dtype=np.float64), 0.0)
    new_arr = np.maximum(np.asarray(new, dtype=np.float64), 0.0)
    if old_arr.shape != new_arr.shape:
        raise ValueError("holub_vonka_a34_relative_change: shape mismatch")
    denom = max(float(np.sum(old_arr)), float(floor))
    return float(np.sum(np.abs(new_arr - old_arr)) / denom)


def hamel_a40_initial_guess(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    *,
    max_iter: int = 25,
    tol: float = 1e-10,
) -> tuple[np.ndarray, float]:
    """Hamel (1999) (A.35)–(A.40) 分步初值：y_j^0 → 线性化 (A.40) → λ / ln N。

    (A.35) 由 (A.15) 在 λ=0 时得 y_j^0 ∝ exp(−c_j)。
    (A.38) 摩尔份额守恒 Σ_j a_ij y_j = b_i/N^0，N^0 由 y^0 估计。
    (A.36)–(A.40) 在 y^0 处做一次线性化 Newton（等价于对 n_j=N^0 y_j 求 λ），
    避免对 y_j→0 组分做完整非线性 (A.40) 迭代时的 ln(y) 奇点。

    ``max_iter`` / ``tol`` 保留接口；当前实现为单步线性 (A.40)。
    """
    _ = max_iter, tol
    n_elem, n_sp = a.shape
    if n_elem <= 0 or n_sp <= 0:
        raise ValueError("hamel_a40_initial_guess requires non-empty atom matrix")

    # (A.35)
    weights = np.exp(np.clip(-c, -50.0, 50.0))
    y0 = weights / (float(np.sum(weights)) + 1e-30)

    # (A.38): N^0 与 b_i^* = b_i/N^0
    ay0 = a @ y0
    positive = np.where((b > 0.0) & (ay0 > 1e-300))[0]
    if positive.size == 0:
        n_scale = max(float(np.sum(np.abs(b))), 1e-20)
    else:
        n_scale = float(np.max(b[positive] / ay0[positive]))
    n0 = n_scale * y0

    # (A.40) 线性化一步 → λ（与增广链 (A.18) 在 n0 处 Taylor 一致）
    lhs = a @ (n0[:, None] * a.T) + 1e-10 * np.eye(n_elem)
    rhs = b - a @ n0
    try:
        lam, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)
    except np.linalg.LinAlgError:
        lam = np.zeros(n_elem, dtype=np.float64)
    lam = np.clip(np.asarray(lam, dtype=np.float64).reshape(-1), -80.0, 80.0)
    if lam.size != n_elem:
        lam = np.zeros(n_elem, dtype=np.float64)

    ln_n = float(np.log(max(n_scale, 1e-300)))
    return lam, ln_n


class SpeciesDB(Protocol):
    """组分热力学数据库接口。"""

    def get_mu0(self, species: str, T: float) -> float:
        """标准化学势 μ°(T) [J/mol]。"""
        ...

    def get_atom_count(self, species: str, element: str) -> int:
        """组分中某元素的原子数。"""
        ...


class GibbsMinimizer:
    """Hamel 附录 A1：Gibbs 自由焓最小化求解器。"""

    def __init__(self, species_db: SpeciesDB) -> None:
        self.species_db = species_db

    def solve(
        self,
        T: float,
        P: float,
        elements: Dict[str, float],
        candidates: List[str],
        tol: float = 1e-10,
        max_iter: int = 100,
        omega: float = 1.0,
        lambda0: np.ndarray | None = None,
        ln_N0: float | None = None,
        eta_strategy: str = "feasible_backtracking",
        tol_mode: str = HAMEL_GIBBS_TOL_MODE,
        initial_guess_mode: str = HAMEL_GIBBS_INITIAL_GUESS_MODE,
        return_diag: bool = False,
    ) -> Dict[str, float] | tuple[Dict[str, float], Dict[str, Any]]:
        """主求解入口。

        Parameters
        ----------
        T : float
            温度 [K]
        P : float
            压力 [Pa]
        elements : dict
            元素总量 {元素名: mol/s}
        candidates : list
            候选气相组分名
        tol, max_iter, omega
            收敛容差、最大迭代次数、初始阻尼因子
        lambda0, ln_N0
            可选 warm-start 初值（拉格朗日乘子与 ln(N)）
        eta_strategy
            阻尼策略：`feasible_backtracking`（默认）或 `legacy_half`
        tol_mode
            收敛判据：默认 ``holub_vonka_a34``（(A.34)）；审计对照 ``l2_residual``
        initial_guess_mode
            初值构造：默认 ``hamel_a40``（(A.35)–(A.40)）；审计对照 ``taylor_linearized``
        return_diag
            True 时返回 (result, diagnostics)

        Returns
        -------
        dict
            {组分名: 摩尔流率 [mol/s]}
        """
        elem_names = list(elements.keys())
        N_E = len(elem_names)
        N_s = len(candidates)

        a = self._build_atom_matrix(elem_names, candidates)
        b = np.array([elements[e] for e in elem_names])

        mu0 = np.array([self.species_db.get_mu0(s, T) for s in candidates])
        c = mu0 / (Rg * T) + np.log(P / P0)

        guess_mode_norm = normalize_gibbs_initial_guess_mode(initial_guess_mode)
        if lambda0 is None:
            if guess_mode_norm == "hamel_a40":
                lam, ln_n_guess = hamel_a40_initial_guess(a, b, c)
                lam = np.clip(lam, -80.0, 80.0)
            else:
                lam = self._initial_guess_taylor_linearized(a, b, c, N_E)
                ln_n_guess = None
        else:
            lam = self._initial_guess_taylor_linearized(a, b, c, N_E)
            ln_n_guess = None
        if lambda0 is not None:
            lam0 = np.asarray(lambda0, dtype=np.float64).reshape(-1)
            if lam0.size == N_E and np.all(np.isfinite(lam0)):
                lam = lam0.copy()
        b_sum = max(float(np.sum(np.abs(b))), 1e-20)
        if ln_N0 is not None:
            ln_N = float(np.clip(ln_N0, -50.0, 50.0))
        elif ln_n_guess is not None:
            ln_N = float(np.clip(ln_n_guess, -50.0, 50.0))
        else:
            ln_N = np.log(b_sum)

        tol_mode_norm = normalize_gibbs_tol_mode(tol_mode)
        diag: Dict[str, Any] = {
            "n_iter": 0,
            "converged": False,
            "n_backtracks": 0,
            "eta_strategy": eta_strategy,
            "tol_mode": tol_mode_norm,
            "initial_guess_mode": guess_mode_norm,
            "final_residual": np.inf,
            "holub_vonka_a34": np.inf,
            "lambda": None,
            "ln_N": None,
        }
        omega0 = float(np.clip(omega, 1e-6, 1.0))

        for it in range(max_iter):
            diag["n_iter"] = int(it + 1)
            n, N = self._calc_moles(lam, ln_N, c, a, N_s)
            if not np.all(np.isfinite(n)) or N <= 0 or not np.isfinite(N):
                omega0 *= 0.5
                if omega0 < 1e-3:
                    break
                continue
            f = self._residuals(n, N, a, b, N_E)
            fnorm = float(np.linalg.norm(f))
            diag["final_residual"] = fnorm
            if tol_mode_norm == "l2_residual" and fnorm < tol:
                diag["converged"] = True
                diag["converged_by"] = "l2_residual"
                break
            J = self._jacobian(n, N, a, b, N_E)
            if not np.all(np.isfinite(J)):
                omega0 *= 0.5
                continue
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    dx = np.linalg.solve(J, -f)
            except np.linalg.LinAlgError:
                omega0 *= 0.5
                continue

            if eta_strategy == "legacy_half":
                lam_new = lam + omega0 * dx[:N_E]
                ln_N_new = float(np.clip(ln_N + omega0 * dx[N_E], -50.0, 50.0))
                if tol_mode_norm == "holub_vonka_a34":
                    n_new, _ = self._calc_moles(lam_new, ln_N_new, c, a, N_s)
                    eta_a34 = holub_vonka_a34_relative_change(n_new, n)
                    diag["holub_vonka_a34"] = eta_a34
                    if eta_a34 < tol:
                        diag["converged"] = True
                        diag["converged_by"] = "holub_vonka_a34"
                        lam, ln_N = lam_new, ln_N_new
                        break
                lam, ln_N = lam_new, ln_N_new
                continue

            # Holub/Vonka 风格：在可行域内回溯 eta，确保残差下降。
            eta = omega0
            accepted = False
            for _ in range(12):
                lam_try = lam + eta * dx[:N_E]
                ln_N_try = float(np.clip(ln_N + eta * dx[N_E], -50.0, 50.0))
                n_try, N_try = self._calc_moles(lam_try, ln_N_try, c, a, N_s)
                if np.all(np.isfinite(n_try)) and np.all(n_try > 0.0) and np.isfinite(N_try) and N_try > 0.0:
                    f_try = self._residuals(n_try, N_try, a, b, N_E)
                    fn_try = float(np.linalg.norm(f_try))
                    if np.isfinite(fn_try) and fn_try <= fnorm * (1.0 - 1e-4 * eta):
                        lam, ln_N = lam_try, ln_N_try
                        accepted = True
                        if tol_mode_norm == "holub_vonka_a34":
                            n_new, _ = self._calc_moles(lam, ln_N, c, a, N_s)
                            eta_a34 = holub_vonka_a34_relative_change(n_new, n)
                            diag["holub_vonka_a34"] = eta_a34
                            if eta_a34 < tol:
                                diag["converged"] = True
                                diag["converged_by"] = "holub_vonka_a34"
                        break
                eta *= 0.5
                diag["n_backtracks"] = int(diag["n_backtracks"]) + 1
                if eta < 1e-6:
                    break
            if not accepted:
                omega0 *= 0.5
                if omega0 < 1e-3:
                    break
            elif bool(diag.get("converged", False)):
                break

        if not bool(diag.get("converged", False)):
            n_check, N_check = self._calc_moles(lam, ln_N, c, a, N_s)
            if np.all(np.isfinite(n_check)) and N_check > 0.0:
                fnorm_check = float(np.linalg.norm(self._residuals(n_check, N_check, a, b, N_E)))
                diag["final_residual"] = fnorm_check
                if tol_mode_norm == "holub_vonka_a34":
                    if float(diag.get("holub_vonka_a34", np.inf)) < tol:
                        diag["converged"] = True
                        diag["converged_by"] = "holub_vonka_a34"
                    elif fnorm_check < tol:
                        diag["converged"] = True
                        diag["converged_by"] = "holub_equilibrium_residual"
                elif fnorm_check < tol:
                    diag["converged"] = True
                    diag["converged_by"] = "l2_residual"

        n_final, _ = self._calc_moles(lam, ln_N, c, a, N_s)
        result = {species: float(n_final[j]) for j, species in enumerate(candidates)}
        diag["lambda"] = np.array(lam, dtype=np.float64)
        diag["ln_N"] = float(ln_N)
        if return_diag:
            return result, diag
        return result

    def _calc_moles(
        self, lam: np.ndarray, ln_N: float, c: np.ndarray, a: np.ndarray, N_s: int
    ) -> tuple[np.ndarray, float]:
        """由 Lagrange 乘子计算各组分摩尔量。

        Hamel (1999) (A.15)–(A.18) / 增广形式：n_j = N·exp(Σ_i a_ij λ_i − c_j)，c_j=μ°/(RT)+ln(P/P0)。
        """
        ln_N_clip = np.clip(ln_N, -50.0, 50.0)
        N = np.exp(ln_N_clip)
        exponents = np.clip(a.T @ lam - c, -50.0, 50.0)
        n = N * np.exp(exponents)
        return n, float(N)

    def _residuals(
        self, n: np.ndarray, N: float, a: np.ndarray, b: np.ndarray, N_E: int
    ) -> np.ndarray:
        """增广残差：元素守恒 (A.18) + 总摩尔归一 (A.20)。"""
        f = np.zeros(N_E + 1)
        f[:N_E] = a @ n - b
        f[N_E] = np.sum(n) / N - 1.0
        return f

    def _jacobian(
        self, n: np.ndarray, N: float, a: np.ndarray, b: np.ndarray, N_E: int
    ) -> np.ndarray:
        """增广 Jacobian，对应 Newton 步 (A.27)–(A.28)。"""
        J = np.zeros((N_E + 1, N_E + 1))
        J[:N_E, :N_E] = a @ (n[:, None] * a.T)
        J[:N_E, N_E] = b
        J[N_E, :N_E] = (a @ n) / N
        J[N_E, N_E] = 0.0
        return J

    def _initial_guess_taylor_linearized(
        self, a: np.ndarray, b: np.ndarray, c: np.ndarray, N_E: int
    ) -> np.ndarray:
        """Taylor 线性化 n_j 初值 + 守恒 lstsq（工程默认，非 (A.40) 分步）。"""
        exp_neg_c = np.exp(np.clip(-c, -50.0, 50.0))
        N0 = max(float(np.sum(np.abs(b))), 1e-20)
        n0 = N0 * exp_neg_c / (np.sum(exp_neg_c) + 1e-30)
        lhs = a @ (n0[:, None] * a.T) + 1e-10 * np.eye(N_E)
        rhs = b - a @ n0
        try:
            return np.linalg.lstsq(lhs, rhs, rcond=None)[0]
        except np.linalg.LinAlgError:
            return np.zeros(N_E)

    def _initial_guess(
        self, a: np.ndarray, b: np.ndarray, c: np.ndarray, N_E: int
    ) -> np.ndarray:
        """向后兼容别名：``taylor_linearized`` 初值。"""
        return self._initial_guess_taylor_linearized(a, b, c, N_E)

    def _build_atom_matrix(self, elem_names: List[str], candidates: List[str]) -> np.ndarray:
        """原子数矩阵 a[i,j]。"""
        a = np.zeros((len(elem_names), len(candidates)))
        for i, elem in enumerate(elem_names):
            for j, species in enumerate(candidates):
                a[i, j] = self.species_db.get_atom_count(species, elem)
        return a
