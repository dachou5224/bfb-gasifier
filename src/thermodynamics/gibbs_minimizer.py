"""Gibbs 自由焓最小化求解器。

Hamel 附录 A1（pp. 163–168）/ TechSpec §5.6。
拉格朗日乘子法 + 带阻尼 Newton-Raphson 迭代。
"""

from __future__ import annotations

from typing import Any, Dict, List, Protocol

import warnings

import numpy as np

from src.core.constants import Rg, P0


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

        lam = self._initial_guess(a, b, c, N_E)
        if lambda0 is not None:
            lam0 = np.asarray(lambda0, dtype=np.float64).reshape(-1)
            if lam0.size == N_E and np.all(np.isfinite(lam0)):
                lam = lam0.copy()
        b_sum = max(float(np.sum(np.abs(b))), 1e-20)
        ln_N = np.log(b_sum) if ln_N0 is None else float(np.clip(ln_N0, -50.0, 50.0))

        diag: Dict[str, Any] = {
            "n_iter": 0,
            "converged": False,
            "n_backtracks": 0,
            "eta_strategy": eta_strategy,
            "final_residual": np.inf,
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
            if fnorm < tol:
                diag["converged"] = True
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
                        break
                eta *= 0.5
                diag["n_backtracks"] = int(diag["n_backtracks"]) + 1
                if eta < 1e-6:
                    break
            if not accepted:
                omega0 *= 0.5
                if omega0 < 1e-3:
                    break

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
        """由 (A1-8) 计算各组分摩尔量。"""
        ln_N_clip = np.clip(ln_N, -50.0, 50.0)
        N = np.exp(ln_N_clip)
        exponents = np.clip(a.T @ lam - c, -50.0, 50.0)
        n = N * np.exp(exponents)
        return n, float(N)

    def _residuals(
        self, n: np.ndarray, N: float, a: np.ndarray, b: np.ndarray, N_E: int
    ) -> np.ndarray:
        """残差向量 (A1-9, A1-10)。"""
        f = np.zeros(N_E + 1)
        f[:N_E] = a @ n - b
        f[N_E] = np.sum(n) / N - 1.0
        return f

    def _jacobian(
        self, n: np.ndarray, N: float, a: np.ndarray, b: np.ndarray, N_E: int
    ) -> np.ndarray:
        """Jacobian 矩阵 (A1-12, A1-13)。"""
        J = np.zeros((N_E + 1, N_E + 1))
        J[:N_E, :N_E] = a @ (n[:, None] * a.T)
        J[:N_E, N_E] = b
        J[N_E, :N_E] = (a @ n) / N
        J[N_E, N_E] = 0.0
        return J

    def _initial_guess(
        self, a: np.ndarray, b: np.ndarray, c: np.ndarray, N_E: int
    ) -> np.ndarray:
        """Taylor 展开初始猜测 (A1-14)。"""
        exp_neg_c = np.exp(np.clip(-c, -50.0, 50.0))
        N0 = max(float(np.sum(np.abs(b))), 1e-20)
        n0 = N0 * exp_neg_c / (np.sum(exp_neg_c) + 1e-30)
        lhs = a @ (n0[:, None] * a.T) + 1e-10 * np.eye(N_E)
        rhs = b - a @ n0
        try:
            return np.linalg.lstsq(lhs, rhs, rcond=None)[0]
        except np.linalg.LinAlgError:
            return np.zeros(N_E)

    def _build_atom_matrix(self, elem_names: List[str], candidates: List[str]) -> np.ndarray:
        """原子数矩阵 a[i,j]。"""
        a = np.zeros((len(elem_names), len(candidates)))
        for i, elem in enumerate(elem_names):
            for j, species in enumerate(candidates):
                a[i, j] = self.species_db.get_atom_count(species, elem)
        return a
