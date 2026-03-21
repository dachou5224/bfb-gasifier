"""Gibbs 自由焓最小化求解器。

Hamel 附录 A1（pp. 163–168）/ TechSpec §5.6。
拉格朗日乘子法 + 带阻尼 Newton-Raphson 迭代。
"""

from __future__ import annotations

from typing import Dict, List, Protocol

import warnings

import numpy as np
from scipy.linalg import solve

P0 = 101325.0  # [Pa] 标准压力
R = 8.314     # [J/(mol·K)] 气体常数


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
    ) -> Dict[str, float]:
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
            收敛容差、最大迭代次数、阻尼因子

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
        c = mu0 / (R * T) + np.log(P / P0)

        lam = self._initial_guess(a, b, c, N_E)
        b_sum = max(float(np.sum(np.abs(b))), 1e-20)
        ln_N = np.log(b_sum)

        for _ in range(max_iter):
            n, N = self._calc_moles(lam, ln_N, c, a, N_s)
            if not np.all(np.isfinite(n)) or N <= 0 or not np.isfinite(N):
                omega *= 0.5
                if omega < 0.01:
                    break
                continue
            f = self._residuals(n, N, a, b, N_E)
            if np.linalg.norm(f) < tol:
                break
            J = self._jacobian(n, N, a, b, N_E)
            if not np.all(np.isfinite(J)):
                omega *= 0.5
                continue
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    dx = solve(J, -f)
            except np.linalg.LinAlgError:
                omega *= 0.5
                continue
            lam_new = lam + omega * dx[:N_E]
            ln_N_new = ln_N + omega * dx[N_E]
            ln_N_new = np.clip(ln_N_new, -50.0, 50.0)
            lam, ln_N = lam_new, ln_N_new

        n_final, _ = self._calc_moles(lam, ln_N, c, a, N_s)
        return {species: float(n_final[j]) for j, species in enumerate(candidates)}

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
