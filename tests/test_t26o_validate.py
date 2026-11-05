"""T26o 验证的判据判定：O4 须同时满足 Corella 范围拒绝与最底床格 900 °C（PR #152 审查 P2）。"""
from tools.t26o_validate import _is_o4

_REASON = (
    "雅可比不能求解：差分点不可评价：KineticsInputError: "
    "Corella 参数只覆盖 600–900 °C，ADR 0010 第 9 条没有批准外推"
)


def test_o4_requires_the_corella_rejection_and_the_boundary_temperature():
    out = _is_o4({"status": "NOT_CONVERGED", "reason": _REASON}, [790.0, 899.999998])
    assert out["ok"]
    assert out["corella_rejection"] and out["bottom_within_tolerance"]
    assert out["difference_point"]


def test_generic_jacobian_failure_is_not_o4():
    """任何差分点不可评价都带“雅可比不能求解”前缀，不能据此判 O4。"""
    out = _is_o4({"reason": "雅可比不能求解：Singular matrix"}, [790.0, 900.0])
    assert not out["ok"]
    assert not out["corella_rejection"]


def test_corella_rejection_away_from_the_boundary_is_not_o4():
    out = _is_o4({"reason": _REASON}, [790.0, 870.0])
    assert not out["ok"]
    assert out["corella_rejection"] and not out["bottom_within_tolerance"]


def test_missing_temperature_is_not_o4():
    assert not _is_o4({"reason": _REASON}, [])["ok"]
