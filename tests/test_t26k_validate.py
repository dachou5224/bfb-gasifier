"""T26k 验证报告：交替迭代记录里的标量字段也要能排版。不重跑长工况。"""

from tools import t26k_validate as mod


def _payload() -> dict:
    return {
        "elapsed_s": 12.5, "status": "DOMAIN_STOP", "reason": "带出灰多于燃料带入灰",
        "rounds": 1, "iterations": 0, "hover_drift": None, "rounds_detail": None, "t26j": None,
        "initial": {"note": "停止", "loops": [
            {"stage": "预启动自由空间焦炭块（ADR 0049 第 3 步）", "volume_ratio": [0.1, 0.2],
             "char_free_kg": [[0.01, 0.0]], "char_hover_kg": [[0.02, 0.0]],
             "freeboard_residual": 1.0e-11, "bed_residual": 1.2},
            {"loop": 1, "solid_residual": 3.5e-4, "minimum_char_kg": 2.0e-5, "bed_char_kg": 3.58,
             "change": {"up": 1.0e-3}, "up_kg_s": [0.1, 0.2], "down_kg_s": [0.3, 0.0],
             "return_kg_s": [0.01, 0.0]},
        ]},
    }


def test_markdown_accepts_scalar_fields_in_alternating_records():
    """固相门、最小 m_C、床区焦炭是标量，_fmt 不能按序列展开。审查意见第 3 条。"""
    text = mod._markdown(_payload())
    assert "固相门 0.00035" in text
    assert "最小 m_C 2e-05 kg" in text
    assert "床区焦炭 3.58 kg" in text
    assert "上行界面流 kg/s：0.1, 0.2" in text
    assert "预启动自由空间焦炭块（ADR 0049 第 3 步）" in text
    assert "自由空间固相门 1.000e-11" in text
    assert "相对变化 up 1.00e-03" in text


def test_markdown_reports_the_inner_newton_when_present():
    payload = _payload()
    payload["rounds_detail"] = [{
        "round": 1, "inner": "CONVERGED", "reason": "",
        "newton": [{
            "iteration": 1, "merit": 3.0,
            "blocks": {"gas": 1.0, "solid": 0.5, "energy": 0.1},
            "trials": [{"i": 0, "step": 1.0, "status": "ok", "reason": "", "accepted": True}],
        }],
    }]
    text = mod._markdown(payload)
    assert "外层第 1 轮内层：`CONVERGED`" in text
    assert "迭代 1：归一残差 3，" in text
    assert "试探 0（1/1 步）：ok" in text


def test_markdown_marks_the_inner_run_as_absent():
    text = mod._markdown(_payload())
    assert "（未运行：初值构造未正常结束）" in text
