"""流程图映射契约：静态路径可导入 + 主链路顺序固定。"""

from __future__ import annotations

import importlib


def _resolve_dotted(qualname: str) -> object:
    parts = qualname.split(".")
    for k in range(len(parts), 0, -1):
        modname = ".".join(parts[:k])
        try:
            mod = importlib.import_module(modname)
        except ImportError:
            continue
        obj: object = mod
        for name in parts[k:]:
            obj = getattr(obj, name)
        return obj
    raise ImportError(qualname)


def test_describe_hamel_flowchart_mapping_importable():
    from src.workflow.simulation_runner import describe_hamel_flowchart_mapping

    for dotted in describe_hamel_flowchart_mapping().values():
        _resolve_dotted(dotted)


def test_hamel_global_nr_pipeline_order_importable():
    from src.workflow.simulation_runner import hamel_global_nr_pipeline_order

    for dotted in hamel_global_nr_pipeline_order():
        _resolve_dotted(dotted)


def test_pipeline_order_has_core_stages():
    from src.workflow.simulation_runner import hamel_global_nr_pipeline_order

    order = hamel_global_nr_pipeline_order()
    assert any("init_precalc_step" in p for p in order)
    assert any("nr_inner_step" in p for p in order)
    assert "src.solvers.outer_loop.run_global_nr_outer_abgleich" in order
    assert "src.solvers.result_builder.finalize_global_nr_result" in order
