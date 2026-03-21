"""BFB 气化炉一维模型 — Streamlit 可视化界面。

运行方式: streamlit run app.py

Source: docs/CLAUDE.md Phase 6.2
"""

from __future__ import annotations

import json
import sys
import os
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

CHARTS_DIR = Path(__file__).parent / "docs" / "charts"
DATA_DIR = Path(__file__).parent / "data"

st.set_page_config(
    page_title="BFB Gasifier 1D",
    page_icon="🔥",
    layout="wide",
    initial_sidebar_state="expanded",
)

# BFB_ModelArchitecture.html 主题色与样式
st.markdown("""
<style>
  @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600;700&display=swap');

  :root {
    --bfb-blue:     #2563eb;
    --bfb-green:    #059669;
    --bfb-amber:    #d97706;
    --bfb-purple:   #7c3aed;
    --bfb-cyan:     #0891b2;
    --bfb-rose:     #db2777;
    --bfb-bg:       #f8fafc;
    --bfb-bg2:      #f1f5f9;
    --bfb-border:   #cbd5e1;
    --bfb-muted:    #64748b;
    --bfb-text:     #1e293b;
  }

  .block-container { padding-top: 1.2rem; max-width: 1400px; }
  .stMetric {
    background: linear-gradient(135deg, #e0f2fe 0%, #f0f9ff 100%);
    border: 1px solid var(--bfb-border);
    border-radius: 8px;
    padding: 12px 16px;
    font-family: 'IBM Plex Sans', sans-serif;
  }
  div[data-testid="stMetricValue"] { font-size: 1.5rem; color: var(--bfb-blue); }
  [data-testid="stSidebar"] {
    background: linear-gradient(180deg, #f8fafc 0%, #f1f5f9 100%);
    border-right: 1px solid var(--bfb-border);
  }
  .stTabs [data-baseweb="tab-list"] {
    background: var(--bfb-bg2);
    border-radius: 6px;
    padding: 4px;
    gap: 2px;
  }
  .stTabs [data-baseweb="tab"] {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    letter-spacing: 0.05em;
    color: var(--bfb-muted);
  }
  .stTabs [aria-selected="true"] { color: var(--bfb-blue) !important; }
  .bfb-header {
    display: flex;
    align-items: flex-end;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: 16px;
    padding-bottom: 16px;
    border-bottom: 1px solid var(--bfb-border);
    margin-bottom: 20px;
  }
  .bfb-eyebrow {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    color: var(--bfb-muted);
    letter-spacing: 0.15em;
    text-transform: uppercase;
    margin-bottom: 4px;
  }
  .bfb-badge {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    padding: 4px 10px;
    border-radius: 4px;
    border: 1px solid var(--bfb-purple);
    color: var(--bfb-purple);
    background: rgba(124, 58, 237, 0.08);
  }
  .bfb-layer {
    border-radius: 6px;
    border: 1px solid;
    padding: 14px 18px;
    margin: 8px 0;
    color: var(--bfb-text);
  }
  .bfb-layer.L1 { border-color: var(--bfb-blue); background: rgba(37, 99, 235, 0.06); }
  .bfb-layer.L2 { border-color: var(--bfb-green); background: rgba(5, 150, 105, 0.06); }
  .bfb-layer.L3 { border-color: var(--bfb-purple); background: rgba(124, 58, 237, 0.06); }
  .bfb-layer.L4 { border-color: var(--bfb-amber); background: rgba(217, 119, 6, 0.06); }
  .bfb-coupling {
    background: var(--bfb-bg2);
    border: 1px solid var(--bfb-border);
    border-left: 3px solid var(--bfb-amber);
    border-radius: 6px;
    padding: 14px 18px;
    font-size: 12px;
    color: var(--bfb-text);
    margin-top: 16px;
  }
  .bfb-coupling strong { color: var(--bfb-text); }
  .bfb-coupling code {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    color: var(--bfb-amber);
    background: rgba(217, 119, 6, 0.12);
    padding: 2px 6px;
    border-radius: 3px;
  }
</style>
""", unsafe_allow_html=True)

# ── 占位值标记集合 ──────────────────────────────────────────────
_MISSING_PREFIXES = ("__MISSING__", "__FROM_PLOT__", "__CALCULATED__")


def _is_valid(v: Any) -> bool:
    """判断字段是否为有效数值（非占位符、非 None、非字符串）。"""
    if v is None:
        return False
    if isinstance(v, str):
        return not any(v.startswith(p.split()[0]) for p in _MISSING_PREFIXES)
    return True


# ═══════════════════════════════════════════════════════════════
# 辅助函数：加载验证工况数据
# ═══════════════════════════════════════════════════════════════

@st.cache_data
def _load_validation_cases() -> dict[str, dict]:
    """读取 data/validation_cases.json，返回以 CASE_* 为键的字典。

    占位值（__MISSING__ / __FROM_PLOT__ / __CALCULATED__）在读取时保留原字符串，
    在消费端通过 _is_valid() 过滤。
    """
    json_path = DATA_DIR / "validation_cases.json"
    if not json_path.exists():
        return {}
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    return {k: v for k, v in data.items() if k.startswith("CASE_")}


def _get_case_inputs(case: dict) -> dict:
    """从工况字典中提取边界输入参数，用于侧边栏自动填入。

    返回字段（与 build_sidebar 返回键一致）。缺失的字段用 None 占位。
    """
    inp = case.get("inputs", {})
    if isinstance(inp, str):
        return {}
    fuel = inp.get("fuel", {})
    if isinstance(fuel, str):
        fuel = {}
    prox = fuel.get("proximate_analysis", {}) if isinstance(fuel, dict) else {}
    ult = fuel.get("ultimate_analysis_dry_wt_pct", {}) if isinstance(fuel, dict) else {}
    oc = inp.get("operating_conditions", {}) if isinstance(inp, dict) else {}
    reactor = inp.get("reactor", {}) if isinstance(inp, dict) else {}

    def _get(d, *keys):
        for k in keys:
            v = d.get(k) if isinstance(d, dict) else None
            if _is_valid(v):
                return float(v) if isinstance(v, (int, float)) else v
        return None

    # 颗粒直径从 particle_diameter_mm.mean_range 取均值
    dp_raw = fuel.get("particle_diameter_mm", {}) if isinstance(fuel, dict) else {}
    dp = None
    if isinstance(dp_raw, dict):
        mr = dp_raw.get("mean_range")
        if isinstance(mr, list) and len(mr) == 2:
            dp = float(np.mean(mr))
    elif isinstance(dp_raw, (int, float)):
        dp = float(dp_raw)

    return {
        "H_bed": _get(reactor, "bed_height_m"),
        "D_bed": _get(reactor, "diameter_m"),
        "P_MPa": _get(oc, "pressure_MPa"),
        "fuel_feed_kg_h": _get(fuel, "feed_rate_kg_h"),
        "ER": _get(oc, "ER"),
        "moisture": _get(prox, "moisture_wt_pct"),
        "C_dry": _get(ult, "C"),
        "H_dry": _get(ult, "H"),
        "O_dry": _get(ult, "O"),
        "VM_daf": _get(prox, "VM_daf_pct"),
        "ash_dry": _get(prox, "ash_dry_wt_pct"),
        "d_p_mm": dp,
    }


def _get_case_axial_profile(case: dict) -> dict | None:
    """提取工况的轴向剖面实验数据。

    返回 {"xi": [...], "species": {"CO": [...], ...}, "T_K": [...]}
    其中 None 占位符已被过滤（保留 null 值以便 Plotly 自动断线）。
    若该工况无可用轴向数据则返回 None。
    """
    out = case.get("outputs", {})
    if isinstance(out, str):
        return None
    ap = out.get("axial_profiles", {})
    if isinstance(ap, str) or not isinstance(ap, dict):
        return None
    if not ap.get("available", False):
        return None
    if "data" in ap and isinstance(ap["data"], str):
        return None

    xi = ap.get("xi")
    if not xi:
        return None

    species_keys = {
        "CO": "CO_mol_wet",
        "H2": "H2_mol_wet",
        "CO2": "CO2_mol_wet",
        "CH4": "CH4_mol_wet",
        "H2O": "H2O_mol_wet",
        "O2": "O2_mol_wet",
    }
    sp_data = {}
    for sp, key in species_keys.items():
        arr = ap.get(key)
        if arr and any(v is not None for v in arr):
            sp_data[sp] = arr

    t_arr = ap.get("T_K")

    return {
        "xi": xi,
        "species": sp_data,
        "T_K": t_arr if t_arr and any(v is not None for v in t_arr) else None,
    }


def _get_case_exit_data(case: dict) -> dict:
    """提取工况出口实验数据（用于 parity plot 和出口组成并排）。

    返回 {"T_K": float|None, "XC_pct": float|None, "CO": float|None, ...}
    """
    out = case.get("outputs", {})
    if isinstance(out, str):
        return {}
    eg = out.get("exit_gas_dry_mol_frac", {})
    if isinstance(eg, str):
        eg = {}
    result = {
        "T_K": out.get("exit_temperature_K") if _is_valid(out.get("exit_temperature_K")) else None,
        "XC_pct": out.get("carbon_conversion_pct") if _is_valid(out.get("carbon_conversion_pct")) else None,
    }
    for sp in ("CO", "CO2", "H2", "CH4", "H2O", "O2", "N2"):
        v = eg.get(sp)
        result[sp] = v if _is_valid(v) else None
    return result


# ═══════════════════════════════════════════════════════════════
# Plotly 通用配置
# ═══════════════════════════════════════════════════════════════

_SPECIES_COLORS = {
    "O2":  "#e74c3c",
    "CO":  "#3498db",
    "CO2": "#2ecc71",
    "H2":  "#9b59b6",
    "H2O": "#1abc9c",
    "CH4": "#f39c12",
    "N2":  "#95a5a6",
    "H2S": "#d35400",
}

# 参考 result-profile.png：各组分线型与实验标记
_SPECIES_LINE_STYLE = {
    "O2":  "dash",
    "CH4": "longdashdot",
    "CO":  "solid",
    "H2":  "dot",
    "CO2": "longdash",
    "H2O": "dashdot",
}
_SPECIES_MARKER = {
    "O2":  "circle",
    "CH4": "cross",
    "CO":  "square",
    "H2":  "triangle-down",
    "CO2": "triangle-up",
    "H2O": "diamond",
}

_PLOTLY_LAYOUT_DEFAULTS = dict(
    font=dict(family="Arial, sans-serif", size=12),
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=60, r=30, t=40, b=50),
)

# 与 validation_cases.json axial_profiles 中合成气主组分一致（湿基 mol/mol）
_SYGAS_PROFILE_SPECIES = ("CO", "H2", "CO2", "CH4")


def _numeric_range_with_padding(values: list, pad_frac: float = 0.08) -> tuple[float, float] | None:
    """从数值列表得到带边界的坐标轴范围。"""
    v = [float(x) for x in values if x is not None and np.isfinite(x)]
    if not v:
        return None
    lo, hi = min(v), max(v)
    span = max(hi - lo, 1e-9)
    return lo - pad_frac * span, hi + pad_frac * span


def _build_temperature_profile_fig(
    xi_arr: np.ndarray,
    T_profile_K: list,
    exp_profile: dict | None = None,
) -> go.Figure:
    """轴向温度剖面 T(ξ)，与 validation_cases `outputs.axial_profiles.T_K` 一致（单位 K，图中显示 °C）。"""
    T_C = [float(t) - 273.15 for t in T_profile_K]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=T_C,
            y=list(xi_arr),
            mode="lines",
            name="T 模拟",
            line=dict(color="#e67e22", width=2.5),
        )
    )
    x_all = list(T_C)
    if exp_profile and exp_profile.get("T_K"):
        xi_exp = exp_profile["xi"]
        t_exp = exp_profile["T_K"]
        pts = [(float(t) - 273.15, xi_exp[i]) for i, t in enumerate(t_exp) if t is not None]
        if pts:
            x_pts, y_pts = zip(*pts)
            x_all.extend(x_pts)
            fig.add_trace(
                go.Scatter(
                    x=list(x_pts),
                    y=list(y_pts),
                    mode="markers",
                    name="T 实验 (validation_cases)",
                    marker=dict(
                        color="#e67e22",
                        symbol="circle-open",
                        size=10,
                        line=dict(width=2, color="#c0392b"),
                    ),
                )
            )
    rng = _numeric_range_with_padding(x_all)
    x_range = rng if rng else (600.0, 1050.0)
    fig.update_layout(
        title=dict(text="温度剖面 T(ξ)", font=dict(size=14)),
        xaxis_title="温度 [°C]",
        yaxis_title="无量纲高度 ξ（0=床底，1=出口）",
        yaxis=dict(range=[0, 1]),
        xaxis=dict(range=list(x_range)),
        height=440,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        **_PLOTLY_LAYOUT_DEFAULTS,
    )
    fig.update_xaxes(showgrid=True, gridcolor="#e8e8e8")
    fig.update_yaxes(showgrid=True, gridcolor="#e8e8e8")
    return fig


def _build_syngas_concentration_fig(
    xi_arr: np.ndarray,
    comp_profiles: dict[str, list],
    exp_profile: dict | None = None,
) -> go.Figure:
    """合成气主组分湿基摩尔分数沿 ξ 分布，与 validation_cases `CO_mol_wet`、`H2_mol_wet` 等一致。"""
    fig = go.Figure()
    x_all: list[float] = []
    for sp in _SYGAS_PROFILE_SPECIES:
        vals = comp_profiles.get(sp, [])
        if not vals:
            continue
        x_all.extend(float(v) for v in vals)
        dash = _SPECIES_LINE_STYLE.get(sp, "solid")
        fig.add_trace(
            go.Scatter(
                x=vals,
                y=list(xi_arr),
                mode="lines",
                name=f"{sp} 模拟",
                line=dict(
                    color=_SPECIES_COLORS.get(sp, "#333"),
                    width=2,
                    dash=dash,
                ),
                legendgroup=sp,
            ),
        )
    if exp_profile and exp_profile.get("species"):
        xi_exp = exp_profile["xi"]
        for sp in _SYGAS_PROFILE_SPECIES:
            arr = exp_profile["species"].get(sp)
            if not arr:
                continue
            pts = [(float(arr[i]), xi_exp[i]) for i in range(len(arr)) if arr[i] is not None]
            if not pts:
                continue
            x_pts, y_pts = zip(*pts)
            x_all.extend(x_pts)
            symbol = _SPECIES_MARKER.get(sp, "circle")
            fig.add_trace(
                go.Scatter(
                    x=list(x_pts),
                    y=list(y_pts),
                    mode="markers",
                    name=f"{sp} 实验",
                    marker=dict(
                        color=_SPECIES_COLORS.get(sp, "#333"),
                        symbol=symbol,
                        size=9,
                        line=dict(width=1.5, color="#333"),
                    ),
                    legendgroup=sp,
                ),
            )
    hi = max(x_all) if x_all else 0.2
    x_max = max(0.12, hi * 1.12)
    fig.update_layout(
        title=dict(text="合成气主组分剖面（湿基 mol/mol）", font=dict(size=14)),
        xaxis_title="摩尔分数 y_j（湿基，与 CO_mol_wet、H2_mol_wet 等字段一致）",
        yaxis_title="无量纲高度 ξ",
        yaxis=dict(range=[0, 1]),
        xaxis=dict(range=[0, x_max]),
        height=440,
        legend=dict(orientation="v", yanchor="top", y=1, xanchor="left", x=1.02),
        **_PLOTLY_LAYOUT_DEFAULTS,
    )
    fig.update_xaxes(showgrid=True, gridcolor="#e8e8e8")
    fig.update_yaxes(showgrid=True, gridcolor="#e8e8e8")
    return fig


# ═══════════════════════════════════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════════════════════════════════

def build_sidebar():
    cases = _load_validation_cases()

    # 工况标签（用于 selectbox）
    case_labels = {"自定义": None}
    label_map = {
        "CASE_HTW_WESSELING_1": "HTW Wesseling #1 (Air/Steam, 2.5MPa)",
        "CASE_HTW_WESSELING_2": "HTW Wesseling #2 (O2/Steam, 2.5MPa)",
        "CASE_HTW_BERRENRATH": "HTW Berrenrath (O2/Steam, 1.0MPa)",
        "CASE_VTT_PRESSURIZED_PEAT_6": "VTT Peat #6 (Steam, 0.5MPa)",
        "CASE_VTT_PRESSURIZED_PEAT_13": "VTT Peat #13 (O2/Steam, 0.5MPa)",
        "CASE_VTT_PRESSURIZED_SAWDUST_14": "VTT Sawdust #14 (O2/Steam, 0.4MPa)",
        "CASE_VTT_ATMOSPHERIC_PEAT_D": "VTT Peat D (Atm.)",
        "CASE_VTT_ATMOSPHERIC_PEAT_E": "VTT Peat E (Atm.)",
        "CASE_WSV400_MAGDEBURG": "WSV400 Magdeburg (Wood, Atm.)",
    }
    for key in cases:
        case_labels[label_map.get(key, key)] = key

    with st.sidebar:
        st.header("操作参数")

        # ── 工况预设选择器 ────────────────────────────────────────
        st.subheader("工况预设")
        selected_label = st.selectbox(
            "选择验证工况",
            list(case_labels.keys()),
            index=0,
            help="选择文献验证工况后，以下参数将自动填入。选择「自定义」可手动设定。",
        )
        selected_case_key = case_labels[selected_label]
        preset = {}
        if selected_case_key and selected_case_key in cases:
            preset = _get_case_inputs(cases[selected_case_key])
            st.info(f"已加载: {cases[selected_case_key].get('plant', selected_case_key)}")

        def pv(field, default):
            """从 preset 取值，缺失时用 default。"""
            v = preset.get(field)
            return float(v) if v is not None else default

        _k = selected_case_key or "custom"

        # ── 反应器参数 ───────────────────────────────────────────
        st.subheader("反应器")
        n_cells = st.slider("轴向 Cell 数", 3, 30, 15, help="将床层离散为串联计算单元")
        H_bed = st.number_input("床层高度 H [m]", 1.0, 20.0, pv("H_bed", 14.5), 0.5, key=f"H_bed_{_k}")
        D_bed = st.number_input("床层直径 D [m]", 0.1, 3.0, pv("D_bed", 0.6), 0.1, key=f"D_bed_{_k}")
        
        heat_loss = st.slider(
            "热损失比例 heat_loss_frac",
            0.0, 0.5, 0.15, 0.01,
            help="基于入流总焓的损失比例。中试炉通常为 0.15-0.25",
            key=f"hlf_{_k}"
        )
        
        recirc = st.slider(
            "固体循环比例 recirculation_frac",
            0.0, 0.5, 0.1, 0.05,
            help="顶部出口固体返回底部的比例",
            key=f"recirc_{_k}"
        )

        # ── 工艺条件 ─────────────────────────────────────────────
        st.subheader("工艺条件")
        P_MPa = st.number_input("操作压力 [MPa]", 0.1, 5.0, pv("P_MPa", 2.5), 0.1, key=f"P_MPa_{_k}")
        fuel_feed_kg_h = st.number_input("燃料进料 [kg/h]", 100.0, 50000.0, pv("fuel_feed_kg_h", 3377.0), 100.0, key=f"fuel_{_k}")
        ER = st.number_input("当量比 ER [-]", 0.1, 0.6, pv("ER", 0.337), 0.01, key=f"ER_{_k}")

        # ── 燃料工业分析 ─────────────────────────────────────────
        st.subheader("燃料工业分析")
        moisture = st.number_input("含水率 M_ar [wt%]", 0.0, 40.0, pv("moisture", 16.9), 0.5, key=f"moisture_{_k}")
        C_dry = st.number_input("碳 C_d [wt%]", 30.0, 90.0, pv("C_dry", 61.5), 0.5, key=f"C_dry_{_k}")
        H_dry = st.number_input("氢 H_d [wt%]", 1.0, 10.0, pv("H_dry", 4.1), 0.1, key=f"H_dry_{_k}")
        O_dry = st.number_input("氧 O_d [wt%]", 5.0, 40.0, pv("O_dry", 21.8), 0.5, key=f"O_dry_{_k}")
        VM_daf = st.number_input("挥发分 VM_daf [wt%]", 20.0, 70.0, pv("VM_daf", 53.42), 0.5, key=f"VM_daf_{_k}")
        ash_dry = st.number_input("灰分 A_d [wt%]", 1.0, 30.0, pv("ash_dry", 11.41), 0.5, key=f"ash_dry_{_k}")
        S_dry = st.number_input(
            "硫 S_d [wt%] 干基（用于化学计量氧）",
            0.0,
            5.0,
            pv("S_dry", 0.0),
            0.05,
            key=f"S_dry_{_k}",
            help="用于 ER→O2/H2O/N2 自动计算；CASE_LU 约 0.51%",
        )
        primary_agent = st.selectbox(
            "气化剂类型",
            ["air_steam", "o2_steam", "air"],
            index=0,
            key=f"agent_{_k}",
            help="air_steam：空气+蒸汽（N2 按 79/21 配平）；o2_steam：纯氧+蒸汽",
        )
        steam_to_o2 = st.number_input(
            "蒸汽/氧 摩尔比 H2O/O2 [-]",
            0.0,
            3.0,
            pv("steam_to_o2", 0.8),
            0.05,
            key=f"sto2_{_k}",
        )

        # ── 颗粒参数 ─────────────────────────────────────────────
        st.subheader("颗粒")
        d_p_mm = st.number_input("Sauter 直径 d_p [mm]", 0.1, 5.0, pv("d_p_mm", 2.0), 0.1, key=f"d_p_{_k}")
        rho_s = st.number_input("颗粒密度 ρ_s [kg/m³]", 500, 3000, 1000, 50, key=f"rho_s_{_k}")

        # ── 求解参数 ─────────────────────────────────────────────
        st.subheader("求解")
        solver_type = st.selectbox(
            "求解器类型",
            ["gauss_seidel", "global_nr"],
            index=0,
            help="gauss_seidel: 稳健的轴向扫描；global_nr: 高精度 Newton-Raphson",
            key=f"solver_{_k}"
        )
        max_iter = st.slider("迭代步数上限", 5, 100, 30)

        run = st.button("运行模型", type="primary", use_container_width=True)

    return dict(
        n_cells=n_cells, H_bed=H_bed, D_bed=D_bed, P_MPa=P_MPa,
        fuel_feed_kg_h=fuel_feed_kg_h, ER=ER, moisture=moisture,
        C_dry=C_dry, H_dry=H_dry, O_dry=O_dry, VM_daf=VM_daf,
        ash_dry=ash_dry, S_dry=S_dry, primary_agent=primary_agent,
        steam_to_o2=steam_to_o2, d_p_mm=d_p_mm, rho_s=rho_s,
        heat_loss_frac=heat_loss, recirculation_frac=recirc,
        solver_type=solver_type,
        max_iter=max_iter, run=run,
        selected_case_key=selected_case_key,
    )


# ═══════════════════════════════════════════════════════════════
# TAB: 模型架构（BFB_ModelArchitecture 风格）
# ═══════════════════════════════════════════════════════════════

def tab_architecture():
    """求解器层级、两相 Cell、流程、模块图 — 参考 BFB_ModelArchitecture.html"""
    arch_tab = st.radio(
        "视图",
        ["求解器层级", "两相 Cell 结构", "逐步计算流程", "Python 模块图"],
        horizontal=True,
        label_visibility="collapsed",
    )

    if arch_tab == "求解器层级":
        st.markdown("#### L1–L4 求解器层级")
        st.markdown("""
<div class="bfb-layer L1">
  <div class="bfb-eyebrow">L1 外层 / Outer Loop</div>
  <strong>Newton-Raphson 全局迭代</strong><br>
  <span style="font-size:12px;color:var(--bfb-muted)">收敛判据：所有 cell 的摩尔守恒残差 ‖f‖ &lt; ε_tol</span>
</div>
<div class="bfb-layer L2">
  <div class="bfb-eyebrow">L2 内层 · 逐 cell</div>
  <strong>动力学计算</strong><br>
  <span style="font-size:12px;color:var(--bfb-muted)">流体力学 → 相间交换 → 异相 R1–R4 → 均相 R5–R11</span>
</div>
<div class="bfb-layer L3">
  <div class="bfb-eyebrow">L3 子程序 · 每次迭代调用</div>
  <strong>Gibbs 自由焓最小化（附录 A1）</strong><br>
  <span style="font-size:12px;color:var(--bfb-muted)">K_eq(T,P) 驱动力修正；(H₂S, NH₃, COS) 平衡分布</span>
</div>
<div class="bfb-layer L4">
  <div class="bfb-eyebrow">L4 输出</div>
  <strong>产品气组成 · T(h) · X_c · ΔP</strong>
</div>
<div class="bfb-coupling">
  <strong>核心耦合机制：</strong> 动力学决定碳转化"走多快"；Gibbs 决定气体"最终形态"。<br>
  耦合方程：<code>R_net = R_kinetic × (1 − Qp/Keq)</code>
</div>
        """, unsafe_allow_html=True)
        st.markdown("---")
        st.markdown("**Mermaid 流程图**")
        st.markdown("""
```mermaid
flowchart TB
    L1[L1 Newton-Raphson 全局迭代] --> L2[L2 动力学计算]
    L2 --> L3[L3 Gibbs 自由焓最小化]
    L3 --> L4[L4 输出产品气]
```
        """)

    elif arch_tab == "两相 Cell 结构":
        st.markdown("#### 两相 Cell 结构")
        st.markdown("""
<div style="display:grid;grid-template-columns:1fr auto 1fr;gap:0;max-width:700px;margin:0 auto;border:1px solid var(--bfb-border);border-radius:8px;overflow:hidden;">
  <div style="padding:16px;background:rgba(37,99,235,0.06);border-right:1px solid var(--bfb-border);color:var(--bfb-text);">
    <strong style="color:var(--bfb-blue)">🫧 气泡相</strong><br>
    <span style="font-size:12px">纯气体 · 无固体 · 均相 R5–R11 · u_b 快</span>
  </div>
  <div style="padding:16px;background:var(--bfb-bg2);display:flex;flex-direction:column;align-items:center;justify-content:center;min-width:100px;color:var(--bfb-text);">
    <span style="font-size:9px;color:var(--bfb-amber)">相间传质</span>
    <code style="font-size:10px">Ṅ_ex,bd = K_bd·V_b·(C_b−C_d)</code>
  </div>
  <div style="padding:16px;background:rgba(5,150,105,0.06);color:var(--bfb-text);">
    <strong style="color:var(--bfb-green)">🌫️ 悬浮相</strong><br>
    <span style="font-size:12px">气体+固体 · 异相+均相 R1–R11 · u_mf 最小流化</span>
  </div>
</div>
        """, unsafe_allow_html=True)
        st.markdown("---")
        st.markdown("**Mermaid 流程图**")
        st.markdown("""
```mermaid
flowchart LR
    subgraph B["🫧 气泡相"]
        B1[纯气体 · 均相 R5-R11]
    end
    subgraph E["相间传质"]
        E1["K_bd·V_b·(C_b−C_d)"]
    end
    subgraph D["🌫️ 悬浮相"]
        D1[气体+固体 · R1-R11]
    end
    B <--> E <--> D
```
        """)

    elif arch_tab == "逐步计算流程":
        st.markdown("#### 逐步计算流程")
        st.markdown("""
```mermaid
flowchart TD
    S1[1. 初始化] --> S2[2. 流体力学]
    S2 --> S3[3. 相间气体交换]
    S3 --> S4[4. 异相反应 R1-R4]
    S4 --> S5[5. 均相反应 R5-R11]
    S5 --> S6[6. Gibbs 微量组分]
    S6 --> S7[7. 能量守恒+收敛]
    S7 -->|未收敛| S2
    S7 -->|收敛| OUT[输出]
```
        """)
        st.markdown("""
| 步骤 | 内容 |
|------|------|
| 1 | T, P，燃料元素分析，气化剂流量 |
| 2 | Ergun→u_mf, Hilligardt→u_b/d_b, Sit-Grace→K_bd |
| 3 | K_bd·V_b·(C_b−C_d) |
| 4 | 炭燃烧/气化 SPM/SCM |
| 5 | Gibbs 驱动力 (1−Qp/Keq) |
| 6 | H₂S, SO₂, COS, NH₃, HCN, NO 平衡分布 |
| 7 | 更新 T，‖f‖ &lt; ε 则输出 |
        """)

    else:  # Python 模块图
        st.markdown("#### Python 模块图")
        st.markdown("""
```mermaid
flowchart LR
    subgraph core["core"]
        C1[cell.py]
        C2[species.py]
        C3[reactor.py]
    end
    subgraph physics["physics"]
        P1[hydrodynamics]
        P2[mass_transfer]
    end
    subgraph kinetics["kinetics"]
        K1[char_reactions]
        K2[gas_reactions]
    end
    subgraph thermo["thermodynamics"]
        T1[equilibrium]
        T2[gibbs_minimizer]
        T3[minor_species]
    end
    subgraph solvers["solvers"]
        S2[cell_solver]
    end
    core --> physics
    core --> kinetics
    kinetics --> thermo
    core --> solvers
```
        """)
        st.markdown("---")
        st.caption("完整架构可视化见 [docs/BFB_ModelArchitecture.html](docs/BFB_ModelArchitecture.html)")


# ═══════════════════════════════════════════════════════════════
# TAB: 模型说明
# ═══════════════════════════════════════════════════════════════

def tab_model_description():
    st.info("💡 求解器层级、两相 Cell、计算流程、模块图见 **🏗️ 架构** 标签页。")

    st.header("模型架构")

    col_img, col_txt = st.columns([1, 1])

    with col_img:
        img = CHARTS_DIR / "model-scheme.png"
        if img.exists():
            st.image(str(img), caption="Fig. 1  Cell 模型示意图 (Hamel & Krumm 2001)")

    with col_txt:
        st.markdown("""
**两相理论（Bubble–Emulsion Model）**

每个 Cell 分为 **气泡相** 和 **悬浮相（乳化相）**：

| 相态 | 含固 | 反应类型 |
|------|------|----------|
| 气泡相 (b) | 无 | 仅均相气相反应 |
| 悬浮相 (d) | 全部固体 | 均相 + 非均相气固反应 |

**相间气体交换**（Sit & Grace，Gleichung 3.50）：

$$\\dot{N}_{ex,bd,j,i} = K_{bd,i} \\cdot V_{b,i} \\cdot (C_{j,b,i} - C_{j,d,i})$$

物质从气泡相扩散/对流进入悬浮相，驱动力为两相浓度差。

**Gibbs–动力学耦合**（v11）：$R_{net} = R_{kinetic} \\times (1 - Q_p/K_{eq})$
        """)

    st.divider()

    col_scale, col_eq = st.columns([1, 1])

    with col_scale:
        img2 = CHARTS_DIR / "different-scale-of-HTW.png"
        if img2.exists():
            st.image(str(img2), caption="Fig. 2  不同规模 HTW 气化炉对比")

    with col_eq:
        st.markdown("### 守恒方程体系")
        st.markdown("""
**气相摩尔守恒** — 对每个 Cell $i$, 组分 $j$：

悬浮相 (d)：
$$0 = \\dot{N}_{zu,d} + \\dot{N}_{rez,d} + \\dot{N}_{d,i+1} + \\dot{N}_{r,d} - \\dot{N}_{d,i} - \\dot{N}_{ex}$$

气泡相 (b)：
$$0 = \\dot{N}_{zu,b} + \\dot{N}_{rez,b} + \\dot{N}_{b,i+1} + \\dot{N}_{r,b} - \\dot{N}_{b,i} + \\dot{N}_{ex}$$

> $\\dot{N}_{ex}$ 在悬浮相为负（流出），在气泡相为正（流入）。

**全局能量守恒**（焓值包含标准生成焓）：
$$0 = \\dot{H}_{in} - \\dot{H}_{out} - \\dot{Q}_W$$
        """)

    st.divider()

    st.markdown("### 化学反应网络（R1–R11）")

    rxn_data = pd.DataFrame([
        {"编号": "R1", "反应": "C + αO₂ → CO/CO₂", "类型": "非均相 (SPM)", "Arrhenius": "k_hobbs", "来源": "Hobbs (1992)"},
        {"编号": "R2", "反应": "C + H₂O → CO + H₂", "类型": "非均相 (SCM)", "Arrhenius": "k_hobbs", "来源": "Hobbs (1992)"},
        {"编号": "R3", "反应": "C + 2H₂ → CH₄", "类型": "非均相 (SCM)", "Arrhenius": "k_hobbs", "来源": "Dutta-Wen"},
        {"编号": "R4", "反应": "C + CO₂ → 2CO", "类型": "非均相 (L-H)", "Arrhenius": "k_hobbs", "来源": "Weeda (1995)"},
        {"编号": "R5", "反应": "2CO + O₂ → 2CO₂", "类型": "均相 (两相不同)", "Arrhenius": "k_standard", "来源": "Dryer & Glassman"},
        {"编号": "R6", "反应": "CH₄ + 2O₂ → CO₂ + 2H₂O", "类型": "均相", "Arrhenius": "k_standard", "来源": "Dryer (1972)"},
        {"编号": "R7", "反应": "CH₄ + H₂O → CO + 3H₂", "类型": "均相", "Arrhenius": "k_jensen", "来源": "Jensen et al."},
        {"编号": "R8", "反应": "CO + H₂O ⇌ CO₂ + H₂", "类型": "均相 (可逆)", "Arrhenius": "k_standard", "来源": "Chen (1987)"},
        {"编号": "R9", "反应": "H₂S + 1.5O₂ → SO₂ + H₂O", "类型": "均相", "Arrhenius": "k_standard", "来源": "估值"},
        {"编号": "R10", "反应": "Tar + O₂ → CO + H₂", "类型": "均相", "Arrhenius": "k_standard", "来源": "代理组分"},
        {"编号": "R11", "反应": "Tar + H₂O → CO + H₂ + CH₄", "类型": "均相", "Arrhenius": "k_standard", "来源": "代理组分"},
    ])
    st.dataframe(rxn_data, use_container_width=True, hide_index=True)

    with st.expander("Arrhenius 速率常数形式"):
        st.markdown("""
| 形式 | 公式 | 适用 |
|------|------|------|
| `k_hobbs` | $k = k_0 \\cdot T \\cdot \\exp(-E/(R_gT))$ | R1–R4 炭反应 |
| `k_standard` | $k = k_0 \\cdot \\exp(-E/(R_gT))$ | R5, R6, R8, R9–R11 |
| `k_jensen` | $k = (A/T) \\cdot \\exp(-E_T/T)$ | R7 专用 |

**注意**：R4 Boudouard 的吸附常数 $k_b$, $k_c$ 的指数 $-E/(R_gT)$ 为正（$E<0$，代表吸附热）。
        """)

    with st.expander("气体物种列表（11 组分）"):
        sp_df = pd.DataFrame([
            {"#": i+1, "组分": sp, "分子量 [g/mol]": mw, "数据来源": src}
            for i, (sp, mw, src) in enumerate([
                ("CO", 28.01, "GRI-Mech 3.0"), ("CO₂", 44.01, "GRI-Mech 3.0"),
                ("H₂", 2.016, "GRI-Mech 3.0"), ("H₂O", 18.015, "GRI-Mech 3.0"),
                ("CH₄", 16.042, "GRI-Mech 3.0"), ("O₂", 32.0, "GRI-Mech 3.0"),
                ("N₂", 28.013, "GRI-Mech 3.0"), ("H₂S", 34.08, "Burcat 2005"),
                ("NH₃", 17.03, "GRI-Mech 3.0"),
                ("TAR1", "—", "代理 (C₆H₆/C₁₀H₈)"),
                ("TAR2", "—", "代理 (C₁₀H₈/C₁₆H₃₄)"),
            ])
        ])
        st.dataframe(sp_df, use_container_width=True, hide_index=True)


# ═══════════════════════════════════════════════════════════════
# TAB: 文献参考
# ═══════════════════════════════════════════════════════════════

def tab_reference_results():
    st.header("文献参考结果")
    st.caption("Hamel & Krumm, Powder Technology 120 (2001) 105–112")

    col1, col2 = st.columns([1, 1])

    with col1:
        img = CHARTS_DIR / "result-profile.png"
        if img.exists():
            st.image(str(img), caption="Fig. 4  轴向组成与温度剖面（实验 vs 模拟）")
        st.markdown("""
**图表解读**：
- 左侧：沿无量纲高度 ξ 的主要气体组分浓度分布
- 右侧：沿高度的温度剖面（圈=实验，虚线=模拟）
- 床层底部 O₂ 迅速消耗（燃烧区），CO/H₂ 在中上部生成
- 温度峰值出现在底部燃烧区（~1000°C），向上逐渐下降
        """)

    with col2:
        img = CHARTS_DIR / "typical-results.png"
        if img.exists():
            st.image(str(img), caption="Fig. 3  多工况校验对角线图（计算 vs 实验）")
        st.markdown("""
**验证目标** (Table 2, 工况 LU)：

| 指标 | 实验值 | 允许误差 |
|------|--------|---------|
| 出口温度 T_freeboard | ~900 °C (1173 K) | < ±10% |
| 碳转化率 X_C | ~85% | < ±10% |
| CO 摩尔分数 y_CO | ~0.25 | < ±15% |

不同标记代表不同燃料：■ 褐煤 ● 泥炭 ▲ 锯屑
        """)

    st.divider()

    st.markdown("### LU 工况输入参数（Table 2，HTW 加压炉）")
    lu_df = pd.DataFrame([
        {"参数": "反应器", "值": "HTW pressurised", "单位": "—"},
        {"参数": "燃料", "值": "Rhenish brown coal", "单位": "—"},
        {"参数": "操作压力 P", "值": "2.5", "单位": "MPa"},
        {"参数": "燃料进料", "值": "3377", "单位": "kg/h"},
        {"参数": "当量比 ER", "值": "0.337", "单位": "—"},
        {"参数": "含水率 M_ar", "值": "16.9", "单位": "wt%"},
        {"参数": "灰分 A_d", "值": "11.41", "单位": "wt%"},
        {"参数": "挥发分 VM_daf", "值": "53.42", "单位": "wt%"},
        {"参数": "碳 C_d", "值": "61.5", "单位": "wt%"},
        {"参数": "氢 H_d", "值": "4.1", "单位": "wt%"},
        {"参数": "氧 O_d", "值": "21.8", "单位": "wt%"},
    ])
    st.dataframe(lu_df, use_container_width=True, hide_index=True)


# ═══════════════════════════════════════════════════════════════
# 图表构建：轴向剖面（Plotly，ξ 为 Y 轴）
# ═══════════════════════════════════════════════════════════════

def _build_axial_profile_fig(
    xi_arr: np.ndarray,
    comp_profiles: dict[str, list],
    T_profile_K: list,
    exp_profile: dict | None = None,
) -> go.Figure:
    """构建轴向剖面图：左=反应器示意，中=气体组分，右=温度。

    参考 result-profile.png：Y 轴为无量纲高度 ξ（0 底部 → 1 顶部），
    各组分使用不同线型与实验标记区分。
    """
    fig = make_subplots(
        rows=1, cols=3,
        column_widths=[0.08, 0.55, 0.37],
        subplot_titles=["", "浓度 mol/mol (湿基)", "温度 °C"],
        shared_yaxes=True,
        horizontal_spacing=0.04,
    )

    species_to_plot = ["O2", "CH4", "CO", "H2", "CO2", "H2O"]

    # 左列：反应器示意（垂直柱，ξ 0→1）
    fig.add_trace(
        go.Scatter(
            x=[0.5, 0.5], y=[0, 1],
            mode="lines",
            line=dict(color="#7f8c8d", width=14),
            showlegend=False,
        ),
        row=1, col=1,
    )
    fig.update_xaxes(showticklabels=False, showgrid=False, range=[0, 1], row=1, col=1)
    fig.update_yaxes(title_text="无量纲高度 ξ", range=[0, 1], row=1, col=1)

    # 中列：气体组分模拟线（各组分不同线型）
    for sp in species_to_plot:
        vals = comp_profiles.get(sp, [])
        if not vals:
            continue
        dash = _SPECIES_LINE_STYLE.get(sp, "solid")
        fig.add_trace(
            go.Scatter(
                x=vals, y=list(xi_arr),
                mode="lines",
                name=f"{sp} sim.",
                line=dict(
                    color=_SPECIES_COLORS.get(sp, "#333"),
                    width=2,
                    dash=dash,
                ),
                legendgroup=sp,
            ),
            row=1, col=2,
        )

    # 中列：叠加实验数据点（各组分不同标记）
    if exp_profile and exp_profile.get("species"):
        xi_exp = exp_profile["xi"]
        for sp, arr in exp_profile["species"].items():
            pts = [(xi_exp[i], v) for i, v in enumerate(arr) if v is not None]
            if not pts:
                continue
            x_pts, y_pts = zip(*pts)
            symbol = _SPECIES_MARKER.get(sp, "circle")
            fig.add_trace(
                go.Scatter(
                    x=list(x_pts), y=list(y_pts),
                    mode="markers",
                    name=f"{sp} mes.",
                    marker=dict(
                        color=_SPECIES_COLORS.get(sp, "#333"),
                        symbol=symbol,
                        size=9,
                        line=dict(width=1.5, color="#333"),
                    ),
                    legendgroup=sp,
                    showlegend=True,
                ),
                row=1, col=2,
            )

    # 右列：温度模拟线
    T_C = [t - 273.15 for t in T_profile_K]
    fig.add_trace(
        go.Scatter(
            x=T_C, y=list(xi_arr),
            mode="lines",
            name="T sim.",
            line=dict(color="#e67e22", width=2, dash="dash"),
        ),
        row=1, col=3,
    )

    # 右列：叠加实验温度点
    if exp_profile and exp_profile.get("T_K"):
        xi_exp = exp_profile["xi"]
        t_exp = exp_profile["T_K"]
        pts = [(xi_exp[i], t - 273.15) for i, t in enumerate(t_exp) if t is not None]
        if pts:
            x_pts, y_pts = zip(*pts)
            fig.add_trace(
                go.Scatter(
                    x=list(x_pts), y=list(y_pts),
                    mode="markers",
                    name="T mes.",
                    marker=dict(
                        color="#e67e22",
                        symbol="circle-open",
                        size=9,
                        line=dict(width=1.5),
                    ),
                ),
                row=1, col=3,
            )

    fig.update_layout(
        height=520,
        legend=dict(
            orientation="v",
            x=1.02,
            y=1.0,
            tracegroupgap=2,
            itemsizing="constant",
        ),
        title=dict(text="全组分湿基剖面（示意 + 全部气体 + 温度）", font=dict(size=13)),
        **_PLOTLY_LAYOUT_DEFAULTS,
    )
    fig.update_yaxes(range=[0, 1], row=1, col=2)
    fig.update_yaxes(range=[0, 1], row=1, col=3)
    # 浓度轴：按模拟与实验数据自适应
    x2_all: list[float] = []
    for sp in species_to_plot:
        vals = comp_profiles.get(sp, [])
        if vals:
            x2_all.extend(float(v) for v in vals)
    if exp_profile and exp_profile.get("species"):
        for arr in exp_profile["species"].values():
            if arr:
                x2_all.extend(float(v) for v in arr if v is not None)
    hi2 = max(x2_all) if x2_all else 0.2
    x2_max = max(0.15, min(0.55, hi2 * 1.15))
    fig.update_xaxes(title_text="mol/mol (湿基)", range=[0, x2_max], row=1, col=2)
    # 温度轴：与 validation_cases T_K 一致，自适应 °C
    t_x_all = list(T_C)
    if exp_profile and exp_profile.get("T_K"):
        t_exp = exp_profile["T_K"]
        t_x_all.extend(float(t) - 273.15 for t in t_exp if t is not None)
    t_rng = _numeric_range_with_padding(t_x_all)
    if t_rng:
        fig.update_xaxes(title_text="温度 [°C]", range=list(t_rng), row=1, col=3)
    else:
        fig.update_xaxes(title_text="温度 [°C]", range=[600, 1050], row=1, col=3)

    for col in (1, 2, 3):
        fig.update_xaxes(showgrid=True, gridcolor="#e8e8e8", row=1, col=col)
        fig.update_yaxes(showgrid=True, gridcolor="#e8e8e8", row=1, col=col)

    return fig


# ═══════════════════════════════════════════════════════════════
# 图表构建：出口组成柱状图（模拟 vs 实验并排）
# ═══════════════════════════════════════════════════════════════

def _build_exit_gas_fig(
    sim_exit: dict[str, float],
    exp_exit: dict | None = None,
) -> go.Figure:
    """出口气体组成 Plotly bar chart（模拟值 + 实验参考值并排）。"""
    bar_sp = ["CO", "CO2", "H2", "H2O", "CH4", "O2", "N2"]
    sim_vals = [sim_exit.get(sp, 0.0) * 100 for sp in bar_sp]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        name="模拟",
        x=bar_sp,
        y=sim_vals,
        marker_color=[_SPECIES_COLORS.get(sp, "#3498db") for sp in bar_sp],
        opacity=0.85,
    ))

    if exp_exit:
        exp_vals = [
            exp_exit.get(sp) * 100 if exp_exit.get(sp) is not None else None
            for sp in bar_sp
        ]
        if any(v is not None for v in exp_vals):
            fig.add_trace(go.Bar(
                name="实验 (文献)",
                x=bar_sp,
                y=exp_vals,
                marker_color="rgba(50,50,50,0.4)",
                marker_pattern_shape="/",
            ))

    fig.update_layout(
        barmode="group",
        height=350,
        xaxis_title="组分",
        yaxis_title="摩尔分数 [mol%]",
        legend=dict(orientation="h", y=1.08),
        **_PLOTLY_LAYOUT_DEFAULTS,
    )
    fig.update_yaxes(showgrid=True, gridcolor="#e0e0e0")
    return fig


# ═══════════════════════════════════════════════════════════════
# 图表构建：流体力学剖面（ξ 为 Y 轴）
# ═══════════════════════════════════════════════════════════════

def _build_hydro_profile_fig(
    xi_arr: np.ndarray,
    db_arr: list,
    eps_b_arr: list,
    kbd_arr: list,
) -> go.Figure:
    """d_b / ε_b / K_bd 三列流体力学剖面，Y 轴为无量纲高度 ξ。"""
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=["气泡直径 d_b [m]", "气泡相分率 ε_b [−]", "传质系数 K_bd [1/s]"],
        shared_yaxes=True,
    )
    traces = [
        (db_arr, "#2980b9", "d_b"),
        (eps_b_arr, "#27ae60", "ε_b"),
        (kbd_arr, "#8e44ad", "K_bd"),
    ]
    for col_idx, (vals, color, name) in enumerate(traces, start=1):
        fig.add_trace(
            go.Scatter(
                x=vals, y=list(xi_arr),
                mode="lines+markers",
                name=name,
                line=dict(color=color, width=2),
                marker=dict(size=5, color=color),
            ),
            row=1, col=col_idx,
        )
        fig.update_xaxes(showgrid=True, gridcolor="#e0e0e0", row=1, col=col_idx)
        fig.update_yaxes(showgrid=True, gridcolor="#e0e0e0", row=1, col=col_idx)

    fig.update_yaxes(title_text="无量纲高度 ξ", range=[0, 1], row=1, col=1)
    fig.update_layout(
        height=350,
        showlegend=False,
        **_PLOTLY_LAYOUT_DEFAULTS,
    )
    return fig


# ═══════════════════════════════════════════════════════════════
# TAB: 模拟运行
# ═══════════════════════════════════════════════════════════════

def tab_run_simulation(params):
    if not params["run"]:
        st.info("在左侧面板设置参数后点击 **运行模型**")
        return

    p = params
    fuel_kg_s = p["fuel_feed_kg_h"] / 3600.0

    from src.core.reactor import Reactor, ReactorConfig
    from src.core.species import configure_tar_components_by_fuel, GAS_SPECIES

    configure_tar_components_by_fuel("coal")

    cfg = ReactorConfig(
        n_cells=p["n_cells"], H_bed=p["H_bed"], D_bed=p["D_bed"],
        P=p["P_MPa"] * 1e6, T_inlet=300.0, fuel_type="coal",
        rho_s=p["rho_s"], d_p=p["d_p_mm"] * 1e-3,
        fuel_feed=fuel_kg_s,
        ER=float(p["ER"]),
        primary_agent=str(p.get("primary_agent", "air_steam")),
        S_dry=float(p.get("S_dry", 0.0)),
        steam_to_o2_molar=float(p.get("steam_to_o2", 0.8)),
        moisture_wt=p["moisture"], C_dry=p["C_dry"], H_dry=p["H_dry"], O_dry=p["O_dry"],
        VM_daf=p["VM_daf"], ash_dry_wt=p["ash_dry"], 
        heat_loss_frac=p["heat_loss_frac"],
        recirculation_frac=p["recirculation_frac"],
    )

    with st.spinner("物理模型演化中 ..."):
        reactor = Reactor(cfg)
        
        # ── 核心标定逻辑：强制物理初始化 ──
        from src.core.species import GAS_SPECIES_INDEX
        idx = GAS_SPECIES_INDEX
        for c in reactor.cells:
            c.T = cfg.T_inlet + 800.0  # 起始温度
            c.N_d[idx["N2"]] = cfg.N2_feed
            c.N_d[idx["H2O"]] = cfg.H2O_feed
            c.N_d[idx["O2"]] = cfg.O2_feed

        result = reactor.solve(
            max_global_iter=p["max_iter"], 
            tol_global=1e-3,
            solver=p["solver_type"]
        )

    # ── KPI 面板 ─────────────────────────────────────────────────
    st.markdown("### 关键结果指标")
    T_exit = result["T_profile"][-1]
    T_exp, XC_exp, yCO_exp = 1173.15, 0.85, 0.25
    T_err = abs(T_exit - T_exp) / T_exp * 100
    XC = result["carbon_conv"]
    XC_err = abs(XC - XC_exp) / XC_exp * 100
    yCO = result["exit_gas"].get("CO", 0.0)
    CO_err = abs(yCO - yCO_exp) / yCO_exp * 100

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("出口温度", f"{T_exit:.0f} K", f"误差 {T_err:.1f}%")
    c2.metric("出口温度", f"{T_exit - 273.15:.0f} °C")
    c3.metric("碳转化率", f"{XC*100:.1f}%", f"误差 {XC_err:.1f}%")
    c4.metric("y_CO", f"{yCO:.4f}", f"误差 {CO_err:.1f}%")
    _conv = result.get("converged_fully")
    if _conv is None:
        _conv = bool(result.get("converged", False))
    _sub = "收敛" if _conv else "未收敛"
    if result.get("converged_fully") is False and result.get("converged_outer", result.get("converged")):
        _sub = "外收敛/内未收敛"
    c5.metric("迭代次数", f"{result['n_iter']}", _sub)

    st.divider()

    # ── 轴向剖面（与 validation_cases.json axial_profiles 对齐）────────────────
    st.markdown("### 轴向剖面")
    st.caption(
        "数据字段与 `data/validation_cases.json` → `outputs.axial_profiles` 一致："
        "`xi` 无量纲高度；`T_K` 温度；`CO_mol_wet`、`H2_mol_wet` 等为湿基摩尔分数。"
        " 选择工况预设时叠加文献实验点。"
    )
    n = len(result["T_profile"])
    h_arr = np.array([(i + 0.5) * cfg.H_bed / n for i in range(n)])
    xi_arr = h_arr / cfg.H_bed

    species_to_plot = ["O2", "CO", "CO2", "H2", "H2O", "CH4"]
    comp_profiles = {sp: [] for sp in species_to_plot}
    for cell in reactor.cells:
        y = cell._mole_fractions("d")
        for sp in species_to_plot:
            from src.core.species import GAS_SPECIES_INDEX
            comp_profiles[sp].append(float(y[GAS_SPECIES_INDEX[sp]]))

    exp_profile = None
    cases = _load_validation_cases()
    case_key = p.get("selected_case_key")
    if case_key and case_key in cases:
        exp_profile = _get_case_axial_profile(cases[case_key])

    col_t, col_s = st.columns(2)
    with col_t:
        temp_fig = _build_temperature_profile_fig(
            xi_arr, result["T_profile"], exp_profile
        )
        st.plotly_chart(temp_fig, use_container_width=True)
    with col_s:
        syngas_fig = _build_syngas_concentration_fig(
            xi_arr, comp_profiles, exp_profile
        )
        st.plotly_chart(syngas_fig, use_container_width=True)

    with st.expander("全组分湿基剖面（含 O₂、H₂O 与反应器示意列）", expanded=False):
        axial_fig = _build_axial_profile_fig(
            xi_arr, comp_profiles, result["T_profile"], exp_profile
        )
        st.plotly_chart(axial_fig, use_container_width=True)

    # ── 新增：炭平衡诊断 ──────────────────────────────────────────
    with st.expander("炭平衡与动力学诊断", expanded=False):
        col_char, col_o2 = st.columns(2)
        
        # 1. 炭质量流剖面
        char_mass_flow = [float(np.sum(c.m_solid[:, 0])) for c in reactor.cells]
        fig_char = go.Figure(go.Scatter(
            x=char_mass_flow, y=list(xi_arr),
            mode="lines+markers",
            line=dict(color="#7f8c8d", width=2.5),
            name="炭质量流率"
        ))
        fig_char.update_layout(
            title="轴向炭质量流率 [kg/s]",
            xaxis_title="m_char [kg/s]",
            yaxis_title="无量纲高度 ξ",
            height=350,
            **_PLOTLY_LAYOUT_DEFAULTS
        )
        col_char.plotly_chart(fig_char, use_container_width=True)
        
        # 2. 氧气消耗速率诊断
        o2_reaction_rate = [float(c.R_gas_d[idx["O2"]]) for c in reactor.cells]
        fig_o2r = go.Figure(go.Scatter(
            x=o2_reaction_rate, y=list(xi_arr),
            mode="lines+markers",
            line=dict(color="#e74c3c", width=2.5),
            name="O2 净反应速率"
        ))
        fig_o2r.update_layout(
            title="轴向 O2 净反应速率 [mol/s]",
            xaxis_title="R_gas_d[O2] [mol/s]",
            yaxis_title="无量纲高度 ξ",
            height=350,
            **_PLOTLY_LAYOUT_DEFAULTS
        )
        col_o2.plotly_chart(fig_o2r, use_container_width=True)

    st.divider()

    # ── 出口气体组成 ─────────────────────────────────────────────
    st.markdown("### 出口气体组成")
    col_bar, col_table = st.columns([1, 1])

    exp_exit = None
    if case_key and case_key in cases:
        exp_exit = _get_case_exit_data(cases[case_key])

    with col_bar:
        exit_fig = _build_exit_gas_fig(result["exit_gas"], exp_exit if exp_exit else None)
        st.plotly_chart(exit_fig, use_container_width=True)

    with col_table:
        from src.core.species import GAS_SPECIES
        table_rows = []
        for sp in GAS_SPECIES:
            y_sim = result["exit_gas"].get(sp, 0.0)
            y_exp = (exp_exit or {}).get(sp)
            row = {
                "组分": sp,
                "模拟 mol%": f"{y_sim*100:.3f}",
            }
            if y_exp is not None:
                row["实验 mol%"] = f"{y_exp*100:.3f}"
                row["误差 %"] = f"{abs(y_sim - y_exp) / max(y_exp, 1e-9) * 100:.1f}"
            table_rows.append(row)
        st.dataframe(pd.DataFrame(table_rows), use_container_width=True, hide_index=True)

    st.divider()

    # ── 流体力学参数 ─────────────────────────────────────────────
    st.markdown("### 流体力学参数")

    db_arr = [c.d_b for c in reactor.cells]
    eps_b_arr = [c.eps_b for c in reactor.cells]
    kbd_arr = [c.K_bd for c in reactor.cells]

    hydro_fig = _build_hydro_profile_fig(xi_arr, db_arr, eps_b_arr, kbd_arr)
    st.plotly_chart(hydro_fig, use_container_width=True)

    hydro_rows = []
    for i, cell in enumerate(reactor.cells):
        hydro_rows.append({
            "Cell": i + 1,
            "h [m]": f"{h_arr[i]:.2f}",
            "u_mf [m/s]": f"{cell.u_mf:.4f}",
            "u0 [m/s]": f"{cell.u0:.4f}",
            "u_b [m/s]": f"{cell.u_b:.3f}",
            "d_b [m]": f"{cell.d_b:.4f}",
            "ε_b": f"{cell.eps_b:.4f}",
            "K_bd [1/s]": f"{cell.K_bd:.2f}",
        })
    with st.expander("各 Cell 流体力学参数详表"):
        st.dataframe(pd.DataFrame(hydro_rows), use_container_width=True, hide_index=True)

    st.caption(
        "⚠ 当前结果基于文献估算参数（部分标注 TODO）。"
        "精确预测需提供 Hamel (1999) Table 5.2 / 6.3 动力学参数。"
    )

    # 保存本次运行结果到 session_state（供 parity plot 使用）
    st.session_state["last_sim_result"] = {
        "T_K": T_exit,
        "XC_pct": XC * 100,
        "CO": result["exit_gas"].get("CO"),
        "CO2": result["exit_gas"].get("CO2"),
        "H2": result["exit_gas"].get("H2"),
        "CH4": result["exit_gas"].get("CH4"),
        "H2O": result["exit_gas"].get("H2O"),
        "case_key": case_key,
    }


# ═══════════════════════════════════════════════════════════════
# 图表构建：动态 Parity Plot（6 宫格）
# ═══════════════════════════════════════════════════════════════

def _build_parity_plot(
    cases: dict[str, dict],
    sim_result: dict | None = None,
) -> go.Figure:
    """动态 Plotly 6 宫格 parity plot。

    面板：(a) T_freeboard [K], (b) X_C [%], (c) CO, (d) CH4, (e) H2, (f) H2O
    每图：背景点（全工况文献值，灰色）+ 对角线 + ±10%/±15% 误差带
    若 sim_result 不为 None，叠加当前运行结果（高亮彩色标记）。
    """
    panels = [
        ("T_freeboard", "T_freeboard [K]", "T_K", "T_K"),
        ("X_C", "碳转化率 X_C [%]", "XC_pct", "XC_pct"),
        ("CO", "CO mol/mol (dry)", "CO", "CO"),
        ("CH4", "CH₄ mol/mol (dry)", "CH4", "CH4"),
        ("H2", "H₂ mol/mol (dry)", "H2", "H2"),
        ("H2O", "H₂O mol/mol (wet)", "H2O", "H2O"),
    ]

    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=[p[1] for p in panels],
        horizontal_spacing=0.08,
        vertical_spacing=0.12,
    )

    fuel_symbols = {
        "coal": "square",
        "peat": "circle",
        "wood": "triangle-up",
        "sawdust": "triangle-up",
    }
    fuel_colors = {
        "coal": "#2c3e50",
        "peat": "#7f8c8d",
        "wood": "#7f8c8d",
        "sawdust": "#95a5a6",
    }

    def _fuel_type(case: dict) -> str:
        ft = case.get("inputs", {}).get("fuel", {}).get("type", "").lower()
        if "coal" in ft:
            return "coal"
        if "peat" in ft:
            return "peat"
        if "sawdust" in ft or "säge" in ft:
            return "sawdust"
        if "wood" in ft or "holz" in ft:
            return "wood"
        return "coal"

    for panel_idx, (panel_id, panel_title, exp_key, sim_key) in enumerate(panels):
        row = panel_idx // 3 + 1
        col = panel_idx % 3 + 1

        # 收集所有工况的实验值（X 轴）
        all_exp = []
        all_sim_bg = []  # 背景点使用实验值作 placeholder（真实模拟结果未知，用实验值散点表示覆盖度）
        # 实际上背景点 x=实验值，y 我们没有其他工况的模拟值
        # 用灰色散点仅显示实验点的范围分布，对角线表示"完美预测"
        case_exp_vals = []
        case_fuel_types = []
        for ck, cv in cases.items():
            ed = _get_case_exit_data(cv)
            if exp_key == "T_K":
                v = ed.get("T_K")
            elif exp_key == "XC_pct":
                v = ed.get("XC_pct")
            else:
                v = ed.get(exp_key)
            if v is not None:
                case_exp_vals.append(float(v))
                case_fuel_types.append(_fuel_type(cv))

        if not case_exp_vals:
            # 无数据，添加占位文本
            fig.add_trace(
                go.Scatter(
                    x=[0.5], y=[0.5],
                    mode="text",
                    text=["数据缺失"],
                    textposition="middle center",
                    showlegend=False,
                ),
                row=row, col=col,
            )
            fig.update_xaxes(range=[0, 1], row=row, col=col)
            fig.update_yaxes(range=[0, 1], row=row, col=col)
            continue

        x_min, x_max = min(case_exp_vals) * 0.85, max(case_exp_vals) * 1.15
        x_diag = np.linspace(x_min, x_max, 100)

        # 对角线（完美预测）
        fig.add_trace(
            go.Scatter(
                x=x_diag, y=x_diag,
                mode="lines",
                line=dict(color="#333333", width=1.5),
                name="完美预测",
                showlegend=(panel_idx == 0),
                legendgroup="diagonal",
            ),
            row=row, col=col,
        )

        # ±10% 误差带
        fig.add_trace(
            go.Scatter(
                x=list(x_diag) + list(x_diag[::-1]),
                y=list(x_diag * 1.10) + list(x_diag[::-1] * 0.90),
                fill="toself",
                fillcolor="rgba(100,180,100,0.15)",
                line=dict(color="rgba(100,180,100,0.4)", dash="dot"),
                name="±10%",
                showlegend=(panel_idx == 0),
                legendgroup="band10",
            ),
            row=row, col=col,
        )

        # ±15% 误差带
        fig.add_trace(
            go.Scatter(
                x=list(x_diag) + list(x_diag[::-1]),
                y=list(x_diag * 1.15) + list(x_diag[::-1] * 0.85),
                fill="toself",
                fillcolor="rgba(200,160,50,0.10)",
                line=dict(color="rgba(200,160,50,0.3)", dash="dot"),
                name="±15%",
                showlegend=(panel_idx == 0),
                legendgroup="band15",
            ),
            row=row, col=col,
        )

        # 文献实验点（灰色，按燃料类型区分符号）
        for ft in set(case_fuel_types):
            mask = [i for i, t in enumerate(case_fuel_types) if t == ft]
            x_pts = [case_exp_vals[i] for i in mask]
            fig.add_trace(
                go.Scatter(
                    x=x_pts, y=x_pts,
                    mode="markers",
                    name=ft,
                    marker=dict(
                        symbol=fuel_symbols.get(ft, "circle"),
                        color=fuel_colors.get(ft, "#999"),
                        size=8,
                        line=dict(color="#666", width=1),
                    ),
                    showlegend=(panel_idx == 0),
                    legendgroup=f"fuel_{ft}",
                ),
                row=row, col=col,
            )

        # 当前运行结果（高亮）
        if sim_result:
            if sim_key == "T_K":
                sim_val = sim_result.get("T_K")
            elif sim_key == "XC_pct":
                sim_val = sim_result.get("XC_pct")
            else:
                sim_val = sim_result.get(sim_key)

            # 对应的实验参考值（若工况已知）
            sim_case_key = sim_result.get("case_key")
            exp_val_for_sim = None
            if sim_case_key and sim_case_key in cases:
                ed = _get_case_exit_data(cases[sim_case_key])
                if exp_key == "T_K":
                    exp_val_for_sim = ed.get("T_K")
                elif exp_key == "XC_pct":
                    exp_val_for_sim = ed.get("XC_pct")
                else:
                    exp_val_for_sim = ed.get(exp_key)

            if sim_val is not None and exp_val_for_sim is not None:
                fig.add_trace(
                    go.Scatter(
                        x=[float(exp_val_for_sim)],
                        y=[float(sim_val)],
                        mode="markers",
                        name="本次运行",
                        marker=dict(
                            symbol="star",
                            color="#e74c3c",
                            size=14,
                            line=dict(color="#fff", width=1.5),
                        ),
                        showlegend=(panel_idx == 0),
                        legendgroup="current_run",
                    ),
                    row=row, col=col,
                )

        # 坐标轴格式
        fig.update_xaxes(
            title_text="实验值", showgrid=True, gridcolor="#e0e0e0",
            range=[x_min, x_max], row=row, col=col,
        )
        fig.update_yaxes(
            title_text="计算值", showgrid=True, gridcolor="#e0e0e0",
            range=[x_min, x_max], row=row, col=col,
        )

    fig.update_layout(
        height=700,
        legend=dict(orientation="v", x=1.02, y=1.0, tracegroupgap=4),
        title_text="多工况校验 Parity Plot（计算值 vs 实验值）",
        **_PLOTLY_LAYOUT_DEFAULTS,
    )
    return fig


# ═══════════════════════════════════════════════════════════════
# TAB: 模块验证
# ═══════════════════════════════════════════════════════════════

def tab_validation():
    st.header("子模块数量级验证（Sanity Checks）")

    st.markdown("""
按 `docs/CLAUDE.md` 要求，每个物理子模块必须通过以下 6 项数量级检验。
以下结果为静态计算（不依赖 Reactor 求解器）。
    """)

    checks = []

    # 1. u_mf
    from src.physics.minimum_fluidization import compute_u_mf
    from src.core.species import gas_viscosity_power_law

    P, T = 2.5e6, 1173.15
    rho_g = P * 0.029 / (8.314 * T)
    mu_g = gas_viscosity_power_law(T, 1.8e-5)
    u_mf = compute_u_mf(rho_g, 1000.0, 0.5e-3, mu_g)
    checks.append({"检验": "u_mf", "条件": "d_p=0.5mm, ρ_s=1000, T=900K, P=2.5MPa",
                    "结果": f"{u_mf:.4f} m/s", "范围": "0.02–0.08 m/s",
                    "状态": "✅" if 0.02 <= u_mf <= 0.08 else "❌"})

    # 2. d_b
    from src.physics.bubble_dynamics import integrate_bubble_diameter
    _, db = integrate_bubble_diameter(0.3, 0.05, 2.5e6, 5.0, D_bed=0.6)
    db_top = db[-1]
    checks.append({"检验": "d_b(H_bed)", "条件": "u0=0.3, u_mf=0.05, P=2.5MPa, H=5m",
                    "结果": f"{db_top:.4f} m", "范围": "0.05–0.30 m",
                    "状态": "✅" if 0.05 <= db_top <= 0.3 else "❌"})

    # 3. K_bd
    from src.physics.mass_transfer import calc_kbd, calc_u_br
    from src.core.species import gas_diffusivity_correlation
    from src.physics.bubble_dynamics import bubble_rise_velocity

    D_g = gas_diffusivity_correlation(1000.0, 2.5e6)
    u_br = calc_u_br(0.05 / 0.45, 2.5e6)
    u_b = bubble_rise_velocity(0.3, 0.05, 0.1)
    K_bd = calc_kbd(u_br, 0.1, D_g, 0.45, u_b)
    checks.append({"检验": "K_bd", "条件": "d_b=0.1m, P=2.5MPa, T=1000K",
                    "结果": f"{K_bd:.2f} 1/s", "范围": "1–15 1/s",
                    "状态": "✅" if 1 <= K_bd <= 15 else "❌"})

    # 4. R4 Boudouard
    from src.kinetics.char_reactions import rate_R4
    r4 = rate_R4(1073.0, 5e4, 1e4)
    checks.append({"检验": "R4 Boudouard", "条件": "T=1073K, P_CO2=5e4Pa",
                    "结果": f"{r4:.2e} mol/(m²·s)", "范围": "1e-7 – 1e-3",
                    "状态": "✅" if 1e-7 < r4 < 1e-3 else "❌"})

    # 5. R8 WGSR
    from src.kinetics.gas_reactions import rate_R8
    r8 = rate_R8(1073.0, 2.5e6, 0.25, 0.10, 0.05, 0.15)
    checks.append({"检验": "R8 WGSR 方向", "条件": "T=1073K, P=2.5MPa, y_CO=0.25",
                    "结果": f"{r8:.4f} mol/(m³·s)", "范围": "> 0 (正向)",
                    "状态": "✅" if r8 > 0 else "❌"})

    # 6. DAEM
    from src.thermal.devolatilization import daem_conversion, _linear_heating_profile
    T_h, t_h = _linear_heating_profile(300.0, 1173.15, 100.0, 500)
    X_VM = daem_conversion(T_h, t_h)
    checks.append({"检验": "DAEM 释放率", "条件": "300→900°C, 100s",
                    "结果": f"{X_VM*100:.1f}%", "范围": "> 80%",
                    "状态": "✅" if X_VM > 0.8 else "❌"})

    df = pd.DataFrame(checks)
    st.dataframe(df, use_container_width=True, hide_index=True)

    n_pass = sum(1 for c in checks if c["状态"] == "✅")
    if n_pass == 6:
        st.success(f"全部 {n_pass}/6 项通过")
    else:
        st.warning(f"{n_pass}/6 项通过")

    st.divider()

    # ── 流体力学演示（Plotly，ξ 为 Y 轴） ─────────────────────────
    st.markdown("### 流体力学子模型演示")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**气泡直径沿高度变化** (Mori-Wen, H=5m)")
        h_demo = np.linspace(0, 5, 100)
        xi_demo = h_demo / 5.0
        _, db_arr = integrate_bubble_diameter(0.3, 0.05, 2.5e6, 5.0, D_bed=0.6, n_points=100)
        fig_db = go.Figure(go.Scatter(
            x=list(db_arr), y=list(xi_demo),
            mode="lines",
            line=dict(color="#2980b9", width=2),
        ))
        fig_db.update_layout(
            height=320,
            xaxis_title="d_b [m]",
            yaxis_title="无量纲高度 ξ",
            yaxis=dict(range=[0, 1]),
            **_PLOTLY_LAYOUT_DEFAULTS,
        )
        fig_db.update_xaxes(showgrid=True, gridcolor="#e0e0e0")
        fig_db.update_yaxes(showgrid=True, gridcolor="#e0e0e0")
        st.plotly_chart(fig_db, use_container_width=True)

    with col2:
        st.markdown("**DAEM 挥发分释放 S 曲线**")
        temps = np.arange(300, 1001, 25)
        x_vals = []
        for T_end in temps:
            T_h2, t_h2 = _linear_heating_profile(300.0, T_end + 273.15, 100.0, 200)
            x_vals.append(daem_conversion(T_h2, t_h2) * 100)
        fig_daem = go.Figure(go.Scatter(
            x=list(temps), y=x_vals,
            mode="lines+markers",
            line=dict(color="#e67e22", width=2),
            marker=dict(size=4),
        ))
        fig_daem.update_layout(
            height=320,
            xaxis_title="终温 [°C]",
            yaxis_title="挥发分转化率 X_VM [%]",
            **_PLOTLY_LAYOUT_DEFAULTS,
        )
        fig_daem.update_xaxes(showgrid=True, gridcolor="#e0e0e0")
        fig_daem.update_yaxes(showgrid=True, gridcolor="#e0e0e0")
        st.plotly_chart(fig_daem, use_container_width=True)

    st.divider()

    # ── 动态 Parity Plot ────────────────────────────────────────
    st.markdown("### 多工况校验 Parity Plot")
    st.caption(
        "背景散点：文献实验值分布 | 红星：本次运行结果（需先在「模拟运行」标签页执行） | "
        "绿带=±10%，黄带=±15%"
    )

    cases = _load_validation_cases()
    sim_result = st.session_state.get("last_sim_result")

    if not cases:
        st.warning("未找到 data/validation_cases.json，无法显示 parity plot。")
    else:
        parity_fig = _build_parity_plot(cases, sim_result)
        st.plotly_chart(parity_fig, use_container_width=True)

        if sim_result is None:
            st.info("提示：在「🚀 模拟运行」标签页执行模型后，本次结果将以红星显示在对应工况的 parity plot 上。")


# ═══════════════════════════════════════════════════════════════
# TAB: TODO 清单
# ═══════════════════════════════════════════════════════════════

def tab_todos():
    st.header("待确认参数 / TODO 清单")
    st.markdown("以下参数使用了文献估算值，需要替换为 Hamel (1999) 实际数据。")

    todos = pd.DataFrame([
        {"优先级": "🔴 高", "模块": "char_reactions.py", "内容": "R1-R4 所有 k₀ / E 参数", "来源": "Table 5.2"},
        {"优先级": "🔴 高", "模块": "gas_reactions.py", "内容": "R5-R9 所有 k₀ / E 参数", "来源": "Table 6.3"},
        {"优先级": "🔴 高", "模块": "tar_reactions.py", "内容": "R10/R11 k₀ / E + 浓度指数", "来源": "文献"},
        {"优先级": "🟡 中", "模块": "devolatilization.py", "内容": "DAEM A / E₀ / σ 按燃料调整", "来源": "实验拟合"},
        {"优先级": "🟡 中", "模块": "drying.py", "内容": "湿颗粒物性 λ / ρ / Cp", "来源": "燃料数据"},
        {"优先级": "🟡 中", "模块": "bubble_dynamics.py", "内容": "A_bed / N_or 分布板参数", "来源": "反应器设计"},
        {"优先级": "🟢 低", "模块": "cell.py", "内容": "粒径类间迁移 m_left / m_right", "来源": "specs Eq.2-3"},
        {"优先级": "🟢 低", "模块": "cell.py", "内容": "干燥/热解与 cell 完整耦合", "来源": "Phase 4 → 5"},
        {"优先级": "🟢 低", "模块": "reactor.py", "内容": "自由板区夹带计算集成", "来源": "freeboard.py"},
        {"优先级": "🟢 低", "模块": "phase_fractions.py", "内容": "Hilligardt ODE 参数校准 (高压)", "来源": "Hilligardt (1986)"},
    ])
    st.dataframe(todos, use_container_width=True, hide_index=True)


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    # BFB_ModelArchitecture 风格 header
    st.markdown("""
    <div class="bfb-header">
      <div>
        <div class="bfb-eyebrow">BFB Gasifier 1D · Hamel & Krumm (2001)</div>
        <h1 style="font-size:1.6rem;font-weight:700;color:#1e293b;margin:0;">鼓泡流化床气化炉<br>一维稳态动力学模型</h1>
        <div class="bfb-eyebrow" style="margin-top:4px;">两相理论 | Gibbs–动力学耦合 | Python 实现</div>
      </div>
      <span class="bfb-badge">v11.0</span>
    </div>
    """, unsafe_allow_html=True)

    params = build_sidebar()

    tab_arch, tab_desc, tab_ref, tab_sim, tab_val, tab_todo = st.tabs([
        "🏗️ 架构",
        "📖 模型说明",
        "📊 文献参考",
        "🚀 模拟运行",
        "✅ 模块验证",
        "📝 TODO",
    ])

    with tab_arch:
        tab_architecture()

    with tab_desc:
        tab_model_description()

    with tab_ref:
        tab_reference_results()

    with tab_sim:
        tab_run_simulation(params)

    with tab_val:
        tab_validation()

    with tab_todo:
        tab_todos()


if __name__ == "__main__":
    main()
