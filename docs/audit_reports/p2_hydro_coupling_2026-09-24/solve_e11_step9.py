"""从已分块通过的 E11 初值，按当前两区焓流重解。不写入 E8 正式快照。"""

from pathlib import Path

import numpy as np
from dataclasses import replace

from src.solvers.wirsum_1997_named_f import coupled_residual_named, outer_P_named
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_e8 import _e8_log_scaled_correct, freeze_e8_q

HERE = Path(__file__).resolve().parent
x = np.load(HERE / "bfb_p2_e11_step9.npz")["x"]
nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
blocks = {"passed": False}
f = np.zeros_like(x)
q = None
for _ in range(4):
    q = outer_P_named(x, nq)
    x = _e8_log_scaled_correct(q, x)
    q = outer_P_named(x, nq)
    f = coupled_residual_named(x, q, rate_scale=1.0)
    blocks = named_endpoint_blocks(f)
    print(blocks, "rows", f[-2:])
    if bool(blocks["passed"]):
        break

out = HERE / "bfb_p3_e11.npz"
np.savez(
    out,
    x=x,
    F=f,
    passed=bool(blocks["passed"]),
    gas_l2=float(blocks["gas_l2"]),
    solid_l2=float(blocks["solid_l2"]),
    tar_l2=float(blocks["tar_l2"]),
    energy_W=float(blocks["energy_W"]),
    M_fb=x[46:48],
    T_bed=float(x[-2]),
    T_fb=float(x[-1]),
    tau_s=np.asarray(q.hydro.tau_s, dtype=np.float64),
    w_return=np.asarray(q.hydro.w_return, dtype=np.float64),
    w_leave=np.asarray(q.hydro.w_leave, dtype=np.float64),
)
print("M_fb", x[46:48], "T", float(x[-2]), float(x[-1]))
print("wrote", out)
