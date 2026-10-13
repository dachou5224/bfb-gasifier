"""Read-only endpoint recertification and single-call timing; no solving or snapshot writes."""
import json
import time
from dataclasses import replace
import numpy as np
from src.solvers.wirsum_1997_named_f_endpoint_io import load_named_e8_snapshot, recertify_named_e8
from src.solvers.wirsum_1997_named_f import outer_P_named, coupled_residual_named
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
s = load_named_e8_snapshot()
t = time.perf_counter()
c = recertify_named_e8(np.asarray(s['e8']['x']))
c.pop('F')
e8_s = time.perf_counter() - t
x = np.load('docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p3_e11.npz')['x']
nq = replace(freeze_e8_q(Q_W=0.0), stage=11)
t = time.perf_counter()
q = outer_P_named(x, nq)
p_s = time.perf_counter() - t
t = time.perf_counter()
f = coupled_residual_named(x, q, rate_scale=1.0)
f_s = time.perf_counter() - t
print(json.dumps({'e8':c, 'e8_recertify_seconds':e8_s, 'p3':named_endpoint_blocks(f), 'p3_energy_rows_W':f[-2:].tolist(), 'single_call_seconds':{'P':p_s, 'frozen_F':f_s}, 'timing_note':'One call, not a solver benchmark.'}, indent=2))
