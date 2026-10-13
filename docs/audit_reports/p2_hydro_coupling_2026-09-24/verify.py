"""Independently verify diagnostic scale 0.962; never write the official snapshot."""
from pathlib import Path
import json
import numpy as np
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q,residual_e7_plus_scaled_e8
from src.solvers.wirsum_1997_named_f_continuation import named_endpoint_blocks
from src.solvers.wirsum_1997_named_f_endpoint_io import named_model_source_hash
p=Path(__file__).resolve().parent
d=np.load(p/'bfb_p2_0924_picard.npz');x=d['x'];s=float(d['hydro_scale'])
h=residual_e7_plus_scaled_e8(x,freeze_e7_q(),freeze_e8_q(),hydro_scale=s,fate_scale=0.)
m=named_endpoint_blocks(h)
assert m['passed'] and np.min(x[:-1])>=0 and x[-1]>0
out=dict(hydro_scale=s,fate_scale=0.,blocks=m,T_K=float(x[-1]),x=x.tolist(),F=h.tolist(),source_hash=named_model_source_hash(),full_E8_certified=False)
(p/'fresh_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(s,x[-1],m)
