from pathlib import Path
_ART = Path(__file__).resolve().parent
exec(open(str(_ART / 'bfb_p2_probe.py')).read().split("print('E7'")[0])
from src.solvers.wirsum_1997_named_f_e8 import residual_e7_plus_scaled_e8, blend_hydro_named_q
from src.solvers.wirsum_1997_named_f_continuation import certify_named_full_rate_endpoint
import time
z=np.load(str(_ART / 'bfb_p2_full0_custom.npz'))['x'];t=time.perf_counter()
nq7=freeze_e7_q();nq8=freeze_e8_q()
ff=residual_e7_plus_scaled_e8(z,nq7,nq8,hydro_scale=0.,fate_scale=0.)
qq=blend_hydro_named_q(outer_P_named(z,nq7),outer_P_named(z,nq8),hydro_scale=0.,fate_scale=0.)
c=certify_named_full_rate_endpoint(z,qq,coupling='frozen_q');c.pop('F')
print('actual_s0',c,'seconds',time.perf_counter()-t,flush=True)
np.savez(str(_ART / 'bfb_p2_verified0.npz'),x=z,F=ff)
