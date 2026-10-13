"""Read-only model evaluation; writes audit evidence, never the E8 endpoint snapshot."""
from pathlib import Path
import json
import numpy as np
from src.solvers.wirsum_1997_named_f_e8 import (
    freeze_e8_q, blend_hydro_named_q, residual_e7_plus_scaled_e8,
)
from src.solvers.wirsum_1997_named_f_e7 import freeze_e7_q
from src.solvers.wirsum_1997_named_f import outer_P_named, coupled_residual_named
from src.solvers.wirsum_1997_named_f_continuation import certify_named_full_rate_endpoint
from src.solvers.wirsum_1997_named_f_endpoint_io import named_model_source_hash

art = Path(__file__).resolve().parent
q7, q8 = freeze_e7_q(), freeze_e8_q()
records = []
for filename, s in [('bfb_p2_verified0.npz', 0.), ('bfb_p2_hydro_accepted.npz', .5)]:
    z = np.load(art / filename)['x']
    h = residual_e7_plus_scaled_e8(z, q7, q8, hydro_scale=s, fate_scale=0.)
    q = blend_hydro_named_q(outer_P_named(z, q7), outer_P_named(z, q8),
                           hydro_scale=s, fate_scale=0.)
    np.testing.assert_array_equal(h, coupled_residual_named(z, q))
    cert = certify_named_full_rate_endpoint(z, q, coupling='frozen_q')
    assert cert['endpoint_certified']
    cert.pop('F')
    records.append(dict(file=filename, hydro_scale=s, fate_scale=0.,
                        T_K=float(z[-1]), x=z.tolist(), F=h.tolist(),
                        intermediate_gate=cert, full_E8_certified=False))
payload = dict(source_hash=named_model_source_hash(),
               note='Intermediate states only. Full E8 and particle fate are not certified.',
               states=records)
(art / 'fresh_verification.json').write_text(json.dumps(payload, indent=2)+'\n')
for r in records:
    print(r['hydro_scale'], r['T_K'], r['intermediate_gate'])
