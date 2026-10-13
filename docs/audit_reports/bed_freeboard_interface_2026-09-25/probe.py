"""Read-only boundary audit; do not solve reactor or change model closures."""
from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace
import json,hashlib
import numpy as np
from src.solvers.wirsum_1997_named_f import outer_P_named,edges_from_x
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_bed_freeboard import eject_bed_top
from src.solvers.wirsum_1997_named_f_p2 import k_eject_from_geo
from src.physics.wirsum_1997_freeboard import (FREEBOARD_MODEL_HAMEL_EXTENDED,calc_roof_mass_flux_eq359,calc_wake_mass_flux_engineering)
from src.core.freeboard.transport import _build_size_class_bundles
from src.core.connectivity.solid_transport import bed_upflow_coefficient
art=Path(__file__).resolve().parent
p=Path('docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e11.npz')
x=np.load(p)['x'];q=outer_P_named(x,replace(freeze_e8_q(),stage=11));sp=q.energy.species;geo=sp.geo_fb
mg=geo.rho_p*(1-geo.eps_d)*(1-geo.eps_b)*geo.area*sp.H_bed
ed=edges_from_x(x[:39],sp);rows=[]
for i,d in enumerate(sp.d_p):
 gi=replace(geo,d_p=float(d));kw=dict(rho_p=gi.rho_p,eps_d=gi.eps_d,eps_b=gi.eps_b,u_b=gi.u_b)
 roof=calc_roof_mass_flux_eq359(d_p=float(d),d_b=gi.d_b,**kw);wake=calc_wake_mass_flux_engineering(**kw)
 rows.append(dict(dp=float(d),M=float(x[14+i]),M_fraction=float(x[14+i]/mg),roof_flux=roof,wake_flux=wake,total_all_solid_flow=(roof+wake)*geo.area,K=k_eject_from_geo(gi,sp.H_bed),char_flow=float(ed['m_eject'][i]),wake_share=wake/(roof+wake),pure_char_equiv_particles_per_s=float(ed['m_eject'][i]/(geo.rho_p*np.pi/6*d**3))))
print('Mgeom',mg,'ash_M',x[39],'char',rows,flush=True)
# Top-cell K is wake only; two-cell mesh test keeps local solid concentration fixed.
mesh=[]
for dh in [1.,.5,.25]:
 cell=SimpleNamespace(eps_b=geo.eps_b,geo=SimpleNamespace(dh=dh),u_b=geo.u_b)
 k=bed_upflow_coefficient(cell,is_top_bed_cell=True)
 mesh.append(dict(dh=dh,K=k,mchar=float(x[15]*dh/sp.H_bed),up=k*float(x[15]*dh/sp.H_bed)))
# Fixed incoming stream: vary hydrodynamics and inspect normalization.
inputs=dict(d_p_classes=np.array([.0002,.0003]),m_char_classes=np.array([.001,.002]),m_ash_classes=np.array([.0001,.0002]),rho_p=geo.rho_p,eps_b=geo.eps_b,eps_d_void=geo.eps_d,d_b=geo.d_b,area=geo.area)
normal=[]
for factor in [.5,1.,2.]:
 _,_,flows,total=_build_size_class_bundles(**inputs,u_b=geo.u_b*factor)
 normal.append(dict(ub_factor=factor,flows=flows.tolist(),total=total))
print('mesh',mesh,'normalized_FB',normal,flush=True)
res=dict(state_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),M_geom=mg,M_ash=float(x[39]),T=x[-2:].tolist(),rows=rows,mesh=mesh,normalized_FB=normal,
 geometry=dict(A=geo.area,H_bed=sp.H_bed,eps_b=geo.eps_b,eps_d=geo.eps_d,ub=geo.u_b,db=geo.d_b))
(art/'evidence.json').write_text(json.dumps(res,indent=2)+'\n')
