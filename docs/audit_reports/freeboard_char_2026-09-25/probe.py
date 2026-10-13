"""Bounded E11 char inventory audit; no solver or production model changes."""
from pathlib import Path
from dataclasses import replace
import hashlib,json
import numpy as np
from scipy.integrate import solve_ivp
from src.core.constants import g
from src.core.species import h_solid_mass,SOLID_H_MODE_INTEGRAL
from src.core.energy_balance_defaults import H_FORMATION_GRAPHITE_CHAR_J_PER_KG
from src.solvers.wirsum_1997_named_f import (outer_P_named,edges_from_x,fb_char_balance,
    coupled_residual_named,energy_zones_named,r234_dense_source,q_at_T)
from src.solvers.wirsum_1997_named_f_e8 import freeze_e8_q
from src.solvers.wirsum_1997_named_f_p2 import finite_hydro_from_x,k_eject_from_geo,eq372_launch_nodes
from src.solvers.wirsum_1997_coupled_f_two_bin_pyro import extents_from_x
from src.solvers.wirsum_1997_named_f_endpoint_io import named_model_source_hash
from src.physics.wirsum_1997_freeboard import drag_z_for_spec,resolve_freeboard_model,FREEBOARD_MODEL_HAMEL_EXTENDED
art=Path(__file__).resolve().parent
path=Path('docs/audit_reports/p2_hydro_coupling_2026-09-24/bfb_p2_e11.npz')
x=np.load(path)['x'];nq=replace(freeze_e8_q(),stage=11);q=outer_P_named(x,nq);sp=q.energy.species;b=nq.reference
ed=edges_from_x(x[:39],sp);f,src,c=fb_char_balance(x,q,rate_scale=1.)
qt=q_at_T(sp,float(x[-2]));ex=extents_from_x(x[:39],qt)
_,c234=r234_dense_source(x[:39],qt,rate_scale=1.,enabled=nq.enabled_rxns)
result=dict(source_hash=named_model_source_hash(),state_hash=hashlib.sha256(path.read_bytes()).hexdigest(),
    T=x[-2:].tolist(),d_p=b.d_p.tolist(),rho_p=b.rho_p,bed_M=x[14:16].tolist(),fb_M=x[46:48].tolist(),
    feed=sp.m_feed.tolist(),m_eject=ed['m_eject'].tolist(),consume_fb=c.tolist(),fb_residual=f.tolist(),
    migration_R1=float(ex['m_migrate_kg_s']),consume_R1=ex['consume_kg_s'].tolist(),consume_R234=c234.tolist(),
    Ktau=(sp.K_eject*q.hydro.tau_s).tolist())
rows=[]
for label,T,n,dt in [('current',x[-1],11,.002),('bed_T_only',x[-2],11,.002),('refine',x[-1],81,.0005)]:
    z=x.copy();z[-1]=T
    snap,geo,kbd=finite_hydro_from_x(z,b,n_bins=n,dt=dt)
    row=dict(label=label,T=float(T),ug=snap.u_g,ub=snap.u_b,rhog=snap.rho_g,kbd=kbd,
        tau=snap.tau_s.tolist(),w_leave=snap.w_leave.tolist(),w_return=snap.w_return.tolist(),
        K=[k_eject_from_geo(replace(geo,d_p=float(d)),b.H_bed) for d in b.d_p])
    rows.append(row);print(row,flush=True)
result['temperature_and_refinement']=rows
snap,geo,_=finite_hydro_from_x(x,b)
drag=drag_z_for_spec(resolve_freeboard_model(FREEBOARD_MODEL_HAMEL_EXTENDED),phi_s=geo.phi_s)
nodes,weights=eq372_launch_nodes(1.53*snap.u_b,.6,81)
rows=[]
for dp in b.d_p:
    for mode in ['constant','ghost_profile_sensitivity_only']:
        times=[];heights=[];fates=[]
        for v0 in nodes:
            if v0==0:
                times.append(0.);heights.append(0.);fates.append('return');continue
            def ode(t,y):
                h,v=y
                ug=snap.u_g if mode=='constant' else snap.u_g+(snap.u_b-snap.u_g)*np.exp(-max(h,0)/(3*b.d_b))
                r=v-ug
                return [v,-18*b.mu_g/(b.rho_p*dp**2)*r-3*drag*snap.rho_g/(4*b.rho_p*dp)*r*abs(r)-g*(1-snap.rho_g/b.rho_p)]
            def down(t,y):return y[0]
            def top(t,y):return y[0]-b.H_freeboard
            down.terminal=True;down.direction=-1;top.terminal=True;top.direction=1
            sol=solve_ivp(ode,[0,5],[0,float(v0)],rtol=1e-9,atol=1e-11,max_step=.01,events=[down,top])
            times.append(sol.t[-1]);heights.append(float(max(sol.y[0])))
            fates.append('return' if len(sol.t_events[0]) else 'leave' if len(sol.t_events[1]) else 'undecided')
        row=dict(dp=float(dp),mode=mode,tau_truncated=float(weights@times),maxheight=float(max(heights)),
            fates={ff:float(sum(w for w,fate in zip(weights,fates) if fate==ff)) for ff in ['return','leave','undecided']})
        rows.append(row);print(row,flush=True)
result['independent_ode']=rows
z=x.copy();z[14]=0.;f0=coupled_residual_named(z,q)[14];z[14]=1e-6;f1=coupled_residual_named(z,q)[14];slope=(f1-f0)/1e-6
hb=h_solid_mass('char',float(x[-2]),h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,mode=SOLID_H_MODE_INTEGRAL)
hf=h_solid_mass('char',float(x[-1]),h_f=H_FORMATION_GRAPHITE_CHAR_J_PER_KG,mode=SOLID_H_MODE_INTEGRAL)
old=float(sum(ed['m_leave'])*hb-sum(ed['m_return'])*hf)
new=float(sum(ed['m_eject'])*hb-sum(ed['m_return']-c*sp.w_return)*hf)
result['fine_and_energy']=dict(fine_source=float(f0),fine_slope=float(slope),fine_linear_root=float(-f0/slope),
    old_char_cross_W=old,consistent_char_cross_W=new,delta_cross_W=new-old,
    zone_energy=list(energy_zones_named(x,q,float(x[-2]),float(x[-1]))),
    coarse_k=float(c[1]/x[47]),coarse_Da=float(c[1]/x[47]*q.hydro.tau_s[1]))
print(result['fine_and_energy'],flush=True)
(art/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
