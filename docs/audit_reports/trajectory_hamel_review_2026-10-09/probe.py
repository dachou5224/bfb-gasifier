"""Independent force-law integration on saved hydrodynamic state; no coupled solve."""
import json,time,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from src.solvers.wesseling_two_zone import q_on_wesseling1_path
from src.solvers.wirsum_1997_named_f_p2 import bed_face_exchange,gas_velocity_and_density
from src.physics.wirsum_1997_freeboard import hamel_eq386_asymptotic_z,integrate_particle_fate_with_time,wirsum_1997_abc,particle_acceleration_cw
from src.core.constants import g
from src.solvers.named_kernel_sources import named_kernel_source_sha256
r=json.load(open('docs/audit_reports/equation_loop_2026-10-08/certified.json'));x=np.array(r['x']);b=q_on_wesseling1_path(r['s']).reference
_,_,geo,*_=bed_face_exchange(x,b,T_bed_K=x[-2]);ug,rho,_,_=gas_velocity_and_density(x[16:23],np.zeros(7),T=x[-1],P=b.P,area=b.area)
d=float(b.d_p[0]);mu=b.mu_g;rp=b.rho_p;H=b.H_freeboard;z=hamel_eq386_asymptotic_z(b.phi_s)
kmu=18*mu/(rp*d*d);kz=.75*z*rho/(rp*d);ge=g*(1-rho/rp);slip=2*ge/(kmu+np.sqrt(kmu*kmu+4*kz*ge))
def calc(u0,shift=0.,profile=False,phi=None):
 zd=z if phi is None else hamel_eq386_asymptotic_z(phi)
 def rhs(t,y):
  gas=ug+shift
  if profile:gas+=(geo.u_b-gas)*np.exp(-max(y[0],0)/(3*geo.d_b))
  rel=y[1]-gas
  return [y[1],-ge-kmu*rel-.75*zd*rho/(rp*d)*abs(rel)*rel]
 def down(t,y):return y[0]
 down.direction=-1;down.terminal=True
 def up(t,y):return y[0]-H
 up.direction=1;up.terminal=True
 start=time.monotonic();sol=solve_ivp(rhs,[0,20000],[0,u0],method='Radau',rtol=1e-9,atol=[1e-11,1e-11],events=[down,up])
 return {'u0':u0,'shift':shift,'profile':profile,'phi':phi,'fate':'returns' if len(sol.t_events[0]) else 'leaves' if len(sol.t_events[1]) else 'undecided','t':float(sol.t[-1]),'h':float(sol.y[0,-1]),'u':float(sol.y[1,-1]),'nfev':sol.nfev,'wall_s':time.monotonic()-start}
rows=[calc(v) for v in (.89,1.34,2.24,4.47)]
rows += [calc(1.34,shift=a) for a in (.0009,.0013,.002)]
# The first audit ran this long integration; omit it on routine replay.
if '--include-long-profile' in sys.argv:
 rows += [calc(1.34,profile=True)]
else:
 rows += [{'profile': True, 'method': 'analytic equilibrium only', 'h_star': float(3*geo.d_b*np.log((geo.u_b-ug)/(slip-ug))), 'long_integration_skipped': True}]
rows += [calc(1.34,phi=a) for a in (.7,.9)]
kwargs=dict(u_g=ug,rho_g=rho,rho_p=rp,d_p=d,mu_g=mu,z=z,cell_height=H,t_max=120.)
code=[]
for dt in (.01,.005):
 start=time.monotonic();ans=integrate_particle_fate_with_time(1.34,dt=dt,**kwargs);code.append({'dt':dt,'answer':ans,'wall_s':time.monotonic()-start})
gaps=[]
for v in (-.2,0.,ug-1e-6,ug+1e-6,2.):
 a,bb,c,_=wirsum_1997_abc(u_p=v,**{k:kwargs[k] for k in ('u_g','rho_g','rho_p','d_p','mu_g','z')});force=-ge-kmu*(v-ug)-kz*abs(v-ug)*(v-ug);gaps.append(abs(a*v*v+bb*v+c-force))
out={'source_sha256':named_kernel_source_sha256(),'params':{'ug':ug,'rho':rho,'rp':rp,'d':d,'mu':mu,'H':H,'ub':geo.u_b,'db':geo.d_b,'phi':b.phi_s,'z':z,'terminal_velocity':ug-slip},'independent':rows,'code_constant':code,'ABC_max_error':max(gaps)}
Path(__file__).with_name('evidence_replay.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
