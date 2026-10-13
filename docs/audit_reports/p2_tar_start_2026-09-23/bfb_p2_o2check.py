from pathlib import Path
_ART = Path(__file__).resolve().parent
exec(open(str(_ART / 'bfb_p2_probe.py')).read().split("print('E7'")[0])
print('dense oxygen',x[5])
for h in [1e-20,1e-26,1e-29]:
 z=x.copy();z[5]+=h;j=(f(z)-F)/h
 print('h',h,'O2diag',j[5],'COcol',j[0])
