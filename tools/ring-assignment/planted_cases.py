#!/usr/bin/env python3
"""tools/ring-assignment/planted_cases.py: planted ring radii for LinearScale.assign (no Python reference computes the assignment).

64 cases: fcc Al, bcc Fe, diamond Si (one with the forbidden {200} and {222} spots a thick [011] crystal shows) and hcp Mg; single zone
axes, three grains and polycrystals; tight (0.3 % radius noise) and broad (3 %) rings; with and without a 9 % elliptical distortion at a
random angle. 300 positions each, seed 7. Each case is {name, radii (pixels, per position), shells (1/A, ascending), labels, dq (the
planted scale, 1/A per pixel)}. A reader of the JSON calls LinearScale.estimateAssigned and compares the scale with dq.

Usage: planted_cases.py out.json      (numpy only)
"""
import itertools, json, math, sys
import numpy as np

# ---- crystals ----
def hkls(n=5):
    return [h for h in itertools.product(range(-n,n+1),repeat=3) if h!=(0,0,0)]
def cubic(a, rule):
    g=[(np.array(h)/a,h) for h in hkls() if rule(*h)]
    return g
fcc=lambda h,k,l: (h+k)%2==0 and (k+l)%2==0
bcc=lambda h,k,l: (h+k+l)%2==0
dia=lambda h,k,l: fcc(h,k,l) and (h%2!=0 or (h+k+l)%4==0)
def hexg(a,c):
    out=[]
    for h,k,l in hkls(4):
        if l%2!=0 and (h+2*k)%3==0: continue
        gx=h/a; gy=(h+2*k)/(a*math.sqrt(3)); gz=l/c
        out.append((np.array([gx,gy,gz]),(h,k,l)))
    return out
def shells_of(g, kmax=2.5):
    L=sorted(np.linalg.norm(v) for v,_ in g if np.linalg.norm(v)<=kmax); out=[]
    for x in L:
        if not out or x>out[-1]*(1+1e-6): out.append(x)
    return out
def zone_spots(g, zone, extra=(), kmax=1.2):
    z=np.array(zone,float); z/=np.linalg.norm(z)
    ref=np.array([1,0,0.]) if abs(z[0])<0.9 else np.array([0,1,0.])
    e1=np.cross(z,ref); e1/=np.linalg.norm(e1); e2=np.cross(z,e1)
    sp=[(v@e1,v@e2) for v,_ in list(g)+list(extra) if abs(v@z)<1e-9 and np.linalg.norm(v)<=kmax]
    return sp
def make(spots_fn, dq, rng, npos=300, sigma=0.003, ell=0.0, det=None):
    phi0=rng.uniform(0,math.pi); radii=[]
    for _ in range(npos):
        sp=spots_fn(rng); p=[]
        for x,y in sp:
            r=math.hypot(x,y)/dq; phi=math.atan2(y,x)
            r*=1+0.5*ell*math.cos(2*(phi-phi0))
            r*=1+sigma*rng.standard_normal()
            if det is None or r<=det: p.append(r)
        radii.append(p)
    return radii
def cases():
    rng=np.random.default_rng(7)
    al=cubic(4.0495,fcc); fe=cubic(2.8665,bcc); si=cubic(5.4309,dia); mg=hexg(3.2094,5.2108)
    sifcc=cubic(5.4309,fcc)
    forb=[(v,h) for v,h in sifcc if sorted(map(abs,h)) in ([0,0,2],[2,2,2])]
    def poly(g,kmax):
        sh=[s for s in shells_of(g) if s<=kmax]
        def f(rng):
            out=[]
            for s in sh:
                for _ in range(rng.integers(1,4)):
                    a=rng.uniform(0,2*math.pi); out.append((s*math.cos(a),s*math.sin(a)))
            return out
        return f
    def single(g,zone,extra=(),kmax=1.2,spin=False):
        sp=zone_spots(g,zone,extra,kmax)
        def f(rng): return sp
        return f
    def multi(g,zones,kmax=1.2):
        sps=[zone_spots(g,z,(),kmax) for z in zones]
        def f(rng): return sps[rng.integers(0,len(sps))]
        return f
    C=[]
    for nm,g,dq,items in (
        ("fcc Al",al,0.027,[("[001]",single(al,(0,0,1))),("[011]",single(al,(0,1,1))),("[111]",single(al,(1,1,1))),("poly",poly(al,1.0)),("3 grains",multi(al,[(0,0,1),(0,1,1),(1,1,1)]))]),
        ("bcc Fe",fe,0.035,[("[001]",single(fe,(0,0,1))),("[011]",single(fe,(0,1,1))),("[111]",single(fe,(1,1,1))),("poly",poly(fe,1.2))]),
        ("diamond Si",si,0.02,[("[011]",single(si,(0,1,1),kmax=0.9)),("[011]+forbidden",single(si,(0,1,1),forb,kmax=0.9)),("[001]",single(si,(0,0,1),kmax=0.9)),("poly",poly(si,0.9))]),
        ("hcp Mg",mg,0.02,[("[0001]",single(mg,(0,0,1),kmax=0.9)),("[2-1-10]",single(mg,(1,0,0),kmax=0.9)),("poly",poly(mg,0.8))]),
    ):
        sh=shells_of(g); labels=[f"s{i}" for i in range(len(sh))]
        for inm,fn in items:
            for width,sig in (("tight",0.003),("broad",0.03)):
                for ell in (0.0,0.09):
                    radii=make(fn,dq,rng,sigma=sig,ell=ell)
                    C.append({"name": f"{nm} {inm} {width} ell{ell}", "radii": [[round(float(x), 4) for x in p] for p in radii], "shells": [float(x) for x in sh], "labels": labels, "dq": dq})
    return C


if __name__ == "__main__":
    with open(sys.argv[1], "w") as fh:
        json.dump(cases(), fh)
