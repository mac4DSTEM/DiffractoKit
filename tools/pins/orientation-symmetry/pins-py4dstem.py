"""Called by pins.py under the same memory cap; execute the pinned plot_orientation_maps itself.

Uses py4DSTEM's auto-sector table/parser via orientation_plan with a tiny no-correlation plan,
and the actual pymatgen PointGroup operations converted to the plan's Cartesian frame.
No implementation of direction-family reduction is proposed for the kit.
"""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import types
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import py4DSTEM
import scipy
from py4DSTEM.process.diffraction import Crystal
from py4DSTEM.process.diffraction.crystal_ACOM import orientation_ranges
from py4DSTEM.process.diffraction.utils import OrientationMap, Orientation
import py4DSTEM.process.diffraction.utils as utils
import py4DSTEM.process.diffraction.crystal as crystal_module
import py4DSTEM.process.diffraction.crystal_viz as viz
import py4DSTEM.process.diffraction.crystal_ACOM as acom

KIT=Path(__file__).resolve().parents[3]
SUITE=Path(os.environ.get('MACEM_SUITE_ROOT',KIT.parent.parent if KIT.parent.name=='_worktrees' else KIT.parent))
LOCK='f050d207851f8accd3aa550365baaf1f38230428'
clone=SUITE/'References/CuboKit/py4DSTEM-f050d207'
assert subprocess.check_output(['git','-C',str(clone),'rev-parse','HEAD'],text=True).strip()==LOCK
hashes={}
for module in [viz,acom,utils,crystal_module]:
    p=Path(module.__file__)
    assert p.read_bytes()==(clone/'py4DSTEM/process/diffraction'/p.name).read_bytes()
    hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()

assert np.__version__ == '2.5.4'
assert scipy.__version__ == '1.18.1'
assert py4DSTEM.__version__ == '0.14.19'
req=json.load(sys.stdin)
records=[]
for symbol in req['symbols']:
    zone,axis,angles=orientation_ranges[symbol]
    # Use py4DSTEM's plan parser with the class's documented automatic range, without needing
    # artificial atoms to fool a space-group detector into every one of the 11 classes.
    hexagonal=symbol in ['-3','-3m','6/m','6/mmm']
    cr=Crystal([[0,0,0]],[14], [1,1,1,90,90,120] if hexagonal else [1,1,1,90,90,90])
    cr.pymatgen_available=False
    with contextlib.redirect_stdout(io.StringIO()):
        cr.orientation_plan(zone_axis_range=np.array(zone) if isinstance(zone,list) else zone,
            fiber_axis=axis,fiber_angles=angles,angle_step_zone_axis=45,angle_step_in_plane=90,
            calculate_correlation_array=False,progress_bar=False)
    # A plot-only probe: raw directions for the monoclinic NNLS branch, and arbitrary family
    # vectors for other classes. Families are intentionally not required to form a rotation.
    system='monoclinic' if symbol=='2/m' else ('hexagonal' if hexagonal else 'cubic')
    cr.pointgroup=types.SimpleNamespace(get_crystal_system=lambda s=system:s)
    cr.pymatgen_available=True
    ops=np.array(next(g['operators'] for g in req['groups'] if g['symbol']==symbol))
    cr.symmetry_operators=np.concatenate((ops,-ops))
    inv=np.linalg.inv(cr.orientation_zone_axis_range)
    cr.symmetry_reduction=np.array([(inv.T@op).T for op in cr.symmetry_operators])
    if cr.orientation_fiber_angles is not None and abs(cr.orientation_fiber_angles[0]-180)<1e-3:
        flipped=cr.orientation_zone_axis_range.copy();flipped[0]*=-1
        inv_flip=np.linalg.inv(flipped)
        cr.symmetry_reduction=np.concatenate((cr.symmetry_reduction,
            np.array([(inv_flip.T@op).T for op in cr.symmetry_operators])))
        cr.symmetry_operators=np.tile(cr.symmetry_operators,(2,1,1))
    families=[]
    for M in req['matrices']:
        o=Orientation(1);o.matrix[0]=np.array(M)
        o=cr.symmetry_reduce_directions(o)
        families.append(o.family[0,:,2])
    ds=np.vstack((np.array(req['directions']),np.array(families)))
    # Include sector corners and positive interior directions as well as random signed inputs.
    ds=np.vstack([ds,cr.orientation_zone_axis_range,
        np.mean(cr.orientation_zone_axis_range,axis=0)])
    om=OrientationMap(1,len(ds),1)
    for i,v in enumerate(ds):
        om.matrix[0,i,0]=np.column_stack((v,v,v))
        om.family[0,i,0]=np.column_stack((v,v,v))
    with contextlib.redirect_stdout(io.StringIO()):
        rgb=cr.plot_orientation_maps(orientation_map=om,mask=np.ones((1,len(ds))),show_legend=False,
            show_axes=False,returnfig=False)[0,:,:,1]
    plt.close('all')
    float_rgb=rgb.astype(np.float32)
    records.append({'symbol':symbol,'range':cr.orientation_zone_axis_range.tolist(),
        'basis':np.linalg.inv(cr.orientation_zone_axis_range).T.tolist(),
        'directions':ds.tolist(),'rgb':rgb.tolist(),'bits':float_rgb.view(np.uint32).tolist(),
        'floatConversionGap':float(np.max(np.abs(rgb-float_rgb.astype(float))))})
json.dump({'lock':LOCK,'version':py4DSTEM.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'python':sys.version.split()[0],'hashes':hashes,'records':records},sys.stdout,allow_nan=False)
