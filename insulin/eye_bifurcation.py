#!/usr/bin/env python3
"""
Rigorous characterization of the local-formate feed-forward in the eye.

The autonomous-eye feed-forward is a self-consistency condition on the local
mitochondrial formate production f:

    f  =  EYE_FMITO  +  gain * mTOR(f) * fb(mTOR(f))                (*)

where mTOR(f) is the eye cell's steady-state mTOR at imposed local formate f
(computed with the ROBUST, polished beta_cell.steady_state -- not a bare
integration), and fb is the S6K1->IRS-1 negative feedback. Fixed points are the
roots of (*). One root -> monostable; three roots (two stable, one unstable) ->
BISTABLE, i.e. a genuine switch with hysteresis. This maps the roots as a
function of blood formate, so we can say whether the switch is real (bistable) or
just a steep monostable response, and locate the fold (saddle-node) points.
"""
import numpy as np
import pandas as pd
from model import BetaCellModel, EyeModel

M = BetaCellModel(); M.calibrate_all(verbose=False)
E = EyeModel(beta=M)

GLU = 5.0


def _mtor_ins_of_f(Fb, fgrid):
    """mTOR and eye insulin as functions of IMPOSED local formate f, at blood
    formate Fb, using the polished steady-state solver (reliable fixed points)."""
    mt = np.empty_like(fgrid); ins = np.empty_like(fgrid); y = None
    for i, f in enumerate(fgrid):
        r, y = M.steady_state(GLU, Fb, f_mito=float(f), y0=y, t_pre=40.0)
        mt[i] = r['mtor']; ins[i] = r['v_ins_syn']
    return mt, ins


def fixed_points(Fb, fgrid, gain=None, ki=None, h=None):
    """All roots of the self-consistency (*), with stability and eye insulin."""
    gain = E.ep.INDUCE_GAIN if gain is None else gain
    ki = E.ep.S6K_KI if ki is None else ki
    h = E.ep.H_S6K if h is None else h
    mt, ins = _mtor_ins_of_f(Fb, fgrid)
    fb = ki ** h / (ki ** h + mt ** h)
    produced = E.EYE_FMITO + gain * mt * fb          # RHS of (*)
    g = produced - fgrid                              # zero at fixed points
    out = []
    for i in range(len(fgrid) - 1):
        if g[i] == 0 or g[i] * g[i + 1] < 0:
            t = -g[i] / (g[i + 1] - g[i])
            f_star = fgrid[i] + t * (fgrid[i + 1] - fgrid[i])
            slope = (produced[i + 1] - produced[i]) / (fgrid[i + 1] - fgrid[i])
            ins_star = ins[i] + t * (ins[i + 1] - ins[i])
            out.append(dict(f=f_star, eye_ins=ins_star, stable=bool(slope < 1.0)))
    return out, (fgrid, produced, ins)


if __name__ == '__main__':
    import os
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    fgrid = np.linspace(0.02, 3.0, 60)
    Fb_mgL = np.array([0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15, 20, 30])
    print('reference params: gain=%.1f  S6K_KI=%.2f  h=%.1f' %
          (E.ep.INDUCE_GAIN, E.ep.S6K_KI, E.ep.H_S6K))
    rows = []
    for mgl in Fb_mgL:
        fps, _ = fixed_points(mgl / 46.03, fgrid)
        n = len(fps)
        tag = 'BISTABLE' if n >= 3 else ('mono' if n == 1 else f'{n}')
        stab = [fp for fp in fps if fp['stable']]
        eyes = ','.join('%.2f' % fp['eye_ins'] for fp in fps)
        print('  blood %5.1f mg/L: %d fixed point(s) [%s]  eye insulin={%s}'
              % (mgl, n, tag, eyes))
        for fp in fps:
            rows.append(dict(blood_formate_mgL=mgl, **fp))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, 'eye_bifurcation.csv'), index=False)
    nb = df.groupby('blood_formate_mgL').size()
    bistable = nb[nb >= 3]
    if len(bistable):
        print('\nBISTABLE (switch, hysteresis) for blood formate in [%.1f, %.1f] mg/L'
              % (bistable.index.min(), bistable.index.max()))
    else:
        print('\nNo bistability found: the response is monostable (steep but single-valued).')
