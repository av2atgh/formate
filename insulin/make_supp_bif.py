#!/usr/bin/env python3
"""Supplementary Figure S1: the regulated eye feed-forward is monostable.
Self-consistency (single crossing), basin test (all initial states converge),
monostability across gains, and the resulting single-valued homeostat.
Uses the formalised model classes (model.py) + eye_bifurcation.py."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from model import BetaCellModel, EyeModel
import eye_bifurcation as B

bc = BetaCellModel(); bc.calibrate_all(verbose=False)
E = EyeModel(beta=bc)
from figstyle import setup, panel
setup()
S = np.linspace(1.4, 5, 3000)
Sopt = S[np.array([E.retinal_outcome(s)[2] for s in S]).argmax()]
GAIN, KI, H = E.ep.INDUCE_GAIN, E.ep.S6K_KI, E.ep.H_S6K

fig, ax = plt.subplots(2, 2, figsize=(9.6, 7.4), constrained_layout=True)

# (A) self-consistency: produced f_mito vs imposed f_mito, several blood formate
fg = np.linspace(0.02, 3.0, 60)
fgc = np.linspace(0.005, 0.40, 60)
for mgl, c in [(1.0, '#1f77b4'), (3.0, '#2ca02c'), (10.0, '#d62728')]:
    mt, _ = B._mtor_ins_of_f(mgl / 46.03, fgc)
    fb = KI ** H / (KI ** H + mt ** H)
    produced = E.EYE_FMITO + GAIN * mt * fb
    ax[0, 0].plot(fgc, produced, '-', color=c, lw=1.8, label='blood formate %g mg/L' % mgl)
    g = produced - fgc
    idx = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
    for i in idx:
        t = -g[i] / (g[i + 1] - g[i]); fs = fgc[i] + t * (fgc[i + 1] - fgc[i])
        ax[0, 0].plot(fs, fs, 'o', color=c, ms=7, mfc='white', mew=1.6, zorder=5)
ax[0, 0].plot(fgc, fgc, 'k--', lw=1, label='self-consistency $f=f$')
ax[0, 0].set(xlabel='imposed local formate $f$ (a.u.)',
             ylabel='produced $f_0+g\\,m(f)\\,\\phi(m)$', xlim=(0, 0.4), ylim=(0, 0.4))
ax[0, 0].legend(frameon=False, fontsize=7, loc='upper left')

# (B) basin test: trajectories from many initial scalings converge to one attractor
f = 1.76 / 46.03
t = np.linspace(0, 250, 400)
for sc, a in zip([0.2, 0.5, 1.0, 2.0, 5.0, 10.0], np.linspace(0.35, 1, 6)):
    y0 = np.maximum(bc.y0 * sc, 1e-9)
    sol = solve_ivp(lambda tt, yy: E._rhs_ff(tt, yy, 5.0, f, GAIN, KI, 1.0, 1.0),
                    (0, 250), y0, method='BDF', rtol=1e-8, atol=1e-10, t_eval=t)
    ins = np.array([bc.fluxes(sol.y[:, k], 5.0, f)['v_ins_syn'] for k in range(sol.y.shape[1])])
    ax[0, 1].plot(t, ins, '-', color=plt.cm.viridis(a), lw=1.2)
ax[0, 1].axhline(1.944, ls=':', color='k', lw=1)
ax[0, 1].text(150, 1.944 + 0.02, 'single stable attractor', fontsize=7.5)
ax[0, 1].set(xlabel='time (h)', ylabel='eye insulin (a.u.)')

# (C) number of fixed points vs feed-forward gain
gains = [2, 3, 4, 6, 8, 10, 14, 20]
nfp = [len(B.fixed_points(3.0 / 46.03, fg, gain=g)[0]) for g in gains]
ax[1, 0].plot(gains, nfp, '-o', ms=6, color='#d62728')
ax[1, 0].axvline(GAIN, ls='--', color='#2ca02c', lw=1)
ax[1, 0].text(GAIN + 0.3, 2.4, 'reference', fontsize=7.5, color='#2ca02c')
ax[1, 0].axhline(3, ls=':', color='gray', lw=0.8)
ax[1, 0].text(15, 3.05, 'bistability threshold', fontsize=7, color='gray')
ax[1, 0].set(xlabel='ATF4 feed-forward gain $g$', ylabel='number of steady states',
             ylim=(0.5, 3.3), yticks=[1, 2, 3])

# (D) eye insulin vs blood formate: single smooth branch (homeostat)
Fb = np.concatenate([np.linspace(0.3, 3, 10), np.linspace(4, 30, 12)])
ins = [E.eye_state_ff(5.0, mgl / 46.03, y0=None)[0]['eye_ins_prod'] for mgl in Fb]
ax[1, 1].plot(Fb, ins, '-^', ms=4, color='#2ca02c')
ax[1, 1].axhline(Sopt, ls=':', color='#2ca02c', lw=1.1)
ax[1, 1].text(0.5, Sopt + 0.02, 'retinal optimum', fontsize=7.5, color='#2ca02c')
ax[1, 1].set_xscale('log')
ax[1, 1].set(xlabel='blood formate (mg/L)', ylabel='eye insulin (a.u.)')

panel(ax[0,0],'a','self-consistency'); panel(ax[0,1],'b','basin test')
panel(ax[1,0],'c','fixed points vs gain'); panel(ax[1,1],'d','homeostat')
fig.savefig('figS_eye_bifurcation.pdf')
print('wrote figS_eye_bifurcation.pdf ; gains', gains, 'nfp', nfp)
