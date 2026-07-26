#!/usr/bin/env python3
"""Gut-microbiome urate cycle (fig22), 2x2: (a) schematic of the delayed
positive-feedback model; (b) steady-state formate recovered (lean/obese, cycle
off/on); (c) blood formate up / urate down as the gut recycling strengthens
(bounded feedback); (d) a urate/purine bolus returns as a delayed formate rise.
Uses model.py GutUrateCycle."""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.integrate import solve_ivp
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from model import BetaCellModel, GutUrateCycle
from figstyle import setup, panel

setup()
GRAY, BLUE, RED, GREEN, PURPLE = '#7f7f7f', '#1f77b4', '#d62728', '#2ca02c', '#9467bd'
MWF, MWU = 46.03, 16.8
bc = BetaCellModel(); bc.calibrate_all(verbose=False)
nc = 30

fig, ax = plt.subplots(2, 2, figsize=(11.6, 8.0), constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.05, h_pad=0.05, wspace=0.12, hspace=0.14)


# (a) schematic of the feedback model
def schematic(a):
    a.axis('off'); a.set_xlim(0, 1); a.set_ylim(0, 1)
    nodes = {'formate': (0.15, 0.56, 'blood\nformate', GREEN),
             'purine': (0.5, 0.89, 'host de novo\npurine synthesis', BLUE),
             'urate': (0.85, 0.56, 'blood\nurate', PURPLE),
             'gut': (0.5, 0.23, 'gut\nmicrobiome', RED)}
    pos = {}
    for k, (xx, yy, lab, col) in nodes.items():
        a.text(xx, yy, lab, ha='center', va='center', fontsize=7.5, color='k',
               bbox=dict(boxstyle='round,pad=0.34', fc='white', ec=col, lw=1.5))
        pos[k] = (xx, yy)

    def arr(u, v, color='#333', rad=-0.18, lw=1.5, ls='-', sB=34):
        (x0, y0), (x1, y1) = pos[u], pos[v]
        a.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                    connectionstyle=f'arc3,rad={rad}', arrowstyle='-|>',
                    mutation_scale=13, color=color, lw=lw, ls=ls,
                    shrinkA=34, shrinkB=sB))
    arr('formate', 'purine', rad=-0.18, sB=40)     # wide purine box
    arr('purine', 'urate', rad=-0.18, sB=30)
    arr('urate', 'gut', rad=-0.18, sB=30)
    arr('gut', 'formate', color=GREEN, rad=-0.18, lw=2.4, ls='--', sB=30)  # delayed return
    a.text(0.24, 0.80, '2 formate\n+ glycine', ha='center', fontsize=6.5, color='#333')
    a.text(0.72, 0.80, 'catabolism', ha='left', fontsize=6.5, color='#333')
    a.text(0.72, 0.30, 'urate\nuptake', ha='left', fontsize=6.5, color='#333')
    a.text(0.24, 0.29, '+1 formate\ndelay $\\tau_{gut}$', ha='center', va='center',
           fontsize=7, color=GREEN, fontweight='bold')
    a.text(0.5, 0.085, '+ acetate\n$\\to$ lipogenesis (obesity)', ha='center',
           va='bottom', fontsize=6.5, color=RED)
    a.text(0.95, 0.50, 'renal\nloss', ha='left', va='center', fontsize=6.5, color='#888')
    a.text(0.95, 0.68, 'obesity:\nadipose XOR$\\uparrow$', ha='left', va='center',
           fontsize=6.5, color=PURPLE)


schematic(ax[0, 0])
panel(ax[0, 0], 'a', 'urate-cycle feedback model')

# (b) steady-state formate: lean/obese x cycle off/on
lean_f, obese_f = [], []
for kg in [0.0, 0.15]:
    g = GutUrateCycle(beta=bc, k_gut=kg)
    Ib, KI, y0 = g._reference(g.DIETS['normal'])
    yl = g._settle(1.0, g.DIETS['normal'], Ib, KI, 1.0, 1.0, g.wp.Si, y0)
    yo = g._settle(1.0, g.DIETS['normal'], Ib, KI, 3.0, 2.0, 0.5 * g.wp.Si, y0)
    lean_f.append(yl[nc - 2] * MWF); obese_f.append(yo[nc - 2] * MWF)
x = np.arange(2); w = 0.36
ax[0, 1].bar(x - w / 2, [lean_f[0], obese_f[0]], w, color=GRAY, label='cycle off')
ax[0, 1].bar(x + w / 2, [lean_f[1], obese_f[1]], w, color=GREEN, label='cycle on')
ax[0, 1].set_xticks(x); ax[0, 1].set_xticklabels(['lean', 'obese'])
ax[0, 1].set_ylabel('blood formate (mg/L)')
ax[0, 1].legend(fontsize=7.5, loc='upper right')
panel(ax[0, 1], 'b', 'urate carbon recovered as formate')

# (c) blood formate up / urate down as gut recycling strengthens (bounded feedback)
kgs = np.linspace(0.0, 0.30, 10)
Fb, Ub = [], []
for kg in kgs:
    g = GutUrateCycle(beta=bc, k_gut=kg)
    _, _, y = g._reference(g.DIETS['normal'])
    Fb.append(y[nc - 2] * MWF); Ub.append(y[nc - 1] * MWU)
ax[1, 0].plot(kgs, Fb, '-o', ms=3, color=GREEN, label='blood formate')
ax[1, 0].set_xlabel('gut urate$\\to$formate recycling  $k_{gut}$ (1/h)')
ax[1, 0].set_ylabel('blood formate (mg/L)', color=GREEN)
ax[1, 0].tick_params(axis='y', labelcolor=GREEN)
axc = ax[1, 0].twinx(); axc.spines['top'].set_visible(False)
axc.plot(kgs, Ub, '-s', ms=3, color=PURPLE)
axc.set_ylabel('blood uric acid (mg/dL)', color=PURPLE)
axc.tick_params(axis='y', labelcolor=PURPLE)
panel(ax[1, 0], 'c', 'bounded positive feedback')

# (d) delay: a urate/purine bolus returns as a delayed formate rise
g = GutUrateCycle(beta=bc, k_gut=0.15, tau_gut=12.0)
Ib, KI, y0 = g._reference(g.DIETS['normal'])
args = (1.0, g.DIETS['normal'], Ib, KI, 1.0, 1.0, g.wp.Si, None, 1.0, None)
yb = y0.copy(); yb[nc - 1] += 0.15
t = np.linspace(0, 60, 200)
sol = solve_ivp(g._rhs, (0, 60), yb, method='BDF', rtol=1e-6, atol=1e-9, args=args, t_eval=t)
Ubt = sol.y[nc - 1] * MWU; Fbt = sol.y[nc - 2] * MWF
ax[1, 1].plot(t, Ubt, '-', color=PURPLE, lw=2, label='blood urate (load)')
ax[1, 1].set_xlabel('time after purine/urate load (h)')
ax[1, 1].set_ylabel('blood uric acid (mg/dL)', color=PURPLE)
ax[1, 1].tick_params(axis='y', labelcolor=PURPLE)
axd = ax[1, 1].twinx(); axd.spines['top'].set_visible(False)
axd.plot(t, Fbt, '-', color=GREEN, lw=2)
axd.set_ylabel('blood formate (mg/L)', color=GREEN)
axd.tick_params(axis='y', labelcolor=GREEN)
kpk = t[Fbt.argmax()]
axd.axvline(kpk, ls=':', color=GREEN, lw=1)
axd.annotate('formate returns\n${\\sim}\\tau_{gut}$ later', (kpk, Fbt.max()),
             xytext=(kpk + 6, Fbt.max()), fontsize=7, color=GREEN, va='center')
panel(ax[1, 1], 'd', 'delayed formate return')

fig.savefig('fig22_urate_cycle.pdf')
print('lean formate off/on %.2f/%.2f ; obese %.2f/%.2f ; delayed peak t=%.0f h'
      % (lean_f[0], lean_f[1], obese_f[0], obese_f[1], kpk))
