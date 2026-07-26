#!/usr/bin/env python3
"""Dietary formate analysis and figures (fig14/fig15/fig16).

Beverage/medication formate contents are taken from Pietzke et al. 2019
(Mol Metab 33:23), Table 2 -- the potential formate generation from
caffeine/aspartame demethylation. That table reports formaldehyde (mmol per
serving) and a formate production rate (umol/h/kg); we take the formate PER
SERVING as formaldehyde(mmol) x 46.03 mg/mmol (full conversion, footnote c).
Whole-food items (fruit juice, red meat, creatine) are not in Table 2 and keep
the book estimates (The Spice of Life, p.53).

A serving's formate is treated as a daily intake and mapped to the coupled
model's dietary intake rate by D_diet = 3e-5 * mg (so 1 g formate/day <-> the
normal endogenous-scale intake 0.030 mM/h). Blood formate, plasma insulin and
glucose come from the lean-healthy fixed point of the three-compartment model
(obesity_diabetes.py); eye insulin from the regulated autonomous eye
(eye.eye_state_ff).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from model import BetaCellModel, CoupledModel, EyeModel

_BC = BetaCellModel(); _BC.calibrate_all(verbose=False)
CM = CoupledModel(beta=_BC)      # three-compartment lean-healthy fixed point
E = EyeModel(beta=_BC)           # regulated autonomous eye

MW = 46.03                      # formate molecular weight (mg/mmol)
D_PER_MG = 3e-5                 # dietary intake rate (mM/h) per mg formate/serving

def fmg(formaldehyde_mmol):     # Table 2 formate per serving (mg), full conversion
    return round(formaldehyde_mmol * MW, 1)

# ---- foods: Table 2 beverages/medications (by formaldehyde mmol) + book foods ----
FOODS = [
    # name, formate mg/serving, category
    ('low-formate diet',        5.0,            'reference'),
    ('Pepsi 12oz',              fmg(0.16),      'soft drink'),
    ('Coke 12oz',               fmg(0.18),      'soft drink'),
    ('Dr Pepper 12oz',          fmg(0.23),      'soft drink'),
    ('Espresso 1oz',            fmg(0.31),      'coffee'),
    ('Analgesic (caffeine)',    fmg(0.51),      'medication'),
    ('Caffeine-free diet Pepsi',fmg(0.58),      'diet soft drink'),
    ('Brewed coffee 8oz',       fmg(0.67),      'coffee'),
    ('Diet Pepsi 12oz',         fmg(0.72),      'diet soft drink'),
    ('Diet Mountain Dew 12oz',  fmg(0.81),      'diet soft drink'),
    ('Diet Coke 12oz',          fmg(0.93),      'diet soft drink'),
    ('Caffeine pill (NoDoz)',   fmg(1.03),      'medication'),
    ('Fruit juice 1L',          300.0,          'whole food'),
    ('Animal meat 200g',        700.0,          'whole food'),
    ('Creatine 1g',             1000.0,         'whole food'),
]

# ---- meta-analysis: items with a per-serving T2D relative risk ----
# formate contents updated to Table 2 (beverages); RR from the meta-analysis lit.
META = [
    ('Coffee (brewed)',   fmg(0.67), 0.91, 0.89, 0.94, 'per cup, Ding 2014'),
    ('Regular cola',      fmg(0.18), 1.13, 1.06, 1.21, 'SSB/serving, Imamura 2015'),
    ('Diet cola',         fmg(0.93), 1.08, 1.02, 1.15, 'ASB/serving, Imamura 2015'),
    ('Fruit juice',       300.0,     1.07, 1.01, 1.14, 'per serving, Imamura 2015'),
    ('Red meat',          700.0,     1.10, 1.06, 1.15, 'per 100 g, InterConnect 2024'),
    ('Creatine',          1000.0,    1.00, 0.90, 1.11, 'glycemic-neutral, Delpino 2022'),
]


def build():
    base = CM.NBC + CM.NAD               # blood pools: Gb, Ip, X, Fb, Ub
    Ib, KI, y0 = CM._reference(0.030)
    rows, y = [], y0
    for name, mg, cat in FOODS:
        D = D_PER_MG * mg
        y = CM._settle(1.0, D, Ib, KI, 1.0, 1.0, CM.wp.Si, y)
        Gb, Ip, Fb = y[base], y[base + 1], y[base + 3]
        s, _ = E.eye_state_ff(5.0, Fb, y0=None)
        rows.append(dict(food=name, mg_formate=mg, category=cat,
                         blood_formate_mgL=Fb * MW, plasma_insulin=Ip,
                         glucose_mM=Gb, eye_ins_prod=s['eye_ins_prod']))
    dm = pd.DataFrame(rows); dm.to_csv('data/diet_model.csv', index=False)

    md = pd.DataFrame([dict(food=n, formate_mg=mg, RR=rr, lo=lo, hi=hi, source=src,
                            logRR=float(np.log(rr))) for n, mg, rr, lo, hi, src in META])
    md.to_csv('data/diet_meta.csv', index=False)
    r = float(np.corrcoef(md.formate_mg, md.logRR)[0, 1])
    return dm, md, r


def figures(dm, md, r):
    from figstyle import setup
    setup()
    S = np.linspace(1.4, 5, 3000)
    Sopt = S[np.array([E.retinal_outcome(s)[2] for s in S]).argmax()]
    CCOL = {'reference': '#bbbbbb', 'soft drink': '#1f77b4', 'diet soft drink': '#17becf',
            'coffee': '#8c564b', 'medication': '#9467bd', 'whole food': '#2ca02c'}

    # fig14: blood formate by food (log y; with insulin-saturation reference)
    # blood-formate level that saturates beta-cell insulin production (Fig. 2)
    grid = np.linspace(0.0, 2.0, 30)
    fs = _BC.formate_scan(11.0, grid, f_mito=0.2 * _BC.p.FMITO)
    v = fs.v_ins_syn.to_numpy()
    sat_mgL = grid[np.argmax(v >= 0.95 * v.max())] * MW      # ~14 mg/L
    d = dm[dm.food != 'low-formate diet'].sort_values('mg_formate').reset_index(drop=True)
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(8.6, 4.4), constrained_layout=True)
    ax.bar(x, d.blood_formate_mgL, color=[CCOL[c] for c in d.category])
    for i in range(len(d)):
        ax.text(i, d.blood_formate_mgL.iloc[i] * 1.06, f'{d.mg_formate.iloc[i]:g}',
                ha='center', va='bottom', fontsize=6.5, color='#333')
    ax.axhline(sat_mgL, ls='--', color='k', lw=1.1)
    ax.text(len(d) - 0.4, sat_mgL * 1.08,
            '$\\beta$-cell insulin-production saturation (Fig. 2)',
            ha='right', va='bottom', fontsize=7.5, color='k')
    ax.set_yscale('log'); ax.set_ylim(0.3, sat_mgL * 2.4)
    ax.set_xticks(x); ax.set_xticklabels([f.replace(' ', '\n', 1) for f in d.food],
                                         rotation=40, ha='right', fontsize=7)
    ax.set_ylabel('blood formate (mg/L)')
    from matplotlib.patches import Patch
    cats = ['soft drink', 'diet soft drink', 'coffee', 'medication', 'whole food']
    ax.legend(handles=[Patch(color=CCOL[c], label=c) for c in cats],
              frameon=False, fontsize=7.5, loc='center left')
    fig.savefig('fig14_diet_formate.pdf')

    # fig15: forest plot, updated formate contents + recomputed r
    m = md.sort_values('formate_mg').reset_index(drop=True)
    y = np.arange(len(m))
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    for i, row in m.iterrows():
        col = '#2ca02c' if row.RR < 1 else ('#d62728' if row.lo > 1 else '#7f7f7f')
        ax.plot([row.lo, row.hi], [i, i], '-', color=col, lw=2)
        ax.plot(row.RR, i, 'o', color=col, ms=7)
    ax.axvline(1.0, ls=':', color='k', lw=1)
    ax.set_yticks(y); ax.set_yticklabels([f'{f} ({int(mg)} mg)' for f, mg in
                                          zip(m.food, m.formate_mg)], fontsize=8.5)
    ax.set_xlabel('relative risk of type 2 diabetes (per serving)')
    ax.text(0.03, 0.03, f'formate vs risk: $r$ = {r:+.2f}', transform=ax.transAxes,
            fontsize=8.5, va='bottom')
    fig.savefig('fig15_diet_forest.pdf')

    # fig16: eye insulin by food (regulated), sorted
    fig, ax = plt.subplots(figsize=(8.6, 4.4), constrained_layout=True)
    ax.plot(x, d.eye_ins_prod, '-o', ms=5, color='#d62728',
            label='eye insulin (predicted; unmeasured)')
    ax.axhline(Sopt, ls=':', color='#2ca02c', lw=1.3)
    ax.text(0.1, Sopt - 0.004, 'retinal optimum', fontsize=8, color='#2ca02c', va='top')
    ax.set_xticks(x); ax.set_xticklabels([f.replace(' ', '\n', 1) for f in d.food],
                                         rotation=40, ha='right', fontsize=7)
    ax.set_ylabel('eye insulin production (a.u.)')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax.set_ylim(1.90, Sopt + 0.01)
    fig.savefig('fig16_eye_insulin.pdf')


if __name__ == '__main__':
    dm, md, r = build()
    figures(dm, md, r)
    print('corr(formate content, log-RR T2D) = %.3f' % r)
    print('blood formate range: %.3f - %.3f mg/L' %
          (dm.blood_formate_mgL.min(), dm.blood_formate_mgL.max()))
    print('plasma insulin range: %.4f - %.4f (flat)' %
          (dm.plasma_insulin.min(), dm.plasma_insulin.max()))
    print('eye insulin range: %.3f - %.3f' %
          (dm.eye_ins_prod.min(), dm.eye_ins_prod.max()))
    print(dm[['food', 'mg_formate', 'blood_formate_mgL', 'eye_ins_prod']].to_string(index=False))
