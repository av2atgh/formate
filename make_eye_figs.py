#!/usr/bin/env python3
"""Eye-compartment figures (fig17 expression, fig18 homeostat, fig10 clinical,
fig13 retinal-health window, fig19 methanol). No in-figure titles; a)/b) panel
labels where multi-panel. fig16 (diet eye insulin) is in make_diet.py."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from model import BetaCellModel, EyeModel
from figstyle import setup, panel

setup()
bc = BetaCellModel(); bc.calibrate_all(verbose=False)
E = EyeModel(beta=bc)
MW = 46.03
S = np.linspace(1.4, 5, 3000)
Sopt = S[np.array([E.retinal_outcome(s)[2] for s in S]).argmax()]
print('retinal optimum Sopt=%.3f' % Sopt)


def fig17_expression():
    d = pd.read_csv('data/eye_onecarbon_expr.csv')
    x = np.arange(len(d)); w = 0.27
    fig, ax = plt.subplots(figsize=(7.4, 3.9), constrained_layout=True)
    ax.bar(x - w, d.SHMT2, w, color='#2ca02c', label='SHMT2 (serine$\\to$1C entry)')
    ax.bar(x, d.MTHFD2, w, color='#1f77b4', label='MTHFD2')
    ax.bar(x + w, d.MTHFD1L, w, color='#d62728', label='MTHFD1L (formate release)')
    ax.set_yscale('log'); ax.set_ylim(0.8, 3000)
    ax.set_xticks(x); ax.set_xticklabels(d.cell, rotation=30, ha='right')
    ax.set_ylabel('single-cell expression (nTPM, HPA)')
    ax.legend(fontsize=8, loc='upper right')
    ri = list(d.cell).index('RPE')
    ax.annotate('insulin-producing RPE:\nentry enzymes,\nlittle MTHFD1L',
                xy=(ri, 13), xytext=(ri - 0.45, 260), ha='left', va='center',
                fontsize=7, color='#333',
                arrowprops=dict(arrowstyle='->', color='#888', lw=0.7))
    ci = list(d.cell).index('Cone')
    ax.annotate('MTHFD1L peaks in cones\n(%.0f, the eye\'s highest)' % d.MTHFD1L[ci],
                xy=(ci + w, d.MTHFD1L[ci]), xytext=(ci + 0.5, 320), ha='left',
                va='center', fontsize=7, color='#d62728',
                arrowprops=dict(arrowstyle='->', color='#d62728', lw=0.7))
    fig.savefig('fig17_eye_onecarbon.pdf')


def fig18_local_production():
    Fb = np.linspace(0.002, 0.45, 24)
    rows, yb, yu, yr = [], None, None, None
    for f in Fb:
        sb, yb = E.eye_state(5.0, f, y0=yb)
        su, yu = E.eye_state_ff(5.0, f, s6k_ki=None, y0=yu)
        sr, yr = E.eye_state_ff(5.0, f, y0=yr)
        rows.append({'mgL': f * MW, 'blood': sb['eye_ins_prod'],
                     'unreg': su['eye_ins_prod'], 'reg': sr['eye_ins_prod']})
    d = pd.DataFrame(rows); d.to_csv('data/eye_local_production.csv', index=False)
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    ax.plot(d.mgL, d.blood, '-o', ms=3, color='#1f77b4', label='blood-dependent (no local production)')
    ax.plot(d.mgL, d.unreg, '-s', ms=3, color='#d62728', label='local, UNregulated (runs to excess)')
    ax.plot(d.mgL, d.reg, '-^', ms=3, color='#2ca02c', label='local, S6K1-IRS1 regulated (homeostat)')
    ax.axhline(Sopt, ls=':', color='#2ca02c', lw=1.1)
    ax.text(0.9, Sopt + 0.02, 'retinal optimum', fontsize=7.5, color='#2ca02c')
    ax.set(xlabel='blood formate (mg/L)', ylabel='eye insulin production (a.u.)')
    ax.legend(fontsize=7.5, loc='center right')
    fig.savefig('fig18_eye_local.pdf')


def fig10_conditions():
    obd = pd.read_csv('data/obesity_diabetes.csv')
    rows = []
    for _, c in obd.iterrows():
        s, _ = E.eye_state_ff(c['glucose_mM'], c['formate_mM'], y0=None)
        rows.append({'cond': c['scenario'], 'Ip': c['insulin'], 'Fb': c['formate_mgL'],
                     'eye': s['eye_ins_prod']})
    ac = pd.DataFrame(rows); ac.to_csv('data/eye_conditions.csv', index=False)
    short = ['lean\nhealthy', 'obese\nnon-diab', 'lean\ndiabetic', 'obese\ndiabetic']
    x = np.arange(4); cols = ['#7f7f7f', '#1f77b4', '#d62728', '#9467bd']
    fig, ax = plt.subplots(figsize=(6.0, 4.2), constrained_layout=True)
    ax.bar(x, ac.eye, color=cols)
    ax.axhline(Sopt, ls=':', color='#2ca02c', lw=1.3)
    ax.text(3.4, Sopt + 0.01, 'retinal optimum', fontsize=7.5, color='#2ca02c', ha='right')
    for i in range(4):
        ax.annotate(f'$I_p$={ac.Ip.iloc[i]:.2f}\nFb={ac.Fb.iloc[i]:.1f}', (i, ac.eye.iloc[i]),
                    textcoords='offset points', xytext=(0, 3), ha='center', fontsize=6.6, color='#333')
    ax.set_xticks(x); ax.set_xticklabels(short, fontsize=8); ax.set_ylim(1.6, 2.72)
    ax.set_ylabel('eye insulin production (a.u.)')
    fig.savefig('fig10_eye.pdf')


def fig13_window():
    d = E.insulin_excess_window()
    Sg = np.linspace(1.5, 2.7, 400)
    m = np.array([E.retinal_outcome(s)[0] for s in Sg])
    p = np.array([E.retinal_outcome(s)[1] for s in Sg])
    h = m - p
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)
    ax[0].plot(d.blood_formate_mgL, d.maintenance, '-', color='#2ca02c', lw=2, label='maintenance (trophic)')
    ax[0].plot(d.blood_formate_mgL, d.pathology, '-', color='#d62728', lw=2, label='pathology (VEGF/gliosis)')
    ax[0].plot(d.blood_formate_mgL, d.retinal_health, '-', color='#1f77b4', lw=2.6, label='net retinal health')
    ax[0].set_xscale('log')
    ax[0].set(xlabel='blood formate (mg/L)', ylabel='index')
    ax[0].legend(fontsize=7, loc='center right')
    ax[1].plot(Sg, h, '-', color='#1f77b4', lw=2.4)
    ax[1].axvline(Sopt, ls=':', color='#2ca02c', lw=1)
    ax[1].fill_between(Sg, h, where=(Sg < Sopt), color='#1f77b4', alpha=0.06)
    ax[1].fill_between(Sg, h, where=(Sg > Sopt), color='#d62728', alpha=0.06)
    ax[1].text(1.62, 0.255, 'deficiency', fontsize=8, color='#1f77b4')
    ax[1].text(2.58, 0.255, 'excess', fontsize=8, color='#d62728', ha='right')
    ax[1].set(xlabel='local eye insulin production', ylabel='net retinal health')
    panel(ax[0], 'a', 'window arms vs blood formate'); panel(ax[1], 'b', 'inverted-U in local insulin')
    fig.savefig('fig13_insulin_excess.pdf')


def fig19_methanol():
    Fb2 = np.geomspace(0.02, 12.0, 80)      # dense log grid -> smooth curves
    rows, yb, yr = [], None, None
    for f in Fb2:
        sb, yb = E.eye_state(5.0, f, y0=yb, cox_ki=None)
        sr, yr = E.eye_state_ff(5.0, f, y0=yr)   # warm-start (monostable)
        rows.append({'mgL': f * MW, 'ins_b': sb['eye_ins_prod'], 'ins_r': sr['eye_ins_prod'],
                     'h_b': E.retinal_outcome(sb['eye_ins_prod'])[2],
                     'h_r': E.retinal_outcome(sr['eye_ins_prod'])[2]})
    df = pd.DataFrame(rows); df.to_csv('data/eye_methanol_local_nocox.csv', index=False)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.7), constrained_layout=True)

    def marks(a):
        for v, c in [(5, '#2ca02c'), (50, '#ff7f0e'), (500, '#d62728')]:
            a.axvline(v, ls=':', color=c, lw=1)

    ax[0].plot(df.mgL, df.ins_b, '-', color='#1f77b4', lw=2, label='blood-dependent eye')
    ax[0].plot(df.mgL, df.ins_r, '-', color='#2ca02c', lw=2, label='regulated autonomous eye')
    ax[0].axhline(Sopt, ls=':', color='#2ca02c', lw=1.1)
    ax[0].text(1.2, Sopt + 0.02, 'optimum', fontsize=7, color='#2ca02c')
    ax[0].set_xscale('log')
    ax[0].set(xlabel='blood formate (mg/L)', ylabel='eye insulin (a.u.)')
    marks(ax[0]); ax[0].legend(fontsize=7, loc='center right')
    ax[1].plot(df.mgL, df.h_b, '-', color='#1f77b4', lw=2, label='blood-dependent eye')
    ax[1].plot(df.mgL, df.h_r, '-', color='#2ca02c', lw=2, label='regulated autonomous eye')
    ax[1].set_xscale('log')
    ax[1].set(xlabel='blood formate (mg/L)', ylabel='net retinal health')
    marks(ax[1]); ax[1].legend(fontsize=7, loc='lower left')
    panel(ax[0], 'a', 'eye insulin'); panel(ax[1], 'b', 'net retinal health')
    fig.savefig('fig19_methanol_local.pdf')


if __name__ == '__main__':
    fig17_expression(); print('fig17')
    fig18_local_production(); print('fig18')
    fig10_conditions(); print('fig10')
    fig13_window(); print('fig13')
    fig19_methanol(); print('fig19')
    print('DONE')
