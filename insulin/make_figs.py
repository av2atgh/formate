#!/usr/bin/env python3
"""Figures for the beta-cell insulin-synthesis model (beta_cell.py).

Run after calibration so the module parameters reflect the fitted values.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import BetaCellModel, WholeBody
from figstyle import setup, panel

M = BetaCellModel()
setup()
BLUE, RED, GREEN, PURPLE = '#1f77b4', '#d62728', '#2ca02c', '#9467bd'


def fig_calibration(path='fig1_calibration.pdf'):
    """GSIS dose-response after calibration to islet-perifusion benchmarks."""
    df = M.gsis_dose_response(G_grid=np.linspace(2.0, 20, 30))
    ec50, si = M._gsis_metrics(df)
    fig, ax = plt.subplots(figsize=(4.2, 3.2), constrained_layout=True)
    ax.plot(df.glucose, df.v_sec, '-o', ms=3, color=BLUE)
    ax.axvline(M.GSIS_EC50, ls=':', color='gray', lw=0.8)
    ax.text(M.GSIS_EC50 + 0.3, ax.get_ylim()[1] * 0.1,
            'target EC50 8 mM', fontsize=7, color='gray')
    ax.set(xlabel='glucose (mM)', ylabel='insulin secretion (a.u./h)')
    ax.text(0.97, 0.05, 'EC$_{50}$ = %.1f mM\nstim. index = %.1f' % (ec50, si),
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
    fig.savefig(path)
    return path


def fig_formate_scan(path='fig2_formate_scan.pdf'):
    """Dietary-formate rescue on a low-endogenous background (glucose 11 mM)."""
    df = M.formate_scan(11.0, np.linspace(0.0, 0.6, 25), f_mito=0.2 * M.p.FMITO)
    fig, ax = plt.subplots(1, 3, figsize=(9.5, 3.0), constrained_layout=True)
    ax[0].plot(df.form_x, df.Ppool, '-o', ms=3, color=BLUE)
    ax[0].set(xlabel='blood formate (mM)', ylabel='purine pool ATP+GXP (mM)')
    ax[1].plot(df.form_x, df.mtor, '-o', ms=3, color=GREEN, label='mTOR')
    ax[1].plot(df.form_x, df.v_ins_syn / df.v_ins_syn.max(), '-s', ms=3,
               color=RED, label='insulin synth (norm.)')
    ax[1].set(xlabel='blood formate (mM)', ylabel='activity (a.u.)')
    ax[1].legend(fontsize=8)
    ax[2].plot(df.form_x, df.v_urate_sec, '-o', ms=3, color=PURPLE)
    ax[2].set(xlabel='blood formate (mM)', ylabel='urate secretion (mM/h)')
    panel(ax[0], 'a', 'purine pool'); panel(ax[1], 'b', 'mTOR & insulin synthesis'); panel(ax[2], 'c', 'urate secretion')
    fig.savefig(path)
    return path


def fig_deficiency(path='fig3_deficiency.pdf'):
    """Mitochondrial-formate deficiency -> impaired insulin (formate theory)."""
    df = M.formate_deficiency_scan(glucose=8.0)
    x = df.f_mito / M.p.FMITO * 100        # % of normal mitochondrial formate
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    ax[0].plot(x, df.form_in / df.form_in.max(), '-^', ms=3, color=BLUE,
               label='[formate]')
    ax[0].plot(x, df.mtor / df.mtor.max(), '-s', ms=4, color=GREEN, alpha=0.5,
               label='mTOR')
    ax[0].plot(x, df.v_ins_syn / df.v_ins_syn.max(), '--', lw=1.8, color=RED,
               label='insulin synthesis')
    ax[0].set(xlabel='mitochondrial formate (% of normal)', ylabel='normalized')
    ax[0].legend(fontsize=8)
    ax[0].invert_xaxis()
    ax[1].plot(x, df.v_urate_sec, '-o', ms=3, color=PURPLE)
    ax[1].set(xlabel='mitochondrial formate (% of normal)',
              ylabel='urate secretion (mM/h)')
    ax[1].invert_xaxis()
    panel(ax[0], 'a', 'formate, mTOR, insulin'); panel(ax[1], 'b', 'urate seesaw')
    fig.savefig(path)
    return path


def fig_gtt(path='fig4_gtt.pdf'):
    """Glucose-tolerance transient: formate-deficient vs -replete cell."""
    def gtt(t):
        return 5.0 + 6.0 * np.exp(-t / 1.0)
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    for fm, c, lab in [(0.2 * M.p.FMITO, RED, 'formate-deficient'),
                       (M.p.FMITO, BLUE, 'formate-replete')]:
        _, y0 = M.steady_state(5.0, 0.1, f_mito=fm)
        tr = M.simulate(5.0, 0.1, glucose_fn=gtt, t_end=6.0, n_out=240, y0=y0,
                        f_mito=fm)
        ax[0].plot(tr.t * 60, tr.glucose, 'k--', lw=0.8)
        ax[0].plot(tr.t * 60, tr.v_sec, color=c, label=lab)
        ax[1].plot(tr.t * 60, tr.mtor, color=c, label=lab)
    ax[0].set(xlabel='time (min)', ylabel='insulin secretion (a.u./h)')
    ax[0].legend(fontsize=8)
    ax[0].text(0.98, 0.95, 'dashed = glucose', transform=ax[0].transAxes,
               ha='right', va='top', fontsize=7, color='gray')
    ax[1].set(xlabel='time (min)', ylabel='mTOR activity')
    panel(ax[0], 'a', 'insulin secretion'); panel(ax[1], 'b', 'mTOR activity')
    fig.savefig(path)
    return path


def fig_wholebody(path='fig5_wholebody.pdf'):
    """Whole-body GTT: central vs peripheral formate action are glucose-matched
    but separated by the plasma-insulin readout."""
    W = WholeBody(beta=M)
    df, trajs = W.scenarios()
    names = ['baseline (formate-deficient)', 'central (beta-cell formate)',
             'peripheral (Carpene, Sg up)']
    cols = {'baseline (formate-deficient)': 'gray',
            'central (beta-cell formate)': RED,
            'peripheral (Carpene, Sg up)': BLUE}
    short = {'baseline (formate-deficient)': 'baseline (deficient)',
             'central (beta-cell formate)': 'central (β-cell)',
             'peripheral (Carpene, Sg up)': 'peripheral (Sg)'}
    fig, ax = plt.subplots(1, 3, figsize=(10.0, 3.1), constrained_layout=True)
    for nm in names:
        tr = trajs[nm]
        ax[0].plot(tr.t * 60, tr.Gb, color=cols[nm], label=short[nm])
        ax[1].plot(tr.t * 60, tr.Ip, color=cols[nm], label=short[nm])
    ax[0].axhline(W.wp.Gb0, ls=':', color='k', lw=0.6)
    ax[0].set(xlabel='time (min)', ylabel='blood glucose (mM)')
    ax[0].legend(fontsize=7.5)
    ax[1].set(xlabel='time (min)', ylabel='plasma insulin (a.u.)')
    # discriminator bars
    d = df[df.scenario.isin(names[1:])]
    xpos = np.arange(len(d))
    ax[2].bar(xpos - 0.19, d['dGlucoseAUC_%'], 0.36, color='#7f7f7f',
              label='Δ glucose AUC')
    ax[2].bar(xpos + 0.19, d['dInsulinAUC_%'], 0.36, color=GREEN,
              label='Δ insulin AUC')
    ax[2].axhline(0, color='k', lw=0.6)
    ax[2].set_xticks(xpos)
    ax[2].set_xticklabels(['central', 'peripheral'], fontsize=8)
    ax[2].set(ylabel='% change vs baseline')
    ax[2].legend(fontsize=7.5)
    panel(ax[0], 'a', 'blood glucose'); panel(ax[1], 'b', 'plasma insulin'); panel(ax[2], 'c', 'AUC change')
    fig.savefig(path)
    return path


if __name__ == '__main__':
    M.calibrate_all()
    for f in (fig_calibration, fig_formate_scan, fig_deficiency, fig_gtt,
              fig_wholebody):
        print('wrote', f())
