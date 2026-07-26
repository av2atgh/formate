#!/usr/bin/env python3
"""Whole-body / systemic figures fig6-fig9 from the model classes, styled with
no in-figure titles (the caption carries the title) and a)/b)/c) panel labels."""
import os
import numpy as np
import matplotlib.pyplot as plt
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from model import BetaCellModel, FormateDiabetes, CoupledModel, DIETS
from figstyle import setup, panel

setup()
GRAY, BLUE, RED, GREEN, PURPLE = '#7f7f7f', '#1f77b4', '#d62728', '#2ca02c', '#9467bd'
bc = BetaCellModel(); bc.calibrate_all(verbose=False)


def fig6_formate_diabetes():
    df, _ = FormateDiabetes(beta=bc).perturbation()
    cols = {'high': RED, 'normal': GREEN, 'low': BLUE}
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.4), constrained_layout=True)
    for diet in ['high', 'normal', 'low']:
        d = df[df.diet == diet]
        ax[0].plot((1 - d.sec_defect) * 100, d.formate_mgL, '-o', ms=3,
                   color=cols[diet], label='%s intake' % diet)
        ax[1].plot(d.glucose_mM, d.formate_mgL, '-o', ms=3, color=cols[diet],
                   label='%s intake' % diet)
    ax[0].set(xlabel='insulin-secretion defect (%)', ylabel='blood formate (mg/L)')
    ax[0].legend(fontsize=8, loc='upper left')
    ax[1].set(xlabel='blood glucose (mM)', ylabel='blood formate (mg/L)')
    panel(ax[0], 'a', 'vs secretion defect'); panel(ax[1], 'b', 'vs blood glucose')
    fig.savefig('fig6_formate_diabetes.pdf')


def fig7_obesity_diabetes():
    df = CoupledModel(beta=bc).scenarios()
    short = ['lean\nnon-diabetic', 'obese\nnon-diabetic', 'lean\ndiabetic', 'obese\ndiabetic']
    x = np.arange(4); cols = [GRAY, BLUE, RED, PURPLE]
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.6), constrained_layout=True)
    for a, col, pct, ylab, pad in [
            (ax[0], 'formate_mgL', 'formate_vs_healthy_%', 'blood formate (mg/L)', 0.06),
            (ax[1], 'urate_mgdL', 'urate_vs_healthy_%', 'blood uric acid (mg/dL)', 0.02)]:
        vals = df[col].to_numpy()
        a.bar(x, vals, color=cols)
        a.axhline(vals[0], ls=':', color='k', lw=0.8)
        for i in range(4):
            a.annotate('%+d%%' % round(df[pct].iloc[i]), (i, vals[i]),
                       textcoords='offset points', xytext=(0, 3), ha='center', fontsize=7.5)
        a.set_xticks(x); a.set_xticklabels(short, fontsize=8); a.set_ylabel(ylab)
        a.set_ylim(0, vals.max() * 1.12 if col == 'formate_mgL' else None)
    ax[1].set_ylim(3.30, df.urate_mgdL.max() * 1.02)
    panel(ax[0], 'a', 'blood formate'); panel(ax[1], 'b', 'blood uric acid')
    fig.savefig('fig7_obesity_diabetes.pdf')


def fig8_intervention():
    df = CoupledModel(beta=bc).intervention()
    styles = {'lean, non-diabetic': (GRAY, '-o', 'lean, non-diabetic'),
              'lean diabetic, formate-replete': (BLUE, '-s', 'lean diabetic, formate-replete'),
              'lean diabetic, formate-deficient': (RED, '-o', 'lean diabetic, formate-deficient'),
              'obese diabetic, formate-deficient': (PURPLE, '--^', 'obese diabetic, formate-deficient')}
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.5), constrained_layout=True)
    for g, (c, st, lab) in styles.items():
        d = df[df.group == g]; xx = d.D_diet / DIETS['normal']
        ax[0].plot(xx, d.glucose_mM, st, ms=3, color=c, label=lab)
        ax[1].plot(xx, d.urate_mgdL, st, ms=3, color=c)
    ax[0].set(xlabel='formate/serine intake ($\\times$ normal)', ylabel='fasting glucose (mM)')
    ax[1].set(xlabel='formate/serine intake ($\\times$ normal)', ylabel='blood uric acid (mg/dL)')
    ax[0].legend(fontsize=6, loc='lower center', bbox_to_anchor=(0.5, 0.06), labelspacing=0.25)
    panel(ax[0], 'a', 'fasting glucose'); panel(ax[1], 'b', 'blood uric acid')
    fig.savefig('fig8_intervention.pdf')


def fig9_metformin():
    df = CoupledModel(beta=bc).metformin()
    pts = {'lean diabetic': (RED, '-o'), 'obese diabetic': (PURPLE, '--^')}
    fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.3), constrained_layout=True)
    for name, (c, st) in pts.items():
        d = df[df.patient == name]
        ax[0].plot(d.metformin_dose, d.glucose_mM, st, ms=3, color=c, label=name)
        ax[1].plot(d.metformin_dose, d.formate_mgL, st, ms=3, color=c)
        ax[2].plot(d.metformin_dose, d.insulin, st, ms=3, color=c)
    ax[0].set(xlabel='metformin dose (a.u.)', ylabel='fasting glucose (mM)')
    ax[1].set(xlabel='metformin dose (a.u.)', ylabel='blood formate (mg/L)')
    ax[2].set(xlabel='metformin dose (a.u.)', ylabel='plasma insulin (a.u.)')
    ax[0].legend(fontsize=7.5, loc='upper right')
    panel(ax[0], 'a', 'glucose falls'); panel(ax[1], 'b', 'formate falls'); panel(ax[2], 'c', 'insulin falls')
    fig.savefig('fig9_metformin.pdf')


if __name__ == '__main__':
    fig6_formate_diabetes(); print('fig6')
    fig7_obesity_diabetes(); print('fig7')
    fig8_intervention(); print('fig8')
    fig9_metformin(); print('fig9')
    print('DONE')
