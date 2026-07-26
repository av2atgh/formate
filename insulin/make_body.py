#!/usr/bin/env python3
"""Regenerate the whole-body / blood-formate / three-compartment CSV outputs
from the formalised model classes (model.py):
  * whole_body_scenarios.csv  -- central vs peripheral formate discriminator
  * formate_diabetes.csv       -- secretion defect vs blood formate
  * obesity_diabetes.csv       -- obesity x diabetes blood formate & urate
"""
import os
import numpy as np
from model import BetaCellModel, WholeBody, FormateDiabetes, CoupledModel

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs('data', exist_ok=True)
bc = BetaCellModel(); bc.calibrate_all(verbose=False)

wb = WholeBody(beta=bc)
df, _ = wb.scenarios()
df.to_csv('data/whole_body_scenarios.csv', index=False)
print('whole_body: central dInsulinAUC=%+.1f%%  peripheral dInsulinAUC=%+.1f%%'
      % (df.iloc[1]['dInsulinAUC_%'], df.iloc[2]['dInsulinAUC_%']))

fd = FormateDiabetes(beta=bc)
dff, _ = fd.perturbation()
dff.to_csv('data/formate_diabetes.csv', index=False)
r = np.corrcoef(dff.glucose_mM, dff.formate_mM)[0, 1]
print('formate_diabetes: corr(glucose, blood formate) = %+.2f' % r)

cm = CoupledModel(beta=bc)
dfo = cm.scenarios()
dfo.to_csv('data/obesity_diabetes.csv', index=False)
print('obesity_diabetes: formate vs healthy %',
      [round(v, 1) for v in dfo['formate_vs_healthy_%']])
