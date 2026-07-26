# A kinetic model couples formate to insulin synthesis in the pancreatic β-cell

*Working manuscript — progress log, formatted for Molecular Systems Biology.
Details (absolute calibration, statistics, full parameter provenance) are
provisional at this stage.*

**Author:** Alexei Vazquez

---

## Abstract

The dietary one-carbon unit formate is required for de novo purine synthesis,
and cellular purine nucleotide levels are sensed by mTORC1, the master regulator
of anabolic protein synthesis. Insulin is a secreted 51-residue peptide whose
production is among the most demanding anabolic tasks of the pancreatic β-cell.
We previously proposed (*The Spice of Life*, 2020) that formate stimulates
insulin synthesis, providing a mechanistic link between one-carbon nutrition and
glucose homeostasis. Here we formalise this hypothesis as a kinetic model of the
β-cell that couples glucose-driven energy metabolism, de novo purine synthesis
and degradation, uric acid secretion, and mTOR-controlled translation of the
insulin peptide. The purine module reuses our published formate/energy model
(Formate-NUDT5). After calibrating the glucose–insulin axis to islet-perifusion
benchmarks (EC50 ≈ 9 mM, stimulation index ≈ 6.6), the model predicts that (i)
raising formate increases the purine nucleotide pool, mTOR activity and insulin
synthesis, with concomitant uric acid secretion; (ii) loss of the mitochondrial
formate source (the MTHFD1L pathway) causes an ~88 % fall in insulin synthesis
and a parallel drop in uric acid — reproducing the low-formate, low-urate state
associated with diabetes; and (iii) dietary formate rescues insulin synthesis on
a formate-deficient background. A global sensitivity analysis shows the
formate→insulin coupling is controlled specifically by purine-synthesis and
mTOR-sensing parameters and is robust to the assumed insulin-secretion kinetics.
The model turns a verbal hypothesis into quantitative, testable predictions.

## Introduction

Formate is the circulating carrier of one-carbon units. It is produced
endogenously by the mitochondrial serine catabolism pathway (SHMT2–MTHFD2–
MTHFD1L) and supplied exogenously from diet and gut microbiota. Two carbons of
every purine ring derive from formate; formate limitation therefore constrains
de novo purine synthesis (Oizel et al, 2020). Independently, mTORC1 was shown to
sense cellular purine nucleotide levels (Hoxhaj et al, 2017). Insulin is a
secreted peptide hormone; its synthesis is an intense, mTOR-dependent anabolic
program triggered acutely by glucose. These facts suggest a chain —
**formate → purine nucleotides → mTORC1 → insulin peptide synthesis** — that
would make formate a nutritional modulator of β-cell function and, by extension,
of diabetes risk. The chain has not been tested directly. We build a kinetic
model to determine whether it is quantitatively coherent and to expose its
critical, experimentally addressable parameters.

## Results

### A kinetic β-cell model coupling formate, purines, mTOR and insulin

The model (`beta_cell.py`) is a system of ordinary differential equations for 14
species in a single β-cell, with glucose as the external input and blood formate
and mitochondrial formate production as levers. Five modules are coupled
(Fig 1A, schematic):

1. **Energy metabolism.** Glucose is committed by a glucokinase-like sensor
   (Hill, half-saturation 8 mM) feeding glycolytic and oxidative ATP production;
   ATP, ADP and AMP are linked by adenylate kinase. The ATP/ADP ratio is both
   the anabolic fuel and, through a Hill trigger, the K-ATP/Ca²⁺ secretion
   signal.
2. **Purine synthesis.** The formate → 10-CHO-THF → PRPP → PPAT → GART → ATIC
   block of the Formate-NUDT5 model, including the optional AMP/PRPP-gated
   PPAT–NUDT5 metabolite glue, terminating in IMP.
3. **Purine degradation.** IMP, AMP and guanine nucleotides converge on
   hypoxanthine/xanthine and are oxidised to urate.
4. **Uric acid secretion.** Urate is exported from the cell.
5. **Insulin synthesis.** mTORC1 activity is a Hill function of the purine
   nucleotide pool (weighted to the guanine-nucleotide, de-novo-controlled pool)
   gated by the adenylate energy charge; it drives translation of pro-insulin,
   which matures and is secreted under the glucose trigger.

Every species has an explicit rate law (mass-action / Michaelis–Menten / Hill).
The model is integrated in time to simulate transients and solved at steady
state for dose-response analysis.

### Calibration to islet-perifusion GSIS benchmarks

We calibrated the glucose–insulin axis to standard dynamic-perifusion benchmarks
of the concentration-dependent (second-phase) response: a sigmoidal dependence
on glucose with half-maximal secretion near 8 mM (human islets 7.9 mM;
PMC6783504) and a stimulation index of ~5–6 between basal (2.8 mM) and
stimulatory (16.7 mM) glucose (Nunemaker et al, 2006). Fitting the glucose
half-saturation and basal drive of the insulin gene reproduced these targets
(**EC50 = 8.9 mM, stimulation index = 6.6**; Fig 1B). The transient first phase
(peak ~1 min) is below the model's time resolution and was not fitted. At steady
state the calibrated model gives physiological adenine-nucleotide levels and
ATP/ADP ratios rising from 2.9 (2.8 mM glucose) to 6.8 (16.7 mM).

### Formate stimulates insulin synthesis

At stimulatory glucose (11 mM) on a formate-limited background, raising blood
formate from 0 to 0.3 mM increased the purine nucleotide pool, mTOR activity
(0.09 → 0.31) and the insulin-synthesis rate ~3.5-fold, saturating thereafter
(Fig 2). Uric acid secretion rose in parallel, from 0.029 to 0.088 mM/h. The
response is steep across the physiological blood-formate range (≤0.3 mM),
supporting the hypothesis that formate availability sets an anabolic ceiling on
insulin production.

### Mitochondrial formate deficiency impairs insulin (formate theory of diabetes)

We simulated formate deficiency by suppressing the endogenous mitochondrial
formate source (the MTHFD1L pathway). Lowering mitochondrial formate production
from normal to zero, at fixed 8 mM glucose, reduced intracellular formate,
collapsed mTOR activity (0.24 → 0.03) and cut insulin synthesis by **88 %**
(Fig 3). Urate secretion fell in parallel (0.089 → 0.018 mM/h). The model thus
predicts that a defect in endogenous formate production — genetic (MTHFD1L
variants), mitochondrial, or environmental (antibiotic depletion of formate-
producing microbiota) — produces a coupled low-insulin, low-urate state,
consistent with the reduced circulating formate reported in diabetes and its
seesaw with glucose (Pietzke et al, 2019). Providing dietary formate rescued
insulin synthesis on the deficient background (Fig 2), the in-silico analog of
sodium-formate supplementation.

In a glucose-tolerance transient, a formate-replete cell entered the glucose
bolus with higher basal mTOR and a fuller granule pool and mounted a roughly
two-fold larger integrated insulin response than a formate-deficient cell,
linking chronic formate status to acute glucose tolerance.

### Sensitivity analysis isolates the controlling parameters

Varying each assumed parameter two-fold (Fig 4 / `data/sensitivity.csv`), we
tracked two outputs: the GSIS stimulation index (SI) and the formate sensitivity
of insulin synthesis (FS, the ratio of insulin synthesis at replete vs deficient
formate). Both outputs were insensitive (relative change < 1 %) to every
insulin-secretion kinetic constant (maturation, secretion rate, secretion
trigger, ATP cost). SI was controlled by the energy/mTOR parameters (basal fuel
floor, purine-sensing Hill coefficient and half-saturation), and FS was
controlled specifically by the purine-synthesis and mTOR-sensing parameters
(purine turnover scale, mitochondrial formate capacity, the IMP-branch and
AMP-degradation rates, and the purine-sensing Hill coefficient). The
formate→insulin prediction therefore rests on the identifiable purine/mTOR
module and not on the loosely-constrained secretion kinetics.

## Discussion

A kinetic model shows that the proposed formate → purine → mTOR → insulin chain
is quantitatively coherent and, once the glucose–insulin axis is calibrated to
islet data, yields non-trivial predictions: a formate dose-response of insulin
synthesis, an 88 % insulin deficit under mitochondrial formate loss, coupled
uric-acid changes, and dietary rescue. The sensitivity analysis pinpoints the
mTOR purine-sensing Hill coefficient and the mitochondrial formate capacity as
the parameters most in need of measurement.

Two caveats bound the current model. First, parameters outside the reused
energy/purine core and the GSIS calibration are scaled, not fitted; absolute
insulin fluxes are in arbitrary units. Second, and most important, the model
routes the entire effect through β-cell **synthesis**. An alternative in which
formate (or its precursor methylamine) acts on **peripheral** glucose uptake
without changing β-cell insulin (Carpene et al, 2019) is not represented and
must be discriminated experimentally — e.g. by measuring proinsulin/insulin
content and INS translation in islets across a formate gradient, and GSIS in
islets from MTHFD1L-deficient models.

## Model and Methods (brief)

Fourteen ODEs (ATP, ADP, AMP, formate, 10-CHO-THF, PRPP, GAR, AICAR, IMP,
guanine nucleotides, hypoxanthine/xanthine, urate, pro-insulin, insulin) with
mass-action/MM/Hill rate laws; time unit hours, concentrations mM (insulin in
arbitrary units). Integrated with an implicit solver (BDF); steady states from a
pre-integration followed by a log-space least-squares polish (residuals < 1e-13
mM/h). The purine/energy core and the PPAT–NUDT5 glue calibration are inherited
from Formate-NUDT5 (Vazquez, submitted); GSIS parameters were fitted here by
least squares to EC50 and stimulation-index targets. Parameters carry provenance
tags [PUB]/[FIT]/[LIT]/[ASM]. Code, calibrated outputs (`data/*.csv`) and
figures are in this repository (`beta_cell.py`, `make_figs.py`).

## Key numbers (this run)

| Quantity | Value |
|---|---|
| GSIS EC50 (calibrated) | 8.9 mM (target 8.0) |
| GSIS stimulation index | 6.6 (target ~5.5) |
| Insulin synthesis, formate 0→0.3 mM (11 mM glucose) | ~3.5× |
| Insulin synthesis loss, mito-formate 100%→0% | 88 % |
| mTOR, mito-formate 100%→0% | 0.24 → 0.03 |
| Urate secretion, mito-formate 100%→0% | 0.089 → 0.018 mM/h |
| GTT integrated insulin, replete vs deficient | ~2× |

## References

- Hoxhaj G et al (2017) The mTORC1 signaling network senses changes in cellular
  purine nucleotide levels. *Cell Rep* 21:1331.
- Oizel K et al (2020) Formate induces a metabolic switch in nucleotide and
  energy metabolism. *Cell Death Dis* 11:310.
- Pietzke M et al (2019) Stratification of cancer and diabetes based on
  circulating levels of formate and glucose. *Cancer Metab* 7:3.
- Carpéné C et al (2019) Methylamine activates glucose uptake in human
  adipocytes without overpassing action of insulin or stimulating its secretion
  in pancreatic islets. *Medicines* 6:89.
- Nunemaker CS et al (2006) Insulin secretion in the conscious mouse is biphasic
  and pulsatile. *Am J Physiol Endocrinol Metab*; and In vivo/in vitro biphasic
  GSIS, *Diabetes* 55:441.
- Vazquez A (2020) *The Spice of Life*.
- Vazquez A, An allosteric metabolite glue couples AICAR inversely to one-carbon
  availability (submitted) — Formate-NUDT5.
