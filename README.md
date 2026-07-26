# A kinetic model of insulin synthesis in the pancreatic β-cell

Tests the *Spice of Life* hypothesis (Ch. 7): **formate stimulates insulin
synthesis**. Proposed chain:

```
formate → 10-CHO-THF → de novo purine synthesis → adenine + guanine nucleotides
        → mTORC1 activation [Hoxhaj 2017] → translation of the insulin peptide
```

`beta_cell.py` is a genuine kinetic model — every species has a `dX/dt` built
from mass-action / Michaelis–Menten / Hill rate laws. It can be **integrated in
time** (`simulate`) or solved for its **steady state** (`steady_state`).
Time unit: hours; concentrations: mM (insulin pools in arbitrary granule units).

## The five coupled modules

| # | Module | Species | Reused from |
|---|--------|---------|-------------|
| 1 | Energy metabolism | ATP, ADP, AMP | Formate-NUDT5 energy core [PUB] |
| 2 | Purine synthesis | formate, 10-CHO-THF, PRPP, GAR, AICAR, IMP | Formate-NUDT5 purine block + glue [PUB]/[FIT] |
| 3 | Purine degradation | IMP→hypoxanthine/xanthine (hx) | new [LIT]/[ASM] |
| 4 | Uric acid secretion | urate | new [LIT]/[ASM] |
| 5 | Protein synthesis | pro-insulin, insulin | new [LIT]/[ASM] |

Glucose is the **external input**; blood formate `form_x` is the **lever**.

## Signal flow

- **Glucose sensing.** Glucokinase gate `gk = G^1.7/(8^1.7+G^1.7)` sets ATP
  production (glycolysis + OxPhos). The `ATP/ADP` ratio is *both* the anabolic
  fuel and, through a Hill trigger, the K‑ATP/Ca²⁺ **secretion** signal.
- **Formate → purines.** FTHFS turns formate into 10‑CHO‑THF; PPAT/GART/ATIC
  (the same rate laws as `Formate-NUDT5/switch_ppat.py`, including the optional
  AMP/PRPP‑gated PPAT–NUDT5 glue) build IMP, which branches to AMP and GMP.
- **Purines → mTOR.** `mtor = f(ATP+GXP) · g(energy charge)` — the Hoxhaj 2017
  purine-sensing arm, gated by the adenylate energy charge (AMPK proxy).
- **mTOR → insulin.** Pro-insulin translation `v_ins_syn = k·mtor·g_tx(G)·aa`,
  maturation, then glucose-triggered secretion.
- **Degradation → urate.** IMP/AMP/GXP → hx → (xanthine oxidase) → urate →
  secreted. Links purine turnover to the book's uric-acid formate marker.

## Running

```
python beta_cell.py       # calibrate + absolute units + steady states + formate scan + deficiency + sensitivity
python whole_body.py      # 2-compartment loop: central vs peripheral formate discrimination
python formate_diabetes.py# secretion-defect perturbation: blood formate vs diabetes x diet (Takase 2025)
python adipocyte.py       # standalone adipocyte metabolic model (formate sink + XOR urate source)
python obesity_diabetes.py# 3-compartment (pancreas+adipose+blood): blood formate/urate vs obesity x diabetes;
                          #   .intervention() = formate/serine supplementation responder analysis;
                          #   .metformin()    = metformin (peripheral + SHMT2/formate) dose response
python eye.py             # eye compartment: local retinal insulin set by blood one-carbon,
                          #   decoupled from plasma insulin (blood-retinal barrier);
                          #   .methanol_toxicity() = biphasic formate dose (spice of life -> complex-IV poison)
python make_figs.py       # fig1..fig5 (calibration, formate scan, deficiency, gtt, whole-body)
```
Outputs are written to `data/*.csv`. Progress write-up: `manuscript.tex` →
`manuscript.pdf` (build: `latexmk -pdf manuscript.tex`); `manuscript.md` is the
earlier plain-text draft.

## Modules & key functions

- `recalibrate_all()` — runs the four data-anchored calibrations in order:
  - `calibrate_energy()` → ATP/ADP to **Detimary 1998** (2.4→11.6 over 1→10 mM).
  - `calibrate_imp_branch()` → IMP/S-AMP to **Gooding 2015** (IMP −77%, S-AMP +3.4×).
  - `calibrate_mtor_content()` → mTOR-independent fraction `B_MTOR` to **Ni 2017**
    (Raptor-KO insulin content −50%).
  - `calibrate()` → glucose→insulin axis (`KG_tx`, `B_TX`) to GSIS (EC50 8 mM, SI 5.5).
- `steady_state()`, `simulate()` — steady-state and time-domain solvers
  (both take `mtor_mult` for the Raptor-KO test).
- `absolute_units()` — per-cell insulin content/secretion (anchored to ~20 pg/cell).
- `formate_scan()`, `formate_deficiency_scan()` — hypothesis and deficiency tests.
- `sensitivity()` — 2× sweep of `[ASM]` params, reporting SI and the formate
  sensitivity of insulin synthesis (FS).

## Key results (data-anchored calibrated model)

Calibration targets (all hit): ATP/ADP **2.4→11.6** at 1→10 mM (Detimary 1998,
ATP 0.6→2.3 mM, physiological); IMP **×0.25** and S-AMP **×3.4** low→high glucose
(Gooding 2015); Raptor-KO content **50%** of WT (Ni 2017); GSIS **EC50 8.0 mM,
SI 5.5**.

- **Formate stimulates insulin synthesis:** at 11 mM glucose, raising blood
  formate 0→0.3 mM raises mTOR and insulin synthesis **~+59%** (sensitive across
  the physiological ≤0.3 mM range); urate rises in parallel.
- **Mitochondrial-formate deficiency impairs insulin:** cutting the MTHFD1L
  source 100%→0% (8 mM glucose) collapses mTOR (0.39→0.03) and lowers insulin
  synthesis by **~46 %** — now *bounded by the Ni-anchored mTOR-independent
  synthesis floor* `B_MTOR` (more conservative than the pre-calibration 88%).
  Urate falls 0.12→0.02 mM/h (the formate–urate seesaw). Dietary formate rescues.
- **Formate priming improves glucose tolerance** (GTT transient) and the
  **whole-body loop** discriminates central from peripheral action (see below).
- **Sensitivity:** both outputs robust to all insulin-secretion kinetic
  constants; SI set by energy/mTOR params, FS by purine-synthesis/mTOR params.

Data landscape and the experiments that would test the model: `literature.md`.
- **Absolute units:** anchoring the granule pool to β-cell insulin content
  (~20 pg/cell) gives stimulated secretion ~16 pg/cell/h (~2.8 fmol/cell/h);
  the scale leaves all calibrated ratios unchanged.
- **Central vs peripheral discrimination (`whole_body.py`):** in a Bergman-type
  glucose–insulin loop, central formate (β-cell synthesis) and peripheral
  formate (Carpéné: glucose effectiveness Sg) both lower glucose AUC ~19%, but
  **only central raises plasma insulin** (+58% vs −2%). → the discriminating
  experiment is to measure plasma insulin during a formate-supplemented GTT.

## Provenance & status

Parameters are tagged `[PUB]` (Formate-NUDT5/Oizel 2020 core), `[FIT]` (Witus
2026 glue; GSIS calibration), `[LIT]` (β-cell physiology / purine chemistry) or
`[ASM]` (assumed). The energy/purine core and the glucose→insulin axis are now
calibrated; the insulin/degradation kinetics remain scaled (absolute insulin
fluxes are arbitrary units). Next: measure the two most sensitive parameters
(mTOR purine-sensing Hill coefficient; mitochondrial formate capacity), and
discriminate the synthesis vs peripheral-uptake alternative experimentally.

## The synthesis vs peripheral-uptake question (now addressed)

The β-cell model routes the effect through **synthesis** (mTOR → translation →
granule refilling), the book's claim. The alternative — formate/methylamine
acting on **peripheral** glucose uptake (Carpéné et al. 2019) — is now
represented in `whole_body.py` as an increase in glucose effectiveness Sg. The
loop shows the two are experimentally separable by plasma insulin (above). See
`00_hypothesis_brief.md` for the original framing.
