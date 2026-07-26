# Literature: data to test the formate → insulin model

Compiled 2026-07-24 from a four-track literature search. Each entry maps to a
specific arm of the β-cell model (`beta_cell.py`) and is tagged **supports**,
**calibrates**, **challenges**, or **gap**. Identifiers are given for the key
papers; a few flagged "provisional" were not opened in full text.

## Verdict

Several datasets can test or calibrate specific arms **now**, and the
mTOR→insulin arm is strongly supported. But the exact chain *mitochondrial
formate → purine → mTOR → insulin* has **never been measured end-to-end in a
β-cell**, and two population-level findings currently point the **wrong way**
(circulating formate ↑ with incident T2D; uric acid ↑ with T2D). Both are
addressed below.

## Evidence by model arm

| Model link | Key data | Tag |
|---|---|---|
| Glucose → ATP/ADP | Detimary 1998 *JBC* 273:33905 — ATP/ADP **2.4→11.6** (1→10 mM); Detimary 1996 *JBC* 271:20559 — GTP & GTP/GDP rise with glucose | **calibrates** |
| Formate → de novo purine | Human tracer PMID 17445548 — oral [¹³C]formate → purine C2/C8 → uric acid; Meiser 2016 formate overflow *Sci Adv* 2:e1601273 | supports (flux); **gap** in β-cell |
| Purine → mTORC1 | Hoxhaj 2017 *Cell Rep* 21:1331 (purine, not pyrimidine, depletion inhibits mTORC1; TSC/Rheb, AMPK-independent); Emmanuel 2017 *Cell Rep* 19:2665 (GTP-loading of Rheb) | supports in general cells; **gap** in β-cell |
| mTORC1 → insulin synthesis / content | Blandino-Rosano 2017 *Nat Commun* 8:16014 (mTORC1→4E-BP2/eIF4E→CPE; ↑proinsulin/insulin ratio); **Ni 2017 *Nat Commun* 8:15755 — β-cell Raptor-KO insulin content −50%**, mass −50%, GSIS −40%; Pende 2000 *Nature* 408:994 (S6K1−/−); Rachdi 2008 *PNAS* 105:9250 & Hamada 2009 *Diabetes* 58:1321 (Tsc2/Rheb gain) | **supports + calibrates** |
| Glucose(ATP/ADP) → secretion | Ashcroft 1984 *Nature* 312:446 (K-ATP) | supports (canonical) |
| Purine pools ↔ secretion | Detimary 1998 (ATP/ADP, GTP); **Gooding 2015 *Cell Rep* 13:157 — glucose drops IMP −77%, raises S-AMP 3.4×; S-AMP infusion → 4.1× exocytosis**; Kibbey 2007 *Cell Metab* 5:253 (mtGTP); Kowluru 2020 review | **calibrates**; note route is K-ATP/exocytosis, not shown via mTOR |
| Purine degradation → urate | Human tracer PMID 17445548 (formate C → uric acid) | supports (flux); see urate tension |
| Mito-formate source / deficiency | MTHFD1L & MTHFD2 expressed in human islet endocrine cells (Human Protein Atlas); Momb 2013 *PNAS* 110:549 (MTHFD1L = formate source, formate-rescuable, embryo-lethal); **Hsu 2013 *PLoS ONE* 8:e77931 — folate deficiency abolishes β-cell insulin biosynthesis + GSIS**; Karampelias 2021 *Nat Commun* 12:3362 (folate→β-cell, methotrexate blocks) | **supports** (deficiency arm) |
| AMPK opposition | da Silva Xavier 2003 *Biochem J* 371:761; Nguyen-Tu 2022 *Diabetologia* 65:997 (chronic AMPK ↓ secretion) | supports (chronic); "formate→AMPK" untested |
| Rapamycin/sirolimus → ↓insulin, diabetes | Barlow 2013 *Diabetes* 62:2674; Lombardi 2017 *Sci Rep* 7:15823; Granata 2023 clinical PTDM | supports mTOR→insulin |
| Whole-body: formate/serine → glucose+insulin (CENTRAL) | **Benatar 2026 *FASEB J* — O-acetyl-serine → +3.7× insulin at 15 min AND −16% OGTT AUC** (matches central signature); Holm 2018 *PLoS ONE* 13:e0194414 (serine → ↓GTT AUC, ↓HOMA-IR, NOD); Altaweel 2009 *J Ocul Pharmacol Ther* (human formate PK, baseline 0.57 mM) | supports central arm |

## Three tensions to confront

**1. A large prospective cohort points the wrong way.** Takase et al.,
*Cardiovasc Diabetol* 2025;24:335 (TMM CommCohort, n=12,461, 354 incident T2D,
4.3 y, non-fasting NMR): **formic acid RR 1.45 (1.11–1.91), P=0.023** per 1-SD —
**higher** formate predicts **more** incident T2D, opposite the core prediction.
Model's main support (Pietzke 2019 *Cancer Metab* 7:3, low formate in
obese/diabetic) is small and retrospective; no independent cohort replicates it.

*Reconciliation (modelled — see `formate_diabetes.py`).* The model has two routes
to diabetes: (A) formate deficiency → low insulin → LOW formate (Pietzke); (B) a
primary insulin-**secretion** defect (the paper's polygenic risk), independent of
formate. Under (B), less insulin → less insulin-driven de novo purine synthesis,
the main formate **sink** → formate **overflows** upward. Making blood formate a
dynamic pool (dietary intake − insulin-driven consumption − renal loss) and
imposing a secretion defect, the model reproduces Takase: as secretion is
impaired, glucose rises **4.8→7.8 mM** and blood formate rises **×1.5–1.6 at
every dietary intake** (0.8→1.3, 1.6→2.6, 4.4→6.4 mg/L for low/normal/high
serine-formate intake) — a positive formate–glucose association from reduced
*utilisation*, not causation. Consistent with the cohort's own data: formate did
**not** differ by genetic risk (11.3 vs 11.2, P=0.77) and was **not** a PRS→T2D
mediator, i.e. it behaves as a downstream consequence marker, exactly route (B).

**2. Knocking the pathway down *increases* secretion.** Pelligra 2023
*Cell Rep* 42:112615: shRNA of SHMT2/MTHFD2 in islets **increased** GSIS, yet
ATF4-driven induction of the same pathway **raised insulin content while lowering
secretion**. This is the *synthesis-vs-secretion dissociation* the model is built
around: the model predicts formate drives *content/synthesis*, and it already
lets acute secretion move independently → the implication is **measure insulin
content / proinsulin translation, not just GSIS**. They never tested MTHFD1L or
measured formate/purines (the missing experiment).

**3. Uric acid runs opposite — but is not causal.** High serum urate predicts
more T2D (RR 1.56; Lv 2013 *PLoS ONE* 8:e56864), tensioning "formate→urate→
healthy". Mendelian randomization is **null** (Pfister 2011 *Diabetologia*
54:2561; Sluijs 2015 *Diabetes* 64:3028): urate is a *marker*, not a cause. Fix:
treat urate in the model as a **flux readout of purine turnover** (tracer-
supported), not a health claim.

Also: glycine *worsens* glucose in obese mice via gluconeogenesis (Alves 2022
*Nutrients* 15:96) — restrict claims to serine/formate, not "one-carbon donors."

## Datasets usable to calibrate the model NOW

- **Energy module** — Detimary 1998: ATP/ADP **2.4 → 11.6** (1→10 mM glucose),
  β-cell-specific; Detimary 1996: GTP/GDP rises in parallel. (Model previously
  under-swung: 2.9→6.8.)
- **mTOR → insulin content** — Ni 2017: β-cell Raptor-KO → insulin content
  **−50%** (mTORC1 loss halves, not abolishes, content → implies an
  mTOR-independent basal synthesis fraction).
- **IMP branch** — Gooding 2015: high glucose drops **IMP −77%** and raises
  **S-AMP (adenylosuccinate) 3.4×**.
- **GSIS** — EC50 ≈ 8 mM, stimulation index ≈ 5–6 (already used).
- **Human formate exposure** — Altaweel 2009: baseline serum formate 0.57 mM;
  3.9 g dose → 30–90 µM rise (for the supplementation arm).

## Decisive experiments the model motivates (none exist)

1. **β-cell-specific MTHFD1L conditional KO** (Ins1-Cre × Mthfd1l^fl): measure
   intracellular formate, purine pools, proinsulin/*INS* translation, insulin
   content, and GSIS. Highest-value test.
2. **Formate/serine supplementation GTT with paired plasma insulin AND urate** —
   runs the `whole_body.py` central-vs-peripheral discriminator directly
   (Benatar's OAS result already leans central).
3. **Does purine depletion inhibit mTORC1 in islets?** (lometrexol/mycophenolate
   + p-S6/p-4E-BP) — the one link with zero β-cell evidence.

## Gaps (explicit "no data found")

- No formate-dosing study with insulin/GTT endpoints.
- No serine/formate study reporting plasma urate as an outcome.
- No MTHFD1L knockdown/KO in any β-cell line/islet with an insulin phenotype.
- No purine→mTORC1 demonstration in β-cells.
- No MTHFD1L GWAS/eQTL signal for T2D (MTHFD1 R653Q is linked; MTHFD1L is not).
