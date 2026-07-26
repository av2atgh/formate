#!/usr/bin/env python3
"""Assemble plos_manuscript.tex from manuscript.tex in PLOS Comp Biol format:
new PLOS preamble + title/abstract/author-summary front matter, body with
author-year citations converted to numbered \\cite{}, PLOS section headings,
and end matter (Supporting information, Acknowledgments, Data availability,
Competing interests, \\bibliography). The original manuscript.tex is untouched."""
import re

src = open('manuscript.tex').read()

# ---- extract the body: Introduction .. end of the Key-numbers table ----
b0 = src.index(r'\section{Introduction}')
b1 = src.index(r'\end{table}', src.index('Key numbers')) + len(r'\end{table}')
body = src[b0:b1]

# ---- citation conversions (author-year -> numbered \cite) ----------------
subs = [
    # grouped citations first
    (r"\(DeFronzo,\s*2009;\s*Roden\s*\\&\s*Shulman,\s*2019\)", r"\\cite{defronzo2009,roden2019}"),
    (r"\(Ducker\s*\\&\s*Rabinowitz,\s*2017;\s*Pietzke,\s*Meiser\s*\\&\s*Vazquez,\s*2020\)", r"\\cite{ducker2017,pietzke2020}"),
    (r"\(Ruberte\s*2004;\s*Villacampa\s*2013\)", r"\\cite{ruberte2004,villacampa2013}"),
    # single parentheticals (full author-year)
    (r"\(Rorsman\s*\\&\s*Ashcroft,\s*2018\)", r"\\cite{rorsman2018}"),
    (r"\(Itoh\s*\\&\s*Okamoto,\s*1980\)", r"\\cite{itoh1980}"),
    (r"\(Saxton\s*\\&\s*Sabatini,\s*2017\)", r"\\cite{saxton2017}"),
    (r"\(Kahn,\s*Hull\s*\\&\s*Utzschneider,\s*2006\)", r"\\cite{kahn2006}"),
    (r"\(Oizel\s+\\textit\{et al\},\s*2020\)", r"\\cite{oizel2020}"),
    (r"\(Brosnan\s+\\textit\{et al\},\s*2018\)", r"\\cite{brosnan2018}"),
    (r"\(Pietzke\s+\\textit\{et al\},\s*2019\)", r"\\cite{pietzke2019}"),
    (r"\(Hoxhaj\s+\\textit\{et al\},\s*2017\)", r"\\cite{hoxhaj2017}"),
    (r"\(\\textit\{The Spice\s+of Life\},\s*2020\)", r"\\cite{spiceoflife}"),
    (r"\(\\textit\{The Spice of Life\},\s*p\.~53\)", r"\\cite{spiceoflife}"),
    (r"\(Detimary\s+\\textit\{et al\},\s*1998\)", r"\\cite{detimary1998}"),
    (r"\(Detimary\s*1996\)", r"\\cite{detimary1996}"),
    (r"\(Ni\s+\\textit\{et al\},\s*2017\)", r"\\cite{ni2017}"),
    (r"\(Nunemaker\s+\\textit\{et al\},\s*2006\)", r"\\cite{nunemaker2006}"),
    (r"\(Carp\\'ene\s+\\textit\{et al\},\s*2019\)", r"\\cite{carpene2019}"),
    (r"\(Tsushima\s+\\textit\{et\s+al\},\s*2013\)", r"\\cite{tsushima2013}"),
    (r"\(Kasahara\s+\\textit\{et al\},\s*2023\)", r"\\cite{kasahara2023}"),
    (r"\(Vazquez,\s*2020\)", r"\\cite{vazquez2020urate}"),
    (r"\(Etchegaray\s+\\textit\{et al\},\s*2023\)", r"\\cite{etchegaray2023}"),
    (r"\(Ben-Sahra\s+\\textit\{et al\},\s*2016\)", r"\\cite{bensahra2016}"),
    (r"\(Pelligra\s+\\textit\{et al\},\s*2023\)", r"\\cite{pelligra2023}"),
    (r"\(Leontieva\s+\\textit\{et al\},\s*2014\)", r"\\cite{leontieva2014}"),
    (r"\(Ruberte\s+\\textit\{et al\},\s*2004\)", r"\\cite{ruberte2004}"),
    (r"\(Villacampa\s+\\textit\{et al\},\s*2013\)", r"\\cite{villacampa2013}"),
    (r",\s*Altaweel\s+\\textit\{et al\},\s*2009", r"~\\cite{altaweel2009}"),
    (r";\s*Ding\s+\\textit\{et al\},\s*2014", r"~\\cite{ding2014}"),
    (r"\(Imamura\s+\\textit\{et al\},\s*2015\)", r"\\cite{imamura2015}"),
    (r"\(Li\s+\\textit\{et al\},\s*2024\)", r"\\cite{li2024}"),
    (r";\s*Delpino\s*\\&\s*Figueiredo,\s*2022", r"~\\cite{delpino2022}"),
    (r"\(Tramonti\s+\\textit\{et al\},\s*2021\)", r"\\cite{tramonti2021}"),
    (r"\(Fraenkel\s+\\textit\{et al\},\s*2008\)", r"\\cite{fraenkel2008}"),
    (r"\(\\textit\{Exp Eye Res\},\s*2025\)", r"\\cite{yao2025}"),
    # narrative citations (name stays, year -> \cite)
    (r"Gooding\s+\\textit\{et al\}\s*\(2015\)", r"Gooding et al.~\\cite{gooding2015}"),
    (r"Takase\s+\\textit\{et al\},\s*2025;", r"Takase et al.~\\cite{takase2025};"),
    (r"Takase\s+\\textit\{et al\}\s*\(2025\)", r"Takase et al.~\\cite{takase2025}"),
    (r"Pietzke,\s*Meiser\s*\\&\s*Vazquez\s*\(2020,\s*Table~2\)", r"Pietzke et al.~\\cite{pietzke2020} (Table~2)"),
    (r";\s*Kawai\s+\\textit\{et al\},\s*2019", r"~\\cite{kawai2019}"),
    (r"\(Vazquez,\s*submitted\)", r"\\cite{vazqueznudt5}"),
    (r"mechanism of \\textit\{The Spice\s+of Life\}", r"mechanism of \\textit{The Spice of Life}~\\cite{spiceoflife}"),
    # short caption/table forms
    (r"\(Pietzke\s*2019\)", r"\\cite{pietzke2019}"),
    (r"Pietzke~2020", r"\\cite{pietzke2020}"),
    (r"\(Takase\s*2025\)", r"\\cite{takase2025}"),
    (r"\(Tsushima\s*2013\)", r"\\cite{tsushima2013}"),
    (r"\(Hsu\s*2013\)", r"\\cite{hsu2013}"),
    (r"\(Pelligra\s*2023\)", r"\\cite{pelligra2023}"),
    (r"\(Detimary\s*1998\)", r"\\cite{detimary1998}"),
    (r"\(Gooding\s*2015\)", r"\\cite{gooding2015}"),
    (r"\(Ni\s*2017\)", r"\\cite{ni2017}"),
    (r";\s*Pietzke\s+2019", r"~\\cite{pietzke2019}"),
    (r";\s*Takase\s+2025", r"~\\cite{takase2025}"),
    (r"reading of Pietzke\s+2019", r"reading of Pietzke et al.~\\cite{pietzke2019}"),
    # ground the "mTOR->insulin arm strongly supported" claim (cites Blandino-Rosano)
    (r"\(Raptor/\\allowbreak TSC2 \$\\beta\$-cell models; rapamycin\)",
     r"(Raptor/\\allowbreak TSC2 $\\beta$-cell models; rapamycin)~\\cite{ni2017,blandino2017,fraenkel2008}"),
    # supporting-information references
    (r"the Supplementary Material gives", r"S1 Appendix gives"),
    (r"\(Supplementary Fig\.~S1\)", r"(Fig A in S1 Appendix)"),
    (r"\(Fig\.~S1\)", r"(Fig A in S1 Appendix)"),
    (r"\s*\(supplementary\s+\\texttt\{lit\\-er\\-a\\-ture\.md\}\)", r""),
    # section cross-references (starred sections are unnumbered -> reword)
    (r"\(\\S\\ref\{ss:dietretina\}\)", r"(below)"),
    (r"\$\\S\$\\ref\{ss:eyecomp\}", r"the eye-compartment model"),
    (r"\(\$\\S\$\\ref\{ss:suppl\}\)", r"(above)"),
    (r"\\label\{ss:eyecomp\}", r""),
    (r"\\label\{ss:suppl\}", r""),
    (r"\\label\{ss:dietretina\}", r""),
    # PLOS style: "Fig" not "Figure."
    (r"Fig\.~\\ref", r"Fig~\\ref"),
]
for pat, rep in subs:
    body, n = re.subn(pat, rep, body)
    if n == 0:
        print('WARN no match:', pat)

# tie preceding space to \cite ("word \cite" -> "word~\cite")
body = re.sub(r"(?<=\S) \\cite\{", r"~\\cite{", body)

# unnumbered PLOS headings + rename Methods section
body = body.replace(r'\subsubsection{', r'\subsubsection*{')
body = body.replace(r'\subsection{', r'\subsection*{')
body = body.replace(r'\section{', r'\section*{')
body = body.replace(r'\section*{Model and Methods}', r'\section*{Materials and methods}')

# ---- new preamble (from PLOS template, essentials + booktabs) ------------
PREAMBLE = r"""\documentclass[10pt,letterpaper]{article}
\usepackage[top=0.85in,left=2.75in,footskip=0.75in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{changepage}
\usepackage{textcomp,marvosym}
\usepackage{cite}
\usepackage{nameref,hyperref}
\usepackage[right]{lineno}
\usepackage[nopatch=eqnum]{microtype}
\DisableLigatures[f]{encoding = *, family = * }
\usepackage[table]{xcolor}
\usepackage{array}
\usepackage{booktabs}
\newcolumntype{+}{!{\vrule width 2pt}}
\newlength\savedwidth
\newcommand\thickcline[1]{%
  \noalign{\global\savedwidth\arrayrulewidth\global\arrayrulewidth 2pt}%
  \cline{#1}%
  \noalign{\vskip\arrayrulewidth}%
  \noalign{\global\arrayrulewidth\savedwidth}%
}
\newcommand\thickhline{\noalign{\global\savedwidth\arrayrulewidth\global\arrayrulewidth 2pt}%
\hline
\noalign{\global\arrayrulewidth\savedwidth}}
\raggedright
\setlength{\parindent}{0.5cm}
\textwidth 5.25in
\textheight 8.75in
\usepackage[aboveskip=1pt,labelfont=bf,labelsep=period,justification=raggedright,singlelinecheck=off]{caption}
\renewcommand{\figurename}{Fig}
\bibliographystyle{plos2025}
\makeatletter
\renewcommand{\@biblabel}[1]{\quad#1.}
\makeatother
\usepackage{lastpage,fancyhdr,graphicx}
\usepackage{epstopdf}
\pagestyle{fancy}
\fancyhf{}
\rfoot{\thepage/\pageref{LastPage}}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrule}{\hrule height 2pt \vspace{2mm}}
\fancyheadoffset[L]{2.25in}
\fancyfootoffset[L]{2.25in}
\lfoot{}

\begin{document}
\vspace*{0.2in}

\begin{flushleft}
{\Large
\textbf{A computational model of the interactions between formate and insulin metabolism}
}
\newline
\\
Alexei Vazquez\textsuperscript{1*}
\\
\bigskip
\textbf{1} Nodes \& Links Ltd, Salisbury House, Station Road, Cambridge, CB1 2LA, United Kingdom
\\
\bigskip
* alexei.vazquez@gmail.com
\end{flushleft}

\section*{Abstract}
Formate, the circulating carrier of one-carbon units, supplies two carbons of
every purine ring, and cellular purine nucleotides are sensed by mTORC1 --- the
regulator of the anabolic translation the pancreatic $\beta$-cell uses to make
insulin. These links suggest an untested chain, formate $\rightarrow$ purine
nucleotides $\rightarrow$ mTORC1 $\rightarrow$ insulin synthesis, connecting
one-carbon nutrition to insulin production and diabetes risk. We formalise this
chain as a computational kinetic model calibrated to independent islet data; the
model itself --- a quantitative substrate for reasoning about formate and insulin
--- is the central contribution. Systemically, raising formate increases the
purine pool, mTOR activity and insulin synthesis, while loss of mitochondrial
formate production cuts insulin synthesis by roughly half with a parallel fall in
uric acid. Embedded in a whole-body glucose--insulin loop, this central action
raises plasma insulin, distinguishing it from a peripheral effect on glucose
uptake. An adipose compartment shows that obesity and an insulin-secretion defect
drive blood formate in opposite directions --- reconciling human cohorts that
disagree in sign --- while food formate content does not predict diabetes risk,
because the $\beta$-cell is already formate-replete. A gut-microbiome compartment
that returns urate carbon as formate closes a delayed positive-feedback urate
cycle. The retina is distinct: it makes its own insulin behind the blood-retinal
barrier, so retinal insulin tracks blood one-carbon rather than plasma insulin.
Self-regulation makes the eye a homeostat; because excess retinal insulin is
pathological, retinal health is an inverted-U in formate, and methanol optic
toxicity emerges as an acute one-carbon overload. The model turns these maps into
testable interventions: formate supplementation benefits only the
formate-deficient diabetic, and in the eye, inhibiting local formate synthesis
versus insulin synthesis is a paired, mechanism-discriminating prediction.

\section*{Author summary}
Insulin is made by specialised pancreatic cells, and too little of it causes
diabetes. I asked whether a small, common nutrient --- formate, the body's
carrier of one-carbon chemical groups --- helps set how much insulin these cells
can build. Formate supplies part of every purine, a building block of DNA and
energy molecules, and cells read their purine levels to switch on the growth
machinery (mTOR) that manufactures insulin. I built a computer model of this
chain, tuned to published measurements from insulin-producing cells, and used it
to explore what happens as formate rises and falls. The model reproduces puzzling
human data --- why blood formate can look high in some people with diabetes and
low in others --- and predicts who might benefit from formate supplements and who
would not. It also speaks to the eye, which makes its own insulin: too little or
too much formate harms the retina, and the model recasts methanol blindness as an
insulin-excess overload rather than direct poisoning of the cell's respiratory
chain. Most usefully, it proposes specific, testable experiments to confirm or
refute each prediction.

\clearpage
\newgeometry{top=0.85in,left=1in,right=1in,footskip=0.75in}
\linenumbers

"""

# ---- end matter ----------------------------------------------------------
ENDMATTER = r"""
\section*{Supporting information}

\paragraph*{S1 Appendix.}
\label{S1_Appendix}
\textbf{Supplementary methods and full model specification.} Complete
specification of every state variable, rate law, ODE balance and parameter value
with provenance; the calibration procedure; and the monostability analysis of the
regulated eye (basin tests and self-consistency/bifurcation analysis, Fig A).

\section*{Acknowledgments}
Nodes \& Links Ltd provided support in the form of salary for Alexei Vazquez, but
did not have any additional role in the conceptualization of the study, analysis,
decision to publish, or preparation of the manuscript. The code and text were
prepared with the assistance of Claude Opus~4.8, an AI assistant developed by
Anthropic.

\section*{Data availability}
All code underlying the results is openly available in the project repository:
\url{https://github.com/av2atgh/formate}. No separate data deposition is required
--- every quantity reported here, including the calibrated model outputs and all
figures, can be regenerated by running the code (the object-oriented model
\texttt{model.py} together with the figure generators \texttt{make\_*.py} and
\texttt{eye\_bifurcation.py}).

\section*{Competing interests}
AV is a paid employee of Nodes \& Links Ltd. This does not alter the author's
adherence to PLOS Computational Biology policies on sharing data and materials.

\nolinenumbers

\bibliography{formate_insulin}

\end{document}
"""

out = PREAMBLE + body + ENDMATTER
open('plos_manuscript.tex', 'w').write(out)
print('wrote plos_manuscript.tex  (%d chars)' % len(out))
# residual author-year check
for m in re.finditer(r".{0,40}\\textit\{et al\}[^)]{0,8}(19|20)\d\d.{0,4}", out):
    print('RESID et-al-year:', repr(m.group(0)))
for m in re.finditer(r".{0,32}(19|20)\d\d\)", out):
    print('RESID YYYY):', repr(m.group(0)))
