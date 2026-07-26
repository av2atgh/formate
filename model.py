#!/usr/bin/env python3
"""
Formalised computational model: formate -> purine -> mTOR -> insulin, across the
pancreatic beta-cell, adipose tissue, the whole body and the eye.

This module is the object-oriented formalisation of the model that is itself a
central contribution of the accompanying manuscript. Every compartment is a class
with an explicit parameter object, rate laws (``fluxes``), balance equations
(``derivatives``) and shared numerical machinery (``simulate`` / ``steady_state``)
inherited from a common :class:`MetabolicCell` base.

Class hierarchy
---------------
``MetabolicCell``            shared numerics (BDF integration, log-space steady
                             state) and rate-law primitives (Hill, Michaelis-
                             Menten, adenylate kinase, formate exchange).
``BetaCellModel``            15-state beta-cell (energy, one-carbon, de novo
                             purine, degradation/urate, mTOR-driven insulin).
``AdipocyteModel``           10-state adipocyte (shared core, but insulin-driven
                             GLUT4 uptake and lipid storage, XOR-rich urate).
``EyeModel``                 a retinal insulin cell (composes a BetaCellModel):
                             blood-retinal barrier, local-formate feed-forward
                             with S6K1 feedback, and the retinal-health window.
``WholeBody``                Bergman minimal-model loop embedding a beta-cell,
                             with an optional dynamic blood-formate pool.
``CoupledModel``             three-compartment (beta-cell + adipocyte + blood)
                             obesity/diabetes/supplementation/metformin model.

Time unit: hours. Concentrations: mM (insulin/lipid pools in arbitrary units).
Parameter provenance tags in the dataclasses: [PUB] Formate-NUDT5 core; [FIT]
fitted here; [LIT] literature; [ASM] assumed.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


# ==========================================================================
# Rate-law primitives (module-level, pure functions)
# ==========================================================================
def hill(x, K, h):
    """Hill activation x^h / (K^h + x^h)."""
    return x ** h / (K ** h + x ** h)


def mm(x, V, K):
    """Michaelis-Menten V x / (K + x)."""
    return V * x / (K + x)


# ==========================================================================
# Base class: shared numerical machinery
# ==========================================================================
class MetabolicCell:
    """Common base for the single-cell metabolic models.

    Subclasses define ``STATE`` (ordered state names), an initial condition
    ``Y0``, and the model-specific ``fluxes`` / ``derivatives``. This base
    supplies the state bookkeeping and the two solvers (time integration and
    log-space steady state) that every compartment reuses.
    """
    STATE: list = []
    Y0: np.ndarray = np.array([])

    def __init__(self):
        self.IDX = {s: i for i, s in enumerate(self.STATE)}
        self.N = len(self.STATE)

    @property
    def y0(self):
        return self.Y0.copy()

    # -- rate-law primitives available to every subclass --------------------
    @staticmethod
    def hill(x, K, h):
        return hill(x, K, h)

    @staticmethod
    def mm(x, V, K):
        return mm(x, V, K)

    @staticmethod
    def adenylate_kinase(adp, amp, atp, k_ak, K_ak):
        """Near-equilibrium 2 ADP <-> ATP + AMP relaxation flux."""
        return k_ak * (adp ** 2 - amp * atp / K_ak)

    # -- shared numerics ----------------------------------------------------
    def _integrate(self, rhs, y0, t_end, n_out=200, rtol=1e-8, atol=1e-10):
        """BDF integration of dy/dt = rhs(t, y)."""
        return solve_ivp(rhs, (0.0, t_end), y0, method='BDF', rtol=rtol,
                         atol=atol, dense_output=True,
                         t_eval=np.linspace(0.0, t_end, n_out))

    def _polish(self, res_log, y_guess):
        """Log-space least-squares steady-state polish; returns the state."""
        sol = least_squares(res_log, np.log(np.maximum(y_guess, 1e-9)),
                            xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=5000)
        return np.exp(sol.x)

    def fluxes(self, y, *args, **kw):        # pragma: no cover - interface
        raise NotImplementedError

    def derivatives(self, y, *args, **kw):   # pragma: no cover - interface
        raise NotImplementedError


# ==========================================================================
# Beta-cell parameters
# ==========================================================================
@dataclass
class BetaParams:
    """Parameters of the beta-cell model. Fields marked ``init=False`` are
    derived (from ``m``, ``v_pur``, ``for0/ffor0`` or the glue fit) in
    ``__post_init__``; call :meth:`recompute` after mutating a base field."""
    # --- energy / adenine nucleotides ---
    Kgk: float = 8.0          # mM glucokinase half-saturation [LIT]
    hgk: float = 1.7          # glucokinase Hill [LIT]
    m: float = 84.0           # mM/h maintenance ATP turnover [PUB]
    BASALFUEL: float = 0.40   # glucose-independent fuel floor [FIT]
    Eg: float = 0.25          # mM glycolysis ADP half-sat [PUB]
    Eo: float = 0.12          # mM OxPhos ADP half-sat [PUB]
    a: float = 168.0          # mM/h max ATPase rate [FIT]
    A: float = 2.0            # mM ATPase ATP half-sat [PUB]
    K_ak: float = 0.8         # adenylate-kinase equilibrium [PUB]
    k_ak: float = 5.0e3       # 1/h adenylate-kinase relaxation [ASM]
    # --- one-carbon / de novo purine ---
    Ftot: float = 0.02        # mM folate pool [ASM]
    H: float = 0.01           # mM FTHFS formate half-sat [PUB]
    Kc: float = 0.002         # mM GART/ATIC CHO-THF half-sat [ASM]
    Kg: float = 0.01          # mM GART GAR half-sat [ASM]
    Kt: float = 0.01          # mM ATIC AICAR half-sat [ASM]
    v_pur: float = 0.30       # mM/h resting purine-turnover scale [ASM]
    for0: float = 0.8         # mM reference intracellular formate [PUB]
    ffor0: float = 0.5        # mM/h reference formate exchange [PUB]
    FMITO: float = 0.5        # mM/h mitochondrial formate capacity [ASM]
    Kprpp: float = 0.3        # mM PPAT PRPP half-sat [ASM]
    ko_prpp: float = 2.0      # 1/h lumped PRPP drain [ASM]
    kdeg_int: float = 0.1     # 1/h GAR/AICAR drain [ASM]
    KdP: float = 0.104        # mM glue PRPP IC50 at zero AMP [FIT]
    # --- purine branch / degradation / urate ---
    k_amp_syn: float = 2.0    # ADSS IMP->S-AMP (per ATP/ADP) [ASM]
    Vmax_samp: float = 0.6    # mM/h ADSL capacity [FIT]
    Ksamp: float = 0.02       # mM ADSL half-sat [ASM]
    k_gmp_syn: float = 2.0    # IMP->GMP (per ATP/ADP) [ASM]
    k_imp_deg: float = 1.0    # 1/h IMP degradation [ASM]
    k_amp_deg: float = 0.5    # 1/h AMP degradation [ASM]
    k_ade_turn: float = 0.03  # 1/h adenylate catabolism [ASM]
    k_gnp_deg: float = 0.5    # 1/h GXP degradation [ASM]
    k_gnp_use: float = 0.3    # 1/h guanine-nucleotide usage [ASM]
    Khx: float = 0.02         # mM xanthine-oxidase half-sat [ASM]
    k_urate: float = 3.0      # 1/h urate secretion [ASM]
    # --- mTOR / insulin synthesis / secretion ---
    W_ATP_MTOR: float = 0.03  # ATP weight in mTOR purine signal [ASM]
    Kmp: float = 0.12         # mM mTOR purine half-sat [ASM]
    hmp: float = 2.0          # Hill purine sensing [ASM]
    Kec: float = 0.85         # energy-charge half-sat [ASM]
    hec: float = 4.0          # Hill energy-charge gate [ASM]
    aa: float = 2.0           # mM amino-acid pool [ASM]
    Kaa: float = 0.5          # mM translation half-sat [ASM]
    KG_tx: float = 7.0        # mM glucose insulin-gene half-max [FIT/LIT]
    hG_tx: float = 2.0        # Hill glucose transcription [LIT]
    B_TX: float = 0.10        # basal insulin-gene drive [FIT]
    B_MTOR: float = 0.29      # mTOR-independent insulin fraction [FIT]
    k_ins: float = 20.0       # a.u./h max pro-insulin synthesis [ASM]
    k_mat: float = 3.0        # 1/h maturation [LIT]
    cost_ins: float = 0.02    # ATP per a.u. insulin [ASM]
    cost_pur: float = 4.0     # ATP per mM purine [LIT]
    Ktrig: float = 12.0       # ATP/ADP secretion-trigger half-sat [ASM]
    htrig: float = 4.0        # Hill secretion trigger [ASM]
    B_TRIG: float = 0.05      # basal secretion [ASM]
    k_sec: float = 6.0        # 1/h max fractional secretion [LIT]
    mu: float = 1.0 / 240.0   # 1/h pool-renewal dilution [ASM]
    # --- derived (computed) ---
    e: float = field(init=False)
    eg: float = field(init=False)
    eo: float = field(init=False)
    h: float = field(init=False)
    Pmax: float = field(init=False)
    S_prpp: float = field(init=False)
    gmax: float = field(init=False)
    tmax: float = field(init=False)
    k_xo: float = field(init=False)
    kF: float = field(init=False)
    Ka: float = field(init=False)
    n_hill: float = field(init=False)

    def __post_init__(self):
        self.recompute()

    def recompute(self):
        """Recompute all derived parameters (energy scales, v_pur-derived
        capacities, formate exchange, glue fit) from the base fields."""
        self.e = 2.0 * self.m
        self.eg = 0.60 * self.e
        self.eo = 0.40 * self.e
        self.rescale_vpur()
        self.kF = self.ffor0 / self.for0
        self.Ka, self.n_hill = self._fit_glue()

    def rescale_vpur(self):
        """v_pur-derived capacities (used by the sensitivity analysis)."""
        self.h = 4.0 * self.v_pur
        self.Pmax = 2.0 * self.v_pur
        self.S_prpp = 2.0 * self.v_pur
        self.gmax = 5.0 * self.v_pur
        self.tmax = 5.0 * self.v_pur
        self.k_xo = 2.0 * self.v_pur

    def _fit_glue(self):
        """Fit the PPAT-NUDT5 glue AMP half-sat Ka and Hill n to the Witus
        et al. PRPP IC50-vs-AMP data (as switch_ppat.py)."""
        GLUE_AMP_MM = np.array([0.0, 0.1, 1.0])
        GLUE_IC50_UM = np.array([104.0, 327.0, 4900.0])

        def resid(p):
            Ka, n = np.exp(p)
            pred = self.KdP * 1000.0 * (1 + (GLUE_AMP_MM / Ka) ** n)
            return np.log(pred[1:]) - np.log(GLUE_IC50_UM[1:])

        f = least_squares(resid, np.log([0.05, 1.3]))
        return tuple(np.exp(f.x))


# ==========================================================================
# Beta-cell model
# ==========================================================================
class BetaCellModel(MetabolicCell):
    """Fifteen-state kinetic beta-cell: glucose-driven energy metabolism, de
    novo purine synthesis and degradation, urate secretion, and mTOR-controlled
    translation of insulin. See :class:`BetaParams` for parameters."""
    STATE = ['atp', 'adp', 'amp', 'form', 'cho', 'prpp', 'gar', 'aic',
             'imp', 'samp', 'gnp', 'hx', 'urate', 'proins', 'ins']
    Y0 = np.array([3.0, 0.3, 0.03, 0.3, 0.01, 0.1, 0.02, 0.02, 0.05, 0.02,
                   1.0, 0.02, 0.05, 1.0, 5.0])

    # calibration targets (independent-dataset anchors)
    GSIS_BASAL_G = 2.8
    GSIS_STIM_G = 16.7
    GSIS_EC50 = 8.0
    GSIS_SI = 5.5
    DET_G = (1.0, 10.0)
    DET_ATP_ADP = (2.4, 11.6)
    NI_CONTENT_RATIO = 0.50
    GOODING_G = (2.8, 16.7)
    GOODING_IMP_RATIO = 0.23
    GOODING_SAMP_RATIO = 3.4
    INS_MW = 5808.0
    INS_CONTENT_PGCELL = 20.0

    ASM_PARAMS = ['v_pur', 'BASALFUEL', 'Kmp', 'hmp', 'Kec', 'hec', 'k_ins',
                  'k_mat', 'B_MTOR', 'cost_ins', 'Ktrig', 'htrig', 'B_TRIG',
                  'k_sec', 'k_amp_syn', 'k_gmp_syn', 'k_imp_deg', 'Vmax_samp',
                  'k_amp_deg', 'k_ade_turn', 'Khx', 'k_urate', 'FMITO']

    def __init__(self, params: BetaParams = None):
        super().__init__()
        self.p = params if params is not None else BetaParams()

    # -- PPAT-NUDT5 glue ----------------------------------------------------
    def glue_theta(self, amp, prpp, fN=1.0):
        """Fraction of PPAT held in the inhibited PPAT-NUDT5 complex."""
        p = self.p
        x = (amp / p.Ka) ** p.n_hill
        return fN * (1.0 + x) / (1.0 + x + prpp / p.KdP)

    # -- rate laws ----------------------------------------------------------
    def fluxes(self, y, glucose, form_x, fN=1.0, mtor_mult=1.0, cox_ki=None):
        """All rates and derived quantities for state y at a given glucose.

        glucose   : extracellular glucose (mM)
        form_x    : blood formate (mM) -- the hypothesis lever
        fN        : PPAT-NUDT5 glue strength in [0,1]
        mtor_mult : mTORC1 activity multiplier (1 = WT, 0 = Raptor KO; also
                    models an mTOR inhibitor)
        cox_ki    : formate half-inhibition (mM) of complex IV (None = off)
        """
        p = self.p
        (atp, adp, amp, form, cho, prpp, gar, aic,
         imp, samp, gnp, hx, urate, proins, ins) = np.maximum(y, 1e-12)

        # energy: glucose sensing and ATP production/consumption
        gk = glucose ** p.hgk / (p.Kgk ** p.hgk + glucose ** p.hgk)
        fuel = p.BASALFUEL + (1.0 - p.BASALFUEL) * gk
        cox = 1.0 if cox_ki is None else 1.0 / (1.0 + (form / cox_ki) ** 2)
        v_glyc = p.eg * fuel * (adp / (p.Eg + adp))
        v_ox = p.eo * fuel * (adp / (p.Eo + adp)) * cox
        v_atp_prod = v_glyc + v_ox
        v_atpase = p.a * (atp / (p.A + atp))
        r_ak = p.k_ak * (adp ** 2 - amp * atp / p.K_ak)

        # purine synthesis (formate block)
        thf = max(p.Ftot - cho, 1e-12)
        v_fthfs = p.h * (form / (p.H + form)) * (thf / p.Ftot)
        onec = cho / (p.Kc + cho)
        theta = self.glue_theta(amp, prpp, fN)
        v_ppat = p.Pmax * (prpp / (p.Kprpp + prpp)) * (1.0 - theta)
        v_gart = p.gmax * (gar / (p.Kg + gar)) * onec
        v_atic = p.tmax * (aic / (p.Kt + aic)) * onec

        # purine branch and degradation
        adr = atp / adp
        v_imp_amp = p.k_amp_syn * imp * adr
        v_samp_amp = p.Vmax_samp * (samp / (p.Ksamp + samp))
        v_imp_gmp = p.k_gmp_syn * imp * adr
        v_imp_deg = p.k_imp_deg * imp
        v_amp_deg = p.k_amp_deg * amp
        v_ade_turn = p.k_ade_turn * atp
        v_gnp_deg = p.k_gnp_deg * gnp
        v_gnp_use = p.k_gnp_use * gnp
        v_xo = p.k_xo * (hx / (p.Khx + hx))
        v_urate_sec = p.k_urate * urate

        # mTOR and insulin synthesis
        Ppool = p.W_ATP_MTOR * atp + gnp
        ec = (atp + 0.5 * adp) / (atp + adp + amp)
        mtor = mtor_mult * (Ppool ** p.hmp / (p.Kmp ** p.hmp + Ppool ** p.hmp)) \
            * (ec ** p.hec / (p.Kec ** p.hec + ec ** p.hec))
        g_tx = p.B_TX + (1.0 - p.B_TX) * glucose ** p.hG_tx \
            / (p.KG_tx ** p.hG_tx + glucose ** p.hG_tx)
        v_ins_syn = p.k_ins * (p.B_MTOR + (1.0 - p.B_MTOR) * mtor) \
            * g_tx * (p.aa / (p.Kaa + p.aa))
        v_mat = p.k_mat * proins
        trig = p.B_TRIG + (1.0 - p.B_TRIG) * (atp / adp) ** p.htrig \
            / (p.Ktrig ** p.htrig + (atp / adp) ** p.htrig)
        v_sec = p.k_sec * ins * trig
        v_atp_cost = p.cost_ins * v_ins_syn + p.cost_pur * v_ppat

        return dict(
            gk=gk, v_glyc=v_glyc, v_ox=v_ox, v_atp_prod=v_atp_prod,
            v_atpase=v_atpase, r_ak=r_ak, v_fthfs=v_fthfs, onec=onec, theta=theta,
            v_ppat=v_ppat, v_gart=v_gart, v_atic=v_atic, v_imp_amp=v_imp_amp,
            v_samp_amp=v_samp_amp, v_imp_gmp=v_imp_gmp, v_imp_deg=v_imp_deg,
            v_amp_deg=v_amp_deg, v_ade_turn=v_ade_turn, v_gnp_deg=v_gnp_deg,
            v_gnp_use=v_gnp_use, v_xo=v_xo,
            v_urate_sec=v_urate_sec, Ppool=Ppool, ec=ec, mtor=mtor, g_tx=g_tx,
            v_ins_syn=v_ins_syn, v_mat=v_mat, trig=trig, v_sec=v_sec,
            v_atp_cost=v_atp_cost, atp_adp=atp / adp, cox=cox,
            atp=atp, adp=adp, amp=amp, form=form, cho=cho, prpp=prpp, gar=gar,
            aic=aic, imp=imp, samp=samp, gnp=gnp, hx=hx, urate=urate,
            proins=proins, ins=ins)

    def derivatives(self, y, glucose, form_x, fN=1.0, f_mito=None,
                    mtor_mult=1.0, cox_ki=None):
        """dX/dt for every species (the kinetic core)."""
        p = self.p
        i = self.IDX
        if f_mito is None:
            f_mito = p.FMITO
        r = self.fluxes(y, glucose, form_x, fN, mtor_mult, cox_ki)
        d = np.zeros(self.N)
        mu = p.mu

        d[i['atp']] = r['v_atp_prod'] - r['v_atpase'] - r['v_atp_cost'] \
            + r['r_ak'] - r['v_ade_turn'] - mu * r['atp']
        d[i['adp']] = -r['v_atp_prod'] + r['v_atpase'] + r['v_atp_cost'] \
            - 2.0 * r['r_ak'] - mu * r['adp']
        d[i['amp']] = r['r_ak'] + r['v_samp_amp'] - r['v_amp_deg'] - mu * r['amp']

        d[i['form']] = f_mito + p.kF * (form_x - r['form']) - r['v_fthfs']
        d[i['cho']] = r['v_fthfs'] - r['v_gart'] - r['v_atic']
        d[i['prpp']] = p.S_prpp - r['v_ppat'] - p.ko_prpp * r['prpp'] - mu * r['prpp']
        d[i['gar']] = r['v_ppat'] - r['v_gart'] - (mu + p.kdeg_int) * r['gar']
        d[i['aic']] = r['v_gart'] - r['v_atic'] - (mu + p.kdeg_int) * r['aic']

        d[i['imp']] = r['v_atic'] - r['v_imp_amp'] - r['v_imp_gmp'] \
            - r['v_imp_deg'] - mu * r['imp']
        d[i['samp']] = r['v_imp_amp'] - r['v_samp_amp'] - mu * r['samp']
        d[i['gnp']] = r['v_imp_gmp'] - r['v_gnp_deg'] - r['v_gnp_use'] - mu * r['gnp']

        d[i['hx']] = r['v_imp_deg'] + r['v_amp_deg'] + r['v_ade_turn'] \
            + r['v_gnp_deg'] - r['v_xo']
        d[i['urate']] = r['v_xo'] - r['v_urate_sec']

        d[i['proins']] = r['v_ins_syn'] - r['v_mat'] - mu * r['proins']
        d[i['ins']] = r['v_mat'] - r['v_sec'] - mu * r['ins']
        return d

    # -- solvers ------------------------------------------------------------
    def simulate(self, glucose, form_x, fN=1.0, t_end=24.0, y0=None, n_out=200,
                 glucose_fn=None, f_mito=None, mtor_mult=1.0, cox_ki=None):
        """Integrate the kinetic model in time; returns a trajectory DataFrame."""
        y0 = self.y0 if y0 is None else np.asarray(y0, float)
        G = glucose_fn if glucose_fn is not None else (lambda t: glucose)

        def rhs(t, yy):
            return self.derivatives(yy, G(t), form_x, fN, f_mito, mtor_mult, cox_ki)

        sol = self._integrate(rhs, y0, t_end, n_out=n_out)
        rows = []
        for k, t in enumerate(sol.t):
            yk = sol.y[:, k]
            r = self.fluxes(yk, G(t), form_x, fN, mtor_mult, cox_ki)
            rows.append({'t': t, 'glucose': G(t),
                         **{s: yk[j] for j, s in enumerate(self.STATE)},
                         'mtor': r['mtor'], 'ec': r['ec'], 'Ppool': r['Ppool'],
                         'v_ins_syn': r['v_ins_syn'], 'v_sec': r['v_sec'],
                         'trig': r['trig'], 'v_urate_sec': r['v_urate_sec'],
                         'v_ppat': r['v_ppat'], 'atp_adp': yk[0] / yk[1]})
        return pd.DataFrame(rows)

    def steady_state(self, glucose, form_x, fN=1.0, y0=None, f_mito=None,
                     t_pre=150.0, mtor_mult=1.0, cox_ki=None):
        """Solve derivatives(y)=0 by pre-integration then a log-space polish.
        Returns (fluxes dict with 'max_residual', state vector)."""
        traj = self.simulate(glucose, form_x, fN, t_end=t_pre, y0=y0, n_out=2,
                              f_mito=f_mito, mtor_mult=mtor_mult, cox_ki=cox_ki)
        y_guess = traj[self.STATE].iloc[-1].to_numpy()

        def res_log(z):
            yy = np.exp(np.clip(z, -30.0, 30.0))
            return self.derivatives(yy, glucose, form_x, fN, f_mito, mtor_mult,
                                    cox_ki) / (yy + 1e-6)

        y = self._polish(res_log, y_guess)
        r = self.fluxes(y, glucose, form_x, fN, mtor_mult, cox_ki)
        r['max_residual'] = float(np.max(np.abs(
            self.derivatives(y, glucose, form_x, fN, f_mito, mtor_mult, cox_ki))))
        return r, y

    # -- analysis scans -----------------------------------------------------
    def formate_scan(self, glucose, form_x_grid, fN=1.0, f_mito=None):
        rows, y0 = [], None
        for fx in form_x_grid:
            r, y = self.steady_state(glucose, fx, fN, y0=y0, f_mito=f_mito, t_pre=40.0)
            y0 = y
            rows.append({'form_x': fx, 'form_in': r['form'], 'atp': r['atp'],
                         'amp': r['amp'], 'gnp': r['gnp'], 'Ppool': r['Ppool'],
                         'ec': r['ec'], 'mtor': r['mtor'], 'v_ins_syn': r['v_ins_syn'],
                         'ins': r['ins'], 'v_sec': r['v_sec'], 'urate': r['urate'],
                         'v_urate_sec': r['v_urate_sec'], 'v_ppat': r['v_ppat'],
                         'max_residual': r['max_residual']})
        return pd.DataFrame(rows)

    def gsis_dose_response(self, form_x=0.3, G_grid=None, fN=1.0, f_mito=None):
        if G_grid is None:
            G_grid = np.array([2.8, 4, 5.5, 7, 8, 10, 12, 16.7, 20])
        rows, y0 = [], None
        for G in G_grid:
            r, y = self.steady_state(G, form_x, fN, y0=y0, f_mito=f_mito, t_pre=40.0)
            y0 = y
            rows.append({'glucose': G, 'v_sec': r['v_sec'], 'v_ins_syn': r['v_ins_syn'],
                         'mtor': r['mtor'], 'atp_adp': r['atp'] / r['adp']})
        return pd.DataFrame(rows)

    def _gsis_metrics(self, df):
        g, s = df.glucose.to_numpy(), df.v_sec.to_numpy()
        basal = np.interp(self.GSIS_BASAL_G, g, s)
        stim = np.interp(self.GSIS_STIM_G, g, s)
        half = basal + 0.5 * (s.max() - basal)
        ec50 = np.nan
        for i in range(1, len(g)):
            if (s[i - 1] - half) * (s[i] - half) <= 0 and s[i] != s[i - 1]:
                ec50 = g[i - 1] + (half - s[i - 1]) * (g[i] - g[i - 1]) / (s[i] - s[i - 1])
                break
        return ec50, stim / max(basal, 1e-9)

    def formate_deficiency_scan(self, glucose=8.0, fmito_grid=None, form_x=0.1, fN=1.0):
        if fmito_grid is None:
            fmito_grid = np.linspace(self.p.FMITO, 0.0, 21)
        rows, y0 = [], None
        for fm in fmito_grid:
            r, y = self.steady_state(glucose, form_x, fN, y0=y0, f_mito=fm, t_pre=40.0)
            y0 = y
            rows.append({'f_mito': fm, 'form_in': r['form'], 'Ppool': r['Ppool'],
                         'mtor': r['mtor'], 'v_ins_syn': r['v_ins_syn'],
                         'v_sec': r['v_sec'], 'ins': r['ins'], 'urate': r['urate'],
                         'v_urate_sec': r['v_urate_sec'], 'max_residual': r['max_residual']})
        return pd.DataFrame(rows)

    # -- calibration (fits parameters in place on self.p) -------------------
    def calibrate_energy(self, verbose=True):
        p = self.p

        def resid(q):
            p.BASALFUEL = 1.0 / (1.0 + np.exp(-q[0]))
            p.a = np.exp(q[1])
            out, y0 = [], None
            for G, tgt in zip(self.DET_G, self.DET_ATP_ADP):
                r, y0 = self.steady_state(G, 0.3, y0=y0, t_pre=60.0)
                out.append(np.log((r['atp'] / r['adp']) / tgt))
            return out

        q0 = [np.log(p.BASALFUEL / (1 - p.BASALFUEL)), np.log(p.a)]
        fit = least_squares(resid, q0)
        p.BASALFUEL = 1.0 / (1.0 + np.exp(-fit.x[0]))
        p.a = np.exp(fit.x[1])
        rlo, y = self.steady_state(self.DET_G[0], 0.3, t_pre=80.0)
        rhi, _ = self.steady_state(self.DET_G[1], 0.3, y0=y, t_pre=80.0)
        if verbose:
            print('energy calib [FIT]: BASALFUEL = %.3f, a = %.0f  ->  ATP/ADP '
                  '%.1f (t %.1f) @%gmM, %.1f (t %.1f) @%gmM; GXP x%.1f with glucose'
                  % (p.BASALFUEL, p.a, rlo['atp'] / rlo['adp'], self.DET_ATP_ADP[0],
                     self.DET_G[0], rhi['atp'] / rhi['adp'], self.DET_ATP_ADP[1],
                     self.DET_G[1], rhi['gnp'] / rlo['gnp']))
        return p.BASALFUEL, p.a

    def calibrate_imp_branch(self, verbose=True):
        p = self.p

        def resid(q):
            p.Vmax_samp = np.exp(q[0])
            rlo, y = self.steady_state(self.GOODING_G[0], 0.3, t_pre=60.0)
            rhi, _ = self.steady_state(self.GOODING_G[1], 0.3, y0=y, t_pre=60.0)
            samp_ratio = rhi['samp'] / max(rlo['samp'], 1e-12)
            return [np.log(samp_ratio / self.GOODING_SAMP_RATIO)]

        fit = least_squares(resid, [np.log(p.Vmax_samp)])
        p.Vmax_samp = np.exp(fit.x[0])
        rlo, y = self.steady_state(self.GOODING_G[0], 0.3, t_pre=80.0)
        rhi, _ = self.steady_state(self.GOODING_G[1], 0.3, y0=y, t_pre=80.0)
        if verbose:
            print('IMP-branch calib [FIT]: Vmax_samp = %.3f  ->  IMP x%.2f (t %.2f), '
                  'S-AMP x%.2f (t %.2f) over %g->%g mM'
                  % (p.Vmax_samp, rhi['imp'] / rlo['imp'], self.GOODING_IMP_RATIO,
                     rhi['samp'] / rlo['samp'], self.GOODING_SAMP_RATIO, *self.GOODING_G))
        return p.Vmax_samp

    def calibrate_mtor_content(self, glucose=8.0, verbose=True):
        p = self.p

        def resid(q):
            p.B_MTOR = 1.0 / (1.0 + np.exp(-q[0]))
            rwt, _ = self.steady_state(glucose, 0.3, t_pre=60.0)
            rko, _ = self.steady_state(glucose, 0.3, t_pre=60.0, mtor_mult=0.0)
            return [np.log((rko['ins'] / rwt['ins']) / self.NI_CONTENT_RATIO)]

        fit = least_squares(resid, [np.log(p.B_MTOR / (1 - p.B_MTOR))])
        p.B_MTOR = 1.0 / (1.0 + np.exp(-fit.x[0]))
        rwt, _ = self.steady_state(glucose, 0.3, t_pre=80.0)
        rko, _ = self.steady_state(glucose, 0.3, t_pre=80.0, mtor_mult=0.0)
        if verbose:
            print('mTOR-content calib [FIT]: B_MTOR = %.3f  ->  Raptor-KO content '
                  '%.0f%% of WT (target %.0f%%)'
                  % (p.B_MTOR, 100 * rko['ins'] / rwt['ins'], 100 * self.NI_CONTENT_RATIO))
        return p.B_MTOR

    def calibrate(self, verbose=True):
        p = self.p
        B_TX_MAX = 0.4

        def objective(q):
            p.KG_tx = np.exp(q[0])
            p.B_TX = B_TX_MAX / (1.0 + np.exp(-q[1]))
            df = self.gsis_dose_response()
            ec50, si = self._gsis_metrics(df)
            if not np.isfinite(ec50):
                ec50 = 100.0
            return [np.log(ec50 / self.GSIS_EC50), np.log(max(si, 1e-6) / self.GSIS_SI)]

        q0 = [np.log(p.KG_tx), 0.0]
        fit = least_squares(objective, q0,
                            bounds=([np.log(4.0), -np.inf], [np.log(12.0), np.inf]))
        p.KG_tx = np.exp(fit.x[0])
        p.B_TX = B_TX_MAX / (1.0 + np.exp(-fit.x[1]))
        df = self.gsis_dose_response()
        ec50, si = self._gsis_metrics(df)
        if verbose:
            print('calibration [FIT]: KG_tx = %.2f mM, B_TX = %.3f  ->  '
                  'EC50 = %.1f mM (target %.1f), SI = %.1f (target %.1f)'
                  % (p.KG_tx, p.B_TX, ec50, self.GSIS_EC50, si, self.GSIS_SI))
        return p.KG_tx, p.B_TX, ec50, si

    def calibrate_all(self, verbose=True):
        """Run the data-anchored calibrations in dependency order, then GSIS."""
        self.calibrate_energy(verbose)
        self.calibrate_imp_branch(verbose)
        self.calibrate_mtor_content(verbose=verbose)
        self.calibrate(verbose)
        return self

    # -- sensitivity and absolute units -------------------------------------
    def _formate_sensitivity(self, glucose=8.0):
        r_rep, _ = self.steady_state(glucose, 0.3, f_mito=self.p.FMITO, t_pre=60.0)
        r_def, _ = self.steady_state(glucose, 0.05, f_mito=0.2 * self.p.FMITO, t_pre=60.0)
        return r_rep['v_ins_syn'] / max(r_def['v_ins_syn'], 1e-9)

    def sensitivity(self, param_names=None, factor=2.0):
        if param_names is None:
            param_names = self.ASM_PARAMS
        p = self.p
        df0 = self.gsis_dose_response()
        _, si0 = self._gsis_metrics(df0)
        fs0 = self._formate_sensitivity()
        rows = [{'param': '(baseline)', 'dir': '-', 'value': np.nan,
                 'SI': si0, 'FS': fs0, 'converged': True}]
        for name in param_names:
            base = getattr(p, name)
            for lab, val in [('down', base / factor), ('up', base * factor)]:
                setattr(p, name, val)
                if name == 'v_pur':
                    p.rescale_vpur()
                try:
                    df = self.gsis_dose_response()
                    _, si = self._gsis_metrics(df)
                    fs = self._formate_sensitivity()
                    rows.append({'param': name, 'dir': lab, 'value': val,
                                 'SI': si, 'FS': fs, 'converged': True})
                except Exception:                       # pragma: no cover
                    rows.append({'param': name, 'dir': lab, 'value': val,
                                 'SI': np.nan, 'FS': np.nan, 'converged': False})
                setattr(p, name, base)
                if name == 'v_pur':
                    p.rescale_vpur()
        df = pd.DataFrame(rows)
        df['SI_rel'] = df.SI / si0
        df['FS_rel'] = df.FS / fs0
        return df

    def absolute_units(self, glucose_stim=16.7, form_x=0.3):
        r, _ = self.steady_state(glucose_stim, form_x)
        ins_unit = self.INS_CONTENT_PGCELL / r['ins']
        out = {'ins_unit_pg_per_au': ins_unit}
        for lab, G in [('basal 2.8mM', 2.8), ('stim 16.7mM', glucose_stim)]:
            rr, _ = self.steady_state(G, form_x)
            content = rr['ins'] * ins_unit
            sec = rr['v_sec'] * ins_unit
            out[lab] = dict(content_pg=content, sec_pg_per_h=sec,
                            sec_pct_content_per_h=100 * sec / max(content, 1e-9),
                            sec_fmol_per_h=sec / self.INS_MW * 1e3)
        return out


# ==========================================================================
# Adipocyte
# ==========================================================================
@dataclass
class AdipoParams:
    """Adipocyte parameters. Energy uses the beta-cell maintenance scale
    (m = 84); K_ak/k_ak/kF are the shared beta-cell constants."""
    m_ref: float = 84.0        # beta-cell maintenance ATP turnover (M.m)
    A_adip: float = 2.0
    Eg: float = 0.25
    Eo: float = 0.12
    KIglut: float = 0.5        # plasma-insulin half-max for GLUT4 [ASM]
    Kglut: float = 5.0         # mM glucose half-max for uptake [LIT]
    BGLUT: float = 0.10        # insulin-independent basal uptake [ASM]
    Ftot: float = 0.02
    H: float = 0.01
    Kc: float = 0.002
    kF: float = 0.5 / 0.8      # formate exchange (= M.kF)
    v_pur_a: float = 0.30
    Kprpp_a: float = 0.3
    prpp_a: float = 0.3
    K_ak: float = 0.8
    k_ak: float = 5.0e3
    k_amp_syn_a: float = 2.0
    k_gmp_syn_a: float = 2.0
    k_imp_deg_a: float = 1.0
    k_amp_deg_a: float = 0.5
    k_ade_turn_a: float = 0.03
    k_gnp_deg_a: float = 0.5
    Khx: float = 0.02
    k_hxexp: float = 6.0       # hypoxanthine export competing with XOR [ASM]
    k_urate_a: float = 3.0
    mu_a: float = 1.0 / 240.0
    KIlipo: float = 0.5
    k_lipo: float = 4.0
    k_lps: float = 0.2
    cost_lip: float = 0.03
    # derived
    e_adip: float = field(init=False)
    egA: float = field(init=False)
    eoA: float = field(init=False)
    aA: float = field(init=False)
    hA: float = field(init=False)
    VpurA: float = field(init=False)
    k_xo_a: float = field(init=False)

    def __post_init__(self):
        self.e_adip = 1.2 * self.m_ref
        self.egA = 0.6 * self.e_adip
        self.eoA = 0.4 * self.e_adip
        self.aA = self.e_adip
        self.hA = 4.0 * self.v_pur_a
        self.VpurA = 2.0 * self.v_pur_a
        self.k_xo_a = 3.0 * self.v_pur_a


class AdipocyteModel(MetabolicCell):
    """Ten-state adipocyte: shared energy/one-carbon/purine core, but with
    insulin-driven GLUT4 uptake and lipid storage, and an XOR-rich urate arm."""
    STATE = ['atp', 'adp', 'amp', 'form', 'cho', 'imp', 'gnp', 'hx', 'urate', 'lipid']
    Y0 = np.array([2.0, 0.2, 0.03, 0.3, 0.01, 0.05, 0.1, 0.02, 0.05, 5.0])

    def __init__(self, params: AdipoParams = None):
        super().__init__()
        self.p = params if params is not None else AdipoParams()

    def fluxes(self, y, glucose, insulin, form_x, mass=1.0, xor_mult=1.0):
        p = self.p
        atp, adp, amp, form, cho, imp, gnp, hx, urate, lipid = np.maximum(y, 1e-12)

        ins_glut = p.BGLUT + (1.0 - p.BGLUT) * insulin / (p.KIglut + insulin)
        uptake = ins_glut * (glucose / (p.Kglut + glucose))
        v_atp_prod = (p.egA + p.eoA) * uptake * (adp / (p.Eg + adp))
        v_atpase = p.aA * (atp / (p.A_adip + atp))
        r_ak = p.k_ak * (adp ** 2 - amp * atp / p.K_ak)

        thf = max(p.Ftot - cho, 1e-12)
        v_fthfs = p.hA * (form / (p.H + form)) * (thf / p.Ftot)
        onec = cho / (p.Kc + cho)
        v_denovo = p.VpurA * (p.prpp_a / (p.Kprpp_a + p.prpp_a)) * onec * (0.3 + 0.7 * uptake)

        adr = atp / adp
        v_imp_amp = p.k_amp_syn_a * imp * adr
        v_imp_gmp = p.k_gmp_syn_a * imp * adr
        v_imp_deg = p.k_imp_deg_a * imp
        v_amp_deg = p.k_amp_deg_a * amp
        v_ade_turn = p.k_ade_turn_a * atp
        v_gnp_deg = p.k_gnp_deg_a * gnp
        v_xo = xor_mult * p.k_xo_a * (hx / (p.Khx + hx))
        v_hx_exp = p.k_hxexp * hx
        v_urate_sec = p.k_urate_a * urate

        v_lipo = p.k_lipo * (insulin / (p.KIlipo + insulin)) * uptake
        v_lps = p.k_lps * lipid
        v_atp_cost = p.cost_lip * v_lipo

        return dict(uptake=uptake, v_atp_prod=v_atp_prod, v_atpase=v_atpase, r_ak=r_ak,
                    v_fthfs=v_fthfs, v_denovo=v_denovo, v_imp_amp=v_imp_amp,
                    v_imp_gmp=v_imp_gmp, v_imp_deg=v_imp_deg, v_amp_deg=v_amp_deg,
                    v_ade_turn=v_ade_turn, v_gnp_deg=v_gnp_deg, v_xo=v_xo,
                    v_hx_exp=v_hx_exp, v_urate_sec=v_urate_sec, v_lipo=v_lipo,
                    v_lps=v_lps, v_atp_cost=v_atp_cost, adr=adr,
                    atp=atp, adp=adp, amp=amp, form=form, cho=cho, imp=imp, gnp=gnp,
                    hx=hx, urate=urate, lipid=lipid)

    def derivatives(self, y, glucose, insulin, form_x, mass=1.0, xor_mult=1.0):
        p = self.p
        i = self.IDX
        r = self.fluxes(y, glucose, insulin, form_x, mass, xor_mult)
        d = np.zeros(self.N)
        mu = p.mu_a
        d[i['atp']] = r['v_atp_prod'] - r['v_atpase'] - r['v_atp_cost'] + r['r_ak'] \
            - r['v_ade_turn'] - mu * r['atp']
        d[i['adp']] = -r['v_atp_prod'] + r['v_atpase'] + r['v_atp_cost'] \
            - 2.0 * r['r_ak'] - mu * r['adp']
        d[i['amp']] = r['r_ak'] + r['v_imp_amp'] - r['v_amp_deg'] - mu * r['amp']
        d[i['form']] = p.kF * (form_x - r['form']) - r['v_fthfs']
        d[i['cho']] = r['v_fthfs'] - 2.0 * r['v_denovo']
        d[i['imp']] = r['v_denovo'] - r['v_imp_amp'] - r['v_imp_gmp'] \
            - r['v_imp_deg'] - mu * r['imp']
        d[i['gnp']] = r['v_imp_gmp'] - r['v_gnp_deg'] - mu * r['gnp']
        d[i['hx']] = r['v_imp_deg'] + r['v_amp_deg'] + r['v_ade_turn'] \
            + r['v_gnp_deg'] - r['v_xo'] - r['v_hx_exp']
        d[i['urate']] = r['v_xo'] - r['v_urate_sec']
        d[i['lipid']] = r['v_lipo'] - r['v_lps'] - mu * r['lipid']
        return d

    def steady_state(self, glucose, insulin, form_x, mass=1.0, xor_mult=1.0,
                     y0=None, t_pre=200.0):
        y0 = self.y0 if y0 is None else np.asarray(y0, float)
        sol = solve_ivp(lambda t, y: self.derivatives(y, glucose, insulin, form_x,
                                                       mass, xor_mult),
                        (0, t_pre), y0, method='BDF', rtol=1e-7, atol=1e-9,
                        t_eval=[t_pre])
        yg = sol.y[:, -1]

        def res(z):
            yy = np.exp(np.clip(z, -30, 30))
            return self.derivatives(yy, glucose, insulin, form_x, mass, xor_mult) / (yy + 1e-6)

        s = least_squares(res, np.log(np.maximum(yg, 1e-9)), xtol=1e-14, ftol=1e-14)
        y = np.exp(s.x)
        r = self.fluxes(y, glucose, insulin, form_x, mass, xor_mult)
        r['formate_uptake'] = self.p.kF * (form_x - r['form'])
        r['max_residual'] = float(np.max(np.abs(
            self.derivatives(y, glucose, insulin, form_x, mass, xor_mult))))
        return r, y


# ==========================================================================
# Eye compartment (composes a beta-cell)
# ==========================================================================
@dataclass
class EyeParams:
    EYE_FMITO_FRAC: float = 0.15   # fraction of beta-cell FMITO [ASM]
    K_regen: float = 3.0
    h_regen: float = 2.0
    EYE_COX_KI: float = 5.0        # mM formate/complex-IV half-inhibition [LIT,weak]
    Km_maint: float = 1.7          # retinal-maintenance insulin half-max [ASM]
    Kp_patho: float = 2.3          # pathology (VEGF/gliosis) half-max, higher [ASM]
    h_ret: float = 6.0             # steepness of both window arms [ASM]
    INDUCE_GAIN: float = 4.0       # ATF4 local-formate feed-forward gain [ASM]
    S6K_KI: float = 0.09           # S6K1->IRS-1 feedback half-max mTOR [FIT]
    H_S6K: float = 6.0             # feedback steepness [ASM]


class EyeModel:
    """Retinal insulin cell behind the blood-retinal barrier: the beta-cell
    machinery driven by blood glucose/formate with strictly local insulin, a low
    endogenous mitochondrial formate, and (optionally) a self-regulating local
    formate feed-forward. Composes a calibrated :class:`BetaCellModel`."""
    MW = 46.03

    def __init__(self, beta: BetaCellModel = None, params: EyeParams = None):
        self.beta = beta if beta is not None else BetaCellModel().calibrate_all(verbose=False)
        self.ep = params if params is not None else EyeParams()

    @property
    def EYE_FMITO(self):
        return self.ep.EYE_FMITO_FRAC * self.beta.p.FMITO

    # -- retinal outcome window --------------------------------------------
    def retinal_outcome(self, S):
        p = self.ep
        maint = S ** p.h_ret / (p.Km_maint ** p.h_ret + S ** p.h_ret)
        patho = S ** p.h_ret / (p.Kp_patho ** p.h_ret + S ** p.h_ret)
        return maint, patho, maint - patho

    def _s6k_feedback(self, mtor, s6k_ki=None, h=None):
        s6k_ki = self.ep.S6K_KI if s6k_ki is None else s6k_ki
        h = self.ep.H_S6K if h is None else h
        return s6k_ki ** h / (s6k_ki ** h + mtor ** h)

    # -- blood-dependent eye ------------------------------------------------
    def eye_state(self, glucose, blood_formate, y0=None, cox_ki=None):
        p = self.ep
        r, y = self.beta.steady_state(glucose, blood_formate, f_mito=self.EYE_FMITO,
                                      y0=y0, cox_ki=cox_ki)
        S = r['v_ins_syn']
        regen = S ** p.h_regen / (p.K_regen ** p.h_regen + S ** p.h_regen)
        return dict(eye_ins_prod=S, eye_mtor=r['mtor'], eye_formate=r['form'],
                    atp_adp=r['atp'] / r['adp'], ec=r['ec'], cox=r['cox'],
                    regen=regen), y

    # -- regulated local-formate feed-forward -------------------------------
    def _rhs_ff(self, t, y, glucose, Fb, gain, s6k_ki, mthfd2, mtor_mult):
        r = self.beta.fluxes(y, glucose, Fb, mtor_mult=mtor_mult)
        fb = 1.0 if s6k_ki is None else self._s6k_feedback(r['mtor'], s6k_ki)
        f_mito = mthfd2 * (self.EYE_FMITO + gain * r['mtor'] * fb)
        return self.beta.derivatives(y, glucose, Fb, f_mito=f_mito, mtor_mult=mtor_mult)

    def eye_state_ff(self, glucose, blood_formate, gain=None, s6k_ki='default',
                     y0=None, t_end=200.0, mthfd2=1.0, mtor_mult=1.0):
        """Regulated autonomous-eye steady state (integration to the stable
        attractor). gain=0 recovers the blood-dependent eye; s6k_ki=None the
        unregulated runaway. mthfd2<1 = intravitreal MTHFD2 inhibitor (local
        formate synthesis); mtor_mult<1 = mTOR (insulin-production) inhibitor."""
        gain = self.ep.INDUCE_GAIN if gain is None else gain
        if s6k_ki == 'default':
            s6k_ki = self.ep.S6K_KI
        y = self.beta.y0 if y0 is None else np.asarray(y0, float)

        def _reldot(yy):
            rr = self.beta.fluxes(yy, glucose, blood_formate, mtor_mult=mtor_mult)
            fbb = 1.0 if s6k_ki is None else self._s6k_feedback(rr['mtor'], s6k_ki)
            fmito = mthfd2 * (self.EYE_FMITO + gain * rr['mtor'] * fbb)
            return self.beta.derivatives(yy, glucose, blood_formate, f_mito=fmito,
                                         mtor_mult=mtor_mult) / (np.abs(yy) + 1e-6)

        resid = np.inf
        for _ in range(12):
            sol = solve_ivp(self._rhs_ff, (0, t_end), y, method='BDF', rtol=1e-8,
                            atol=1e-10, t_eval=[t_end],
                            args=(glucose, blood_formate, gain, s6k_ki, mthfd2, mtor_mult))
            y = sol.y[:, -1]
            resid = float(np.max(np.abs(_reldot(y))))
            if resid < 1e-7:
                break
        r = self.beta.fluxes(y, glucose, blood_formate, mtor_mult=mtor_mult)
        fb = 1.0 if s6k_ki is None else self._s6k_feedback(r['mtor'], s6k_ki)
        return dict(eye_ins_prod=r['v_ins_syn'], eye_mtor=r['mtor'],
                    eye_formate=r['form'],
                    f_mito_local=mthfd2 * (self.EYE_FMITO + gain * r['mtor'] * fb),
                    s6k_fb=fb, max_residual=resid), y

    # -- scans --------------------------------------------------------------
    def insulin_excess_window(self, glucose=5.0, Fb_grid=None):
        if Fb_grid is None:
            Fb_grid = np.concatenate([np.linspace(0.02, 1.0, 16),
                                      np.linspace(1.2, 12.0, 24)])
        rows, y = [], None
        for Fb in Fb_grid:
            s, y = self.eye_state(glucose, Fb, y0=y, cox_ki=None)
            m, p, h = self.retinal_outcome(s['eye_ins_prod'])
            rows.append({'blood_formate_mM': Fb, 'blood_formate_mgL': Fb * self.MW,
                         'eye_ins_prod': s['eye_ins_prod'], 'maintenance': m,
                         'pathology': p, 'retinal_health': h})
        return pd.DataFrame(rows)

    def formate_response(self, glucose=5.0, Fb_grid=None):
        if Fb_grid is None:
            Fb_grid = np.linspace(0.005, 0.20, 20)
        rows, y = [], None
        for Fb in Fb_grid:
            s, y = self.eye_state(glucose, Fb, y0=y)
            rows.append({'blood_formate_mM': Fb, 'blood_formate_mgL': Fb * self.MW,
                         'eye_ins_prod': s['eye_ins_prod'], 'eye_mtor': s['eye_mtor'],
                         'regen': s['regen']})
        return pd.DataFrame(rows)


# ==========================================================================
# Whole-body loop (Bergman minimal model) + blood-formate dynamics
# ==========================================================================
@dataclass
class WholeBodyParams:
    Gb0: float = 5.0     # mM basal glucose [LIT]
    Sg: float = 1.5      # 1/h glucose effectiveness [LIT]
    p2: float = 1.2      # 1/h remote-insulin turnover [LIT]
    Si: float = 1.0      # insulin sensitivity [ASM]
    kI: float = 30.0     # 1/h plasma insulin clearance [LIT]
    beta: float = 10.0   # islet mass/distribution [ASM]


@dataclass
class BloodFormateParams:
    P_endog: float = 0.02   # mM/h endogenous formate production [ASM]
    kren: float = 0.5       # 1/h renal formate clearance [ASM]
    U0: float = 0.01        # mM/h insulin-independent formate sink [ASM]
    Uins: float = 0.12      # mM/h max insulin-driven formate sink [ASM]
    KF: float = 0.05        # mM formate half-sat for consumption [ASM]


DIETS = {'low': 0.010, 'normal': 0.030, 'high': 0.080}   # mM/h dietary formate/serine


class WholeBody:
    """Bergman minimal-model glucose-insulin loop embedding a beta-cell as the
    insulin source. Discriminates central (beta-cell) from peripheral (glucose-
    effectiveness) formate action."""

    def __init__(self, beta: BetaCellModel = None, params: WholeBodyParams = None):
        self.beta = beta if beta is not None else BetaCellModel().calibrate_all(verbose=False)
        self.wp = params if params is not None else WholeBodyParams()
        self.NBC = self.beta.N

    def _basal(self, form_x, f_mito, peri, fN=1.0):
        wp = self.wp
        y_bc = self.beta.y0
        tr = self.beta.simulate(wp.Gb0, form_x, fN=fN, t_end=200.0, y0=y_bc,
                                n_out=2, f_mito=f_mito)
        y_bc = tr[self.beta.STATE].iloc[-1].to_numpy()
        r = self.beta.fluxes(y_bc, wp.Gb0, form_x, fN)
        Ip = wp.beta * r['v_sec'] / wp.kI
        y = np.concatenate([y_bc, [wp.Gb0, Ip, 0.0]])
        return y, Ip

    def _rhs(self, t, y, form_x, f_mito, peri, Ib_ref, fN, Ra):
        wp = self.wp
        n = self.NBC
        y_bc, Gb, Ip, X = y[:n], y[n], y[n + 1], y[n + 2]
        d_bc = self.beta.derivatives(y_bc, max(Gb, 1e-6), form_x, fN, f_mito)
        r = self.beta.fluxes(y_bc, max(Gb, 1e-6), form_x, fN)
        Sg_eff = wp.Sg * (1.0 + peri)
        dGb = Ra(t) - (Sg_eff + X) * Gb + Sg_eff * wp.Gb0
        dIp = wp.beta * r['v_sec'] - wp.kI * Ip
        dX = -wp.p2 * X + wp.p2 * wp.Si * (Ip - Ib_ref)
        return np.concatenate([d_bc, [dGb, dIp, dX]])

    def gtt(self, form_x=0.1, f_mito=None, peri=0.0, fN=1.0, bolus=8.0, t_end=4.0,
            n_out=240, Ib_ref=None):
        wp = self.wp
        n = self.NBC
        if f_mito is None:
            f_mito = self.beta.p.FMITO
        y0, Ib = self._basal(form_x, f_mito, peri, fN)
        if Ib_ref is None:
            Ib_ref = Ib
        y0 = y0.copy()
        y0[n] += bolus

        def Ra(t):
            return 0.0

        sol = solve_ivp(self._rhs, (0.0, t_end), y0, method='BDF', rtol=1e-7, atol=1e-9,
                        args=(form_x, f_mito, peri, Ib_ref, fN, Ra),
                        t_eval=np.linspace(0, t_end, n_out))
        rows = []
        for k, t in enumerate(sol.t):
            y = sol.y[:, k]
            r = self.beta.fluxes(y[:n], max(y[n], 1e-6), form_x, fN)
            rows.append({'t': t, 'Gb': y[n], 'Ip': y[n + 1], 'X': y[n + 2],
                         'v_sec': r['v_sec'], 'mtor': r['mtor'],
                         'ins': y[self.beta.IDX['ins']]})
        return pd.DataFrame(rows), Ib

    @staticmethod
    def _auc_above(t, y, base):
        return float(np.trapz(np.maximum(np.asarray(y) - base, 0.0), t))

    def scenarios(self, peri_val=0.30):
        wp = self.wp
        defic = 0.2 * self.beta.p.FMITO
        repl = self.beta.p.FMITO
        _, Ib_ref = self._basal(0.1, defic, 0.0)
        conds = [
            ('baseline (formate-deficient)', dict(f_mito=defic, peri=0.0)),
            ('central (beta-cell formate)', dict(f_mito=repl, peri=0.0)),
            ('peripheral (Carpene, Sg up)', dict(f_mito=defic, peri=peri_val)),
            ('both', dict(f_mito=repl, peri=peri_val)),
        ]
        rows, trajs = [], {}
        for name, kw in conds:
            tr, Ib = self.gtt(form_x=0.1, Ib_ref=Ib_ref, **kw)
            trajs[name] = tr
            g_auc = self._auc_above(tr.t, tr.Gb, wp.Gb0)
            i_auc = float(np.trapz(tr.Ip, tr.t))
            rows.append({'scenario': name, 'fasting_Gb': tr.Gb.iloc[-1],
                         'peak_Gb': tr.Gb.max(), 'Gb_2h': np.interp(2.0, tr.t, tr.Gb),
                         'glucose_AUC': g_auc, 'insulin_AUC': i_auc,
                         'peak_Ip': tr.Ip.max()})
        df = pd.DataFrame(rows)
        base = df.iloc[0]
        df['dGlucoseAUC_%'] = 100 * (df.glucose_AUC / base.glucose_AUC - 1)
        df['dInsulinAUC_%'] = 100 * (df.insulin_AUC / base.insulin_AUC - 1)
        return df, trajs


class FormateDiabetes:
    """Whole-body loop with a dynamic blood-formate pool: an insulin-secretion
    defect reduces the anabolic formate sink, so blood formate overflows upward
    (the Takase-vs-Pietzke reconciliation)."""
    MW = 46.03
    DIETS = DIETS

    def __init__(self, beta: BetaCellModel = None, wp: WholeBodyParams = None,
                 fp: BloodFormateParams = None):
        self.beta = beta if beta is not None else BetaCellModel().calibrate_all(verbose=False)
        self.wp = wp if wp is not None else WholeBodyParams()
        self.fp = fp if fp is not None else BloodFormateParams()
        self.NBC = self.beta.N

    def _rhs(self, t, y, sec_defect, D_diet, Ib_ref, KI):
        wp, fp, n = self.wp, self.fp, self.NBC
        ybc = y[:n]
        Gb, Ip, X, Fb = y[n], y[n + 1], y[n + 2], y[n + 3]
        Gb = max(Gb, 1e-6)
        Fb = max(Fb, 1e-9)
        d_bc = self.beta.derivatives(ybc, Gb, Fb, fN=1.0)
        r = self.beta.fluxes(ybc, Gb, Fb, 1.0)
        v_sec = r['v_sec']
        dGb = -(wp.Sg + X) * Gb + wp.Sg * wp.Gb0
        dIp = wp.beta * sec_defect * v_sec - wp.kI * Ip
        dX = -wp.p2 * X + wp.p2 * wp.Si * (Ip - Ib_ref)
        ins_frac = Ip / (KI + Ip)
        U_form = (fp.U0 + fp.Uins * ins_frac) * (Fb / (fp.KF + Fb))
        dFb = D_diet + fp.P_endog - U_form - fp.kren * Fb
        return np.concatenate([d_bc, [dGb, dIp, dX, dFb]])

    def _settle(self, sec_defect, D_diet, Ib_ref, KI, y0, t_end=400.0):
        sol = solve_ivp(self._rhs, (0.0, t_end), y0, method='BDF', rtol=1e-7, atol=1e-9,
                        args=(sec_defect, D_diet, Ib_ref, KI), t_eval=[t_end])
        return sol.y[:, -1]

    def _healthy_reference(self, D_diet_normal=None):
        n = self.NBC
        D_diet_normal = self.DIETS['normal'] if D_diet_normal is None else D_diet_normal
        y0 = np.concatenate([self.beta.y0, [self.wp.Gb0, 0.5, 0.0, 0.05]])
        KI = 0.5
        y = self._settle(1.0, D_diet_normal, 0.05, KI, y0)
        Ib = y[n + 1]
        KI = max(Ib, 1e-3)
        y = self._settle(1.0, D_diet_normal, Ib, KI, y)
        Ib = y[n + 1]
        return Ib, KI, y

    def perturbation(self, sec_grid=None):
        n = self.NBC
        if sec_grid is None:
            sec_grid = np.linspace(1.0, 0.15, 16)
        Ib_ref, KI, y_healthy = self._healthy_reference()
        rows = []
        for diet, D in self.DIETS.items():
            y0 = y_healthy.copy()
            for sd in sec_grid:
                y0 = self._settle(sd, D, Ib_ref, KI, y0)
                Gb, Ip, Fb = y0[n], y0[n + 1], y0[n + 3]
                rows.append({'diet': diet, 'D_diet': D, 'sec_defect': sd,
                             'glucose_mM': Gb, 'insulin': Ip, 'formate_mM': Fb,
                             'formate_mgL': Fb * self.MW})
        return pd.DataFrame(rows), Ib_ref


# ==========================================================================
# Three-compartment coupled model (beta-cell + adipocyte + blood)
# ==========================================================================
class CoupledModel:
    """Pancreas + adipose + blood pools. Obesity scales adipose mass, XOR and
    insulin resistance; diabetes is a secretion defect; supports supplementation
    and the dual-action metformin analyses."""
    MW = 46.03
    URATE_MGDL = 16.8
    c_adip = 0.16          # tissue->blood coupling [ASM]
    kren_u = 0.3           # 1/h renal urate clearance [ASM]
    U0_urate = 0.06        # mM/h non-adipose urate production [ASM]
    OBESE_MASS = 3.0
    OBESE_XOR = 2.0
    OBESE_SI = 0.5
    DIAB_DEFECT = 0.40
    METF_PERIPH = 1.0
    METF_FORM = 0.6
    DIETS = DIETS

    def __init__(self, beta: BetaCellModel = None, adipo: AdipocyteModel = None,
                 wp: WholeBodyParams = None, fp: BloodFormateParams = None):
        self.beta = beta if beta is not None else BetaCellModel().calibrate_all(verbose=False)
        self.adipo = adipo if adipo is not None else AdipocyteModel()
        self.wp = wp if wp is not None else WholeBodyParams()
        self.fp = fp if fp is not None else BloodFormateParams()
        self.NBC = self.beta.N
        self.NAD = self.adipo.N

    def _rhs(self, t, y, sec_defect, D_diet, Ib_ref, KI, M_adip, xor_mult, Si,
             f_mito_bc=None, sg_mult=1.0, pendog=None):
        wp, fp = self.wp, self.fp
        nbc, nad = self.NBC, self.NAD
        bc = y[:nbc]
        ad = y[nbc:nbc + nad]
        Gb, Ip, X, Fb, Ub = y[nbc + nad:nbc + nad + 5]
        Gb = max(Gb, 1e-6)
        Fb = max(Fb, 1e-9)
        Ip = max(Ip, 0.0)
        if pendog is None:
            pendog = fp.P_endog

        d_bc = self.beta.derivatives(bc, Gb, Fb, fN=1.0, f_mito=f_mito_bc)
        r_bc = self.beta.fluxes(bc, Gb, Fb, 1.0)
        v_sec = r_bc['v_sec']

        d_ad = self.adipo.derivatives(ad, Gb, Ip, Fb, mass=M_adip, xor_mult=xor_mult)
        r_ad = self.adipo.fluxes(ad, Gb, Ip, Fb, mass=M_adip, xor_mult=xor_mult)
        ad_form_uptake = self.adipo.p.kF * (Fb - r_ad['form'])
        ad_urate_sec = r_ad['v_urate_sec']

        Sg_eff = wp.Sg * sg_mult
        dGb = -(Sg_eff + X) * Gb + Sg_eff * wp.Gb0
        dIp = wp.beta * sec_defect * v_sec - wp.kI * Ip
        dX = -wp.p2 * X + wp.p2 * Si * (Ip - Ib_ref)
        U_syst = (fp.U0 + fp.Uins * Ip / (KI + Ip)) * (Fb / (fp.KF + Fb))
        dFb = D_diet + pendog - U_syst - self.c_adip * M_adip * ad_form_uptake - fp.kren * Fb
        dUb = self.c_adip * M_adip * ad_urate_sec + self.U0_urate - self.kren_u * Ub
        return np.concatenate([d_bc, d_ad, [dGb, dIp, dX, dFb, dUb]])

    def _settle(self, sec_defect, D_diet, Ib_ref, KI, M_adip, xor_mult, Si, y0,
                t_end=600.0, f_mito_bc=None, sg_mult=1.0, pendog=None):
        sol = solve_ivp(self._rhs, (0, t_end), y0, method='BDF', rtol=1e-6, atol=1e-9,
                        args=(sec_defect, D_diet, Ib_ref, KI, M_adip, xor_mult, Si,
                              f_mito_bc, sg_mult, pendog), t_eval=[t_end])
        return sol.y[:, -1]

    def _reference(self, D_diet):
        nbc, nad = self.NBC, self.NAD
        y0 = np.concatenate([self.beta.y0, self.adipo.y0,
                             [self.wp.Gb0, 0.5, 0.0, 0.05, 0.3]])
        KI = 0.5
        y = self._settle(1.0, D_diet, 0.05, KI, 1.0, 1.0, self.wp.Si, y0)
        Ib = y[nbc + nad + 1]
        KI = max(Ib, 1e-3)
        y = self._settle(1.0, D_diet, Ib, KI, 1.0, 1.0, self.wp.Si, y)
        Ib = y[nbc + nad + 1]
        return Ib, KI, y

    def _groups(self):
        S = self.wp.Si
        return [
            ('lean, non-diabetic', 1.0, 1.0, 1.0, S),
            ('obese, non-diabetic', 1.0, self.OBESE_MASS, self.OBESE_XOR, self.OBESE_SI * S),
            ('lean, diabetic', self.DIAB_DEFECT, 1.0, 1.0, S),
            ('obese, diabetic', self.DIAB_DEFECT, self.OBESE_MASS, self.OBESE_XOR, self.OBESE_SI * S),
        ]

    def scenarios(self, D_diet=None):
        nbc, nad = self.NBC, self.NAD
        D_diet = self.DIETS['normal'] if D_diet is None else D_diet
        Ib, KI, y0 = self._reference(D_diet)
        rows, y = [], y0.copy()
        for name, sd, Ma, xm, Si in self._groups():
            y = self._settle(sd, D_diet, Ib, KI, Ma, xm, Si, y)
            Gb, Ip, X, Fb, Ub = y[nbc + nad:nbc + nad + 5]
            rows.append({'scenario': name, 'glucose_mM': Gb, 'insulin': Ip,
                         'formate_mM': Fb, 'formate_mgL': Fb * self.MW,
                         'urate_mM': Ub, 'urate_mgdL': Ub * self.URATE_MGDL})
        df = pd.DataFrame(rows)
        base = df.iloc[0]
        df['formate_vs_healthy_%'] = 100 * (df.formate_mM / base.formate_mM - 1)
        df['urate_vs_healthy_%'] = 100 * (df.urate_mM / base.urate_mM - 1)
        return df

    def intervention(self, D_grid=None):
        nbc, nad = self.NBC, self.NAD
        if D_grid is None:
            D_grid = np.linspace(0.01, 0.16, 9)
        fdef = 0.2 * self.beta.p.FMITO
        types = {
            'lean, non-diabetic': (1.0, 1.0, 1.0, None),
            'lean diabetic, formate-replete': (self.DIAB_DEFECT, 1.0, 1.0, None),
            'lean diabetic, formate-deficient': (self.DIAB_DEFECT, 1.0, 1.0, fdef),
            'obese diabetic, formate-deficient': (self.DIAB_DEFECT, self.OBESE_MASS,
                                                  self.OBESE_XOR, fdef),
        }
        Ib, KI, y0 = self._reference(self.DIETS['normal'])
        rows = []
        for name, (sd, Ma, xm, fm) in types.items():
            Si = (self.OBESE_SI * self.wp.Si) if 'obese' in name else self.wp.Si
            y = y0.copy()
            for D in D_grid:
                y = self._settle(sd, D, Ib, KI, Ma, xm, Si, y, f_mito_bc=fm)
                Gb, Ip, X, Fb, Ub = y[nbc + nad:nbc + nad + 5]
                bc = y[:nbc]
                rows.append({'group': name, 'D_diet': D, 'glucose_mM': Gb,
                             'insulin': Ip, 'formate_mM': Fb, 'formate_mgL': Fb * self.MW,
                             'bc_formate': bc[self.beta.IDX['form']],
                             'urate_mgdL': Ub * self.URATE_MGDL})
        return pd.DataFrame(rows)

    def metformin(self, dose_grid=None):
        nbc, nad = self.NBC, self.NAD
        if dose_grid is None:
            dose_grid = np.linspace(0.0, 1.0, 9)
        pts = {'lean diabetic': (1.0, 1.0, self.wp.Si),
               'obese diabetic': (self.OBESE_MASS, self.OBESE_XOR, self.OBESE_SI * self.wp.Si)}
        Ib, KI, y0 = self._reference(self.DIETS['normal'])
        D = self.DIETS['normal']
        rows = []
        for name, (Ma, xm, Si) in pts.items():
            y = y0.copy()
            for d in dose_grid:
                sg = 1.0 + self.METF_PERIPH * d
                pend = self.fp.P_endog * (1.0 - self.METF_FORM * d)
                fbc = self.beta.p.FMITO * (1.0 - self.METF_FORM * d)
                y = self._settle(self.DIAB_DEFECT, D, Ib, KI, Ma, xm, Si, y,
                                 f_mito_bc=fbc, sg_mult=sg, pendog=pend)
                Gb, Ip, X, Fb, Ub = y[nbc + nad:nbc + nad + 5]
                rows.append({'patient': name, 'metformin_dose': d, 'glucose_mM': Gb,
                             'insulin': Ip, 'formate_mM': Fb, 'formate_mgL': Fb * self.MW,
                             'urate_mgdL': Ub * self.URATE_MGDL})
        return pd.DataFrame(rows)


# ==========================================================================
# Gut-microbiome urate cycle (Vazquez 2020, "The Urate Cycle and Obesity")
# ==========================================================================
class GutUrateCycle(CoupledModel):
    """CoupledModel + a gut-microbiome compartment. The gut microbiome degrades
    blood urate and returns formate to the circulation, closing a positive
    feedback:  formate -> de novo purine synthesis -> urate -> (gut microbiome)
    -> formate.  Per the hypothesis one of the two purine-derived formate carbons
    of each urate is recovered as formate (the glycine carbons return as acetate,
    feeding lipogenesis -- noted, not tracked).  The gut processing/transit delay
    is modelled by a linear chain of n_gut compartments (mean delay
    tau_gut = n_gut / k_g), the linear-chain trick for a gamma-distributed lag.

    k_gut  : 1/h, rate at which blood urate is taken up by the gut route [ASM]
    Y_form : mol formate returned per mol urate degraded (1, per hypothesis) [ASM]
    tau_gut: h, mean gut processing/transit delay [ASM]
    """
    def __init__(self, beta=None, adipo=None, wp=None, fp=None,
                 k_gut=0.2, tau_gut=12.0, Y_form=1.0, n_gut=3):
        super().__init__(beta, adipo, wp, fp)
        self.k_gut = k_gut
        self.tau_gut = tau_gut
        self.Y_form = Y_form
        self.n_gut = n_gut
        self.k_g = n_gut / tau_gut

    def _rhs(self, t, y, *args, **kw):
        nc = self.NBC + self.NAD + 5
        d = np.asarray(super()._rhs(t, y[:nc], *args, **kw), float)
        g = y[nc:nc + self.n_gut]
        Ub = max(y[nc - 1], 1e-9)
        gut_uptake = self.k_gut * Ub                 # blood urate -> gut microbiome
        d[nc - 1] -= gut_uptake                       # urate diverted from renal
        dg = np.empty(self.n_gut)
        dg[0] = gut_uptake - self.k_g * g[0]
        for i in range(1, self.n_gut):
            dg[i] = self.k_g * g[i - 1] - self.k_g * g[i]
        d[nc - 2] += self.Y_form * self.k_g * g[-1]   # delayed formate return
        return np.concatenate([d, dg])

    def _reference(self, D_diet):
        nc = self.NBC + self.NAD + 5
        y0 = np.concatenate([self.beta.y0, self.adipo.y0,
                             [self.wp.Gb0, 0.5, 0.0, 0.05, 0.3], np.zeros(self.n_gut)])
        KI = 0.5
        y = self._settle(1.0, D_diet, 0.05, KI, 1.0, 1.0, self.wp.Si, y0)
        Ib = y[nc - 4]
        KI = max(Ib, 1e-3)
        y = self._settle(1.0, D_diet, Ib, KI, 1.0, 1.0, self.wp.Si, y)
        return y[nc - 4], KI, y
