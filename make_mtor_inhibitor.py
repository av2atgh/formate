#!/usr/bin/env python3
"""Insulin-production (mTORC1) inhibitor in the eye under methanol (fig21_mtor_inhibitor.pdf).

Insulin synthesis in the model is mTOR-driven translation, so an mTORC1 inhibitor
(rapamycin/sirolimus; reduces proinsulin biosynthesis) is modelled by mtor_mult<1
in eye.eye_state_ff. Unlike the MTHFD2 inhibitor (which acts UPSTREAM on local
formate and cannot touch a blood-borne methanol flood), the mTOR inhibitor acts
DOWNSTREAM at the synthesis node and caps insulin regardless of formate source ---
so it rescues methanol optic toxicity. Three panels: (A) net retinal health vs
blood formate for no drug, MTHFD2 inhibitor (upstream, fails), mTOR inhibitor
(downstream, rescues); (B) dose-response at lethal methanol (therapeutic window);
(C) eye insulin -- the mTOR inhibitor clips the excess while sparing the baseline.
"""
import numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
os.chdir('/Users/avazquez/av2atg/formate/insulin')
from model import BetaCellModel, EyeModel
bc = BetaCellModel(); bc.calibrate_all(verbose=False)
E = EyeModel(beta=bc)
from figstyle import setup, panel
setup()
S=np.linspace(1.4,5,3000); Sopt=S[np.array([E.retinal_outcome(s)[2] for s in S]).argmax()]
H=lambda s: E.retinal_outcome(s)[2]

fig,ax=plt.subplots(1,3,figsize=(12.6,4.0),constrained_layout=True)
Fb=np.concatenate([np.linspace(0.3,5,10),np.geomspace(6,600,16)])
def marks(a):
    for v,c in [(5,'#2ca02c'),(50,'#ff7f0e'),(500,'#d62728')]: a.axvline(v,ls=':',color=c,lw=1)

# (A) net health vs blood formate: no drug / MTHFD2i (upstream) / mTORi (downstream)
h_no,h_mthfd2,h_mtor,ins_no,ins_mtor=[],[],[],[],[]
for fbmg in Fb:
    f=fbmg/46.03
    s0,_=E.eye_state_ff(5.0,f,y0=None);                 h_no.append(H(s0['eye_ins_prod'])); ins_no.append(s0['eye_ins_prod'])
    sm,_=E.eye_state_ff(5.0,f,mthfd2=0.0,y0=None);       h_mthfd2.append(H(sm['eye_ins_prod']))
    st,_=E.eye_state_ff(5.0,f,mtor_mult=0.5,y0=None);    h_mtor.append(H(st['eye_ins_prod'])); ins_mtor.append(st['eye_ins_prod'])
ax[0].plot(Fb,h_no,'-',color='#7f7f7f',lw=2,label='no drug')
ax[0].plot(Fb,h_mthfd2,'--',color='#9467bd',lw=2,label='MTHFD2 inhibitor (upstream)')
ax[0].plot(Fb,h_mtor,'-',color='#1f77b4',lw=2.4,label='mTOR inhibitor 50% (downstream)')
marks(ax[0]); ax[0].set_xscale('log')
ax[0].set(xlabel='blood formate (mg/L)  [methanol range]',ylabel='net retinal health')
ax[0].legend(frameon=False,fontsize=7,loc='lower left')

# (B) dose-response at lethal methanol (500 mg/L)
inh=np.linspace(0,0.85,25); hh=[]
for I in inh:
    s,_=E.eye_state_ff(5.0,500/46.03,mtor_mult=1-I,y0=None); hh.append(H(s['eye_ins_prod']))
hh=np.array(hh); k=hh.argmax()
ax[1].axhspan(0.40,0.43,color='#2ca02c',alpha=0.08)
ax[1].plot(inh*100,hh,'-o',ms=3,color='#1f77b4')
ax[1].plot(inh[k]*100,hh[k],'*',ms=14,color='#1f77b4',zorder=5)
ax[1].annotate('therapeutic\nwindow',(inh[k]*100,hh[k]),(inh[k]*100+3,0.33),fontsize=7,color='#1f77b4')
ax[1].annotate('over-inhibition\n-> deficiency',(85,hh[-1]),(52,0.29),fontsize=7,color='#555',
               arrowprops=dict(arrowstyle='->',color='#555',lw=0.8))
ax[1].set(xlabel='mTOR (insulin-synthesis) inhibition (%)',ylabel='net retinal health')

# (C) eye insulin: mTOR inhibitor clips the excess, spares the baseline
ax[2].axhspan(Sopt,2.7,color='#d62728',alpha=0.06); ax[2].axhspan(1.5,Sopt,color='#1f77b4',alpha=0.06)
ax[2].plot(Fb,ins_no,'-',color='#7f7f7f',lw=2,label='no drug')
ax[2].plot(Fb,ins_mtor,'-',color='#1f77b4',lw=2.4,label='mTOR inhibitor 50%')
ax[2].axhline(Sopt,ls=':',color='#2ca02c',lw=1.2); ax[2].text(0.5,Sopt+0.03,'optimum',fontsize=7,color='#2ca02c')
marks(ax[2]); ax[2].set_xscale('log')
ax[2].set(xlabel='blood formate (mg/L)  [methanol range]',ylabel='eye insulin (a.u.)')
ax[2].legend(frameon=False,fontsize=7,loc='upper left')

panel(ax[0],'a','health vs methanol'); panel(ax[1],'b','dose-response'); panel(ax[2],'c','eye insulin')
fig.savefig('fig21_mtor_inhibitor.pdf')
print('methanol 500 mg/L health: no drug %.3f, MTHFD2i %.3f, mTORi50%% %.3f'%(h_no[-1],h_mthfd2[-1],h_mtor[-1]))
print('dose-response peak: %d%% inhibition -> health %.3f'%(int(inh[k]*100),hh[k]))
print('physiological (Fb=5) health: no drug %.3f -> mTORi %.3f (minimal harm)'%(
    H(E.eye_state_ff(5.0,5/46.03,y0=None)[0]['eye_ins_prod']),
    H(E.eye_state_ff(5.0,5/46.03,mtor_mult=0.5,y0=None)[0]['eye_ins_prod'])))
print('wrote fig21_mtor_inhibitor.pdf')
