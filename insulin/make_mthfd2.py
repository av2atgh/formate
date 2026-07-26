#!/usr/bin/env python3
"""Intravitreal MTHFD2 inhibitor in the eye (fig20_mthfd2.pdf).

MTHFD2 supplies the retina's local mitochondrial formate, so an intravitreal
MTHFD2 inhibitor scales down the LOCAL production (eye.eye_state_ff mthfd2 arg)
while blood-delivered formate still crosses the barrier. Three panels: (A) eye
insulin vs inhibition for an insulin-excess (lean secretion-defect) vs a
deficient (obese/low-formate) retina; (B) net retinal health for the two (benefit
is conditional); (C) methanol -- inhibition gives no protection because the toxic
formate is blood-derived and bypasses MTHFD2.
"""
import numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
os.chdir('/Users/avazquez/av2atg/formate/insulin')
from model import BetaCellModel, EyeModel
bc = BetaCellModel(); bc.calibrate_all(verbose=False)
E = EyeModel(beta=bc)
from figstyle import setup, panel
setup()
S=np.linspace(1.4,5,3000); hh=np.array([E.retinal_outcome(s)[2] for s in S]); Sopt=S[hh.argmax()]

fig,ax=plt.subplots(1,3,figsize=(12.6,4.0),constrained_layout=True)
inh=np.linspace(0,1,21)                      # MTHFD2 inhibition fraction (x100 = %)
CEX,CDEF='#d62728','#1f77b4'

# excess (lean-diabetic) and deficient (obese non-diabetic) retinas
COND={'excess (lean secretion-defect)':(6.3,1.76,CEX),
      'deficient (obese, low formate)':(4.9,1.25,CDEF)}
data={}
for lab,(glu,fbmg,c) in COND.items():
    eye=[]; hlt=[]
    for I in inh:
        s,_=E.eye_state_ff(glu,fbmg/46.03,mthfd2=1-I,y0=None)
        eye.append(s['eye_ins_prod']); hlt.append(E.retinal_outcome(s['eye_ins_prod'])[2])
    data[lab]=(np.array(eye),np.array(hlt),c)

# (A) eye insulin vs MTHFD2 inhibition
ax[0].axhspan(Sopt,2.7,color='#d62728',alpha=0.06); ax[0].axhspan(1.3,Sopt,color='#1f77b4',alpha=0.06)
ax[0].text(0.5,2.55,'insulin excess',fontsize=7,color='#d62728',ha='center')
ax[0].text(0.5,1.45,'deficiency',fontsize=7,color='#1f77b4',ha='center')
for lab,(eye,hlt,c) in data.items():
    ax[0].plot(inh*100,eye,'-o',ms=3,color=c,label=lab)
ax[0].axhline(Sopt,ls=':',color='#2ca02c',lw=1.2); ax[0].text(2,Sopt+0.02,'optimum',fontsize=7,color='#2ca02c')
ax[0].set(xlabel='intravitreal MTHFD2 inhibition (%)',ylabel='eye insulin (a.u.)',ylim=(1.35,2.65))
ax[0].legend(frameon=False,fontsize=7,loc='lower left')

# (B) net retinal health vs MTHFD2 inhibition
for lab,(eye,hlt,c) in data.items():
    ax[1].plot(inh*100,hlt,'-o',ms=3,color=c,label=lab)
    if 'excess' in lab:
        k=hlt.argmax(); ax[1].plot(inh[k]*100,hlt[k],'*',ms=13,color=c,zorder=5)
        ax[1].annotate('normalised\n(therapeutic)',(inh[k]*100,hlt[k]),(inh[k]*100-6,hlt[k]-0.055),
                       fontsize=7,color=c,ha='center')
ax[1].annotate('worsens\ndeficiency',(100,data['deficient (obese, low formate)'][1][-1]),
               (72,0.30),fontsize=7,color=CDEF,ha='center')
ax[1].set(xlabel='intravitreal MTHFD2 inhibition (%)',ylabel='net retinal health')
ax[1].legend(frameon=False,fontsize=7,loc='lower left')

# (C) methanol: no-drug vs full MTHFD2 inhibition
Fb=np.concatenate([np.linspace(0.3,5,12),np.geomspace(6,600,16)])
h1=[]; h0=[]
for fbmg in Fb:
    s1,_=E.eye_state_ff(5.0,fbmg/46.03,mthfd2=1.0,y0=None)
    s0,_=E.eye_state_ff(5.0,fbmg/46.03,mthfd2=0.0,y0=None)
    h1.append(E.retinal_outcome(s1['eye_ins_prod'])[2]); h0.append(E.retinal_outcome(s0['eye_ins_prod'])[2])
ax[2].plot(Fb,h1,'-',color='#2ca02c',lw=2,label='no drug')
ax[2].plot(Fb,h0,'--',color='#9467bd',lw=2,label='MTHFD2 inhibitor (100%)')
for v,c in [(5,'#2ca02c'),(50,'#ff7f0e'),(500,'#d62728')]: ax[2].axvline(v,ls=':',color=c,lw=1)
ax[2].annotate('no protection\n(toxin is blood-derived,\nbypasses MTHFD2)',(120,0.245),(20,0.33),
               fontsize=7,color='#555',arrowprops=dict(arrowstyle='->',color='#555',lw=0.8))
ax[2].set_xscale('log')
ax[2].set(xlabel='blood formate (mg/L)  [methanol range]',ylabel='net retinal health')
ax[2].legend(frameon=False,fontsize=7,loc='lower left')

panel(ax[0],'a','eye insulin vs inhibition'); panel(ax[1],'b','net retinal health'); panel(ax[2],'c','vs methanol: no rescue')
fig.savefig('fig20_mthfd2.pdf')
print('excess retina health: %.3f (no drug) -> peak %.3f'%(data['excess (lean secretion-defect)'][1][0],data['excess (lean secretion-defect)'][1].max()))
print('deficient retina health: %.3f (no drug) -> %.3f (full inhib)'%(data['deficient (obese, low formate)'][1][0],data['deficient (obese, low formate)'][1][-1]))
print('methanol 500mg/L health: no-drug %.3f, MTHFD2i %.3f'%(h1[-1],h0[-1]))
print('wrote fig20_mthfd2.pdf')
