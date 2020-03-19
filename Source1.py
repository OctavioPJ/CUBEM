
# coding: utf-8

# In[13]:

import sys
sys.path.append('D:\\')
from PyCAR.PyCIT.FT import MeshFlux,MeshPower
import numpy as np
from PyCAR.PyCIT.FT import GroupFlux
from citvappru.SourceCAREM2_1_1 import Geometry,PrecCalc

BetaDir='D:\\Cube\\'

Flux = MeshFlux(BetaDir+'NoobCube1g.cdb.meshflux',Nodosx,Nodosy,Nodosz,Ngrupos)
    
NPRC = 6
NOGBIB = Flux.shape[-1]
#NuFis=[2.0150E-02]

lmks=[0.0124,0.0306,0.1135,0.3071,1.1905,3.1748]
Betas=[0.0002170,0.0014979,0.0013778,0.0028296,0.0009288,0.0003314,0.0071826]
BetT=Betas[-1]

#valor obtenido de NoobCube2g.cii
keff = 1.00055
keff = 1.0005462

chi_g=[0.999827,0.000173]
#NuFis=[2.8132E-02,3.6287E-02]
NuFis=[2.81322002E-03,3.62874009E-02]
#NuFis=list(map(lambda XS:XS/keff,NuFis))#esto no va multiplicado por (1-b)
NuFis


# In[7]:

BetT


# In[8]:

sum(Betas[:-1])


# In[9]:

def Initialize(iterable):
    IterM={}#matriz de betas
    for state in range(Flux.shape[0]):
        IterM[state]={}
        for nx in range(Flux.shape[1]):
            IterM[state][nx] = {}
            for ny in range(Flux.shape[2]):
                IterM[state][nx][ny] = {}
                for nz in range(Flux.shape[3]):
                    IterM[state][nx][ny][nz]= iterable
    return IterM

BetaM = Initialize(Betas)


# In[10]:

LmkM = Initialize(lmks)


# In[11]:

NuFisM = Initialize({group:NuFis[group] for group in range(Flux.shape[-1])})
NuFisM


# calculo de precursores en equilibrio

# # LA FUENTE SE DA EN NEUTRONES POR SEGUNDO,
# # NO EN NEUTRONES/CM3
# 
# SE SUBESTIMA LA CANTIDAD DE PRODUCCIONES

# In[15]:

GFlux=GroupFlux(BetaDir+'NoobCube2g.cdb.flux')
GFlux.shape


# In[16]:

Vmesh = (6*6*6)
Vcell = 6*Vmesh


# In[17]:

Vmesh


# In[26]:

C0 = {}
for state in range(Flux.shape[0]):
    C0[state]={}
    for nx in range(Flux.shape[1]):
        C0[state][nx] = {}
        for ny in range(Flux.shape[2]):
            C0[state][nx][ny] = {}
            for nz in range(Flux.shape[3]):
                
                Nu_Flux = sum([NuFisM[state][nx][ny][nz][group]*Flux[state,nx,ny,nz,group]                               for group in range(Flux.shape[-1])])
                
                Bet_k = BetaM[state][nx][ny][nz]
                Lamb_k = LmkM[state][nx][ny][nz]
                              
                C0[state][nx][ny][nz] = [Bet_k[prec]*Nu_Flux*Vmesh/Lamb_k[prec]                                              for prec in range(NPRC)]

alf=np.array(list(C0[0][22][11].values()))
lmkr=np.array(lmks)
(lmkr*alf).sum(axis=1)


# FUENTE FOR DUMMIES, confirmado que devuelve los 100 MW

# In[22]:

FissRate = NuFis[0]*Flux[...,0]+NuFis[1]*Flux[...,1]
Fissions = FissRate*Vmesh
Precs = BetT*Fissions
NewP = Fissions.sum()
'{:.5e}'.format(NewP)


# In[19]:

Fissions.shape


# In[217]:

PowS=7.49346E+03

SourceP=4.074395E+12
TotalSL=5.65582E+14
TotalSP=5.61508E+14

#calculado con las nufisiones divididas por keff y multiplicadas por 1-bet
PowE=1.00000E+06
#TotalEP=7.49340E+16#subcritico
TotalEP=7.55175E+16#supercritico
TotalEL=7.54763E+16


# In[218]:

PowS/PowE


# In[219]:

SourceP/TotalSP


# In[220]:

BetT/(1-BetT)


# In[221]:

TotalSL/(TotalSP+SourceP)


# In[222]:

TotalSP/TotalEP


# In[223]:

'{:.5E}'.format(BetT*TotalEP)


# In[27]:

alf


# In[228]:

alf.shape


# AQUI HAY ALGUN KILOMBO, NO COMPENSA LA FUENTE

# In[29]:

import matplotlib.pyplot as plt
alf=lmkr*alf
plt.plot(range(32),alf.sum(axis=1),'-')
plt.plot(range(32),Precs[0,22,11,:],'o')
plt.show()


# In[230]:

Q={}
EQUILIBRIUM=True
if EQUILIBRIUM:
    C = C0
#else:
#    C = C0 + dCdt*dt
#Q0={}
for group in range(Flux.shape[-1]):
    Q[group]={}
    for state in range(Flux.shape[0]):
        Q[group][state]={}
        for nx in range(Flux.shape[1]):
            Q[group][state][nx] = {}
            for ny in range(Flux.shape[2]):
                Q[group][state][nx][ny] = {}
                for nz in range(Flux.shape[3]):
                    _Lmk= LmkM[state][nx][ny][nz]
                    _C  = C[state][nx][ny][nz]
                    #_Q[state][nx][ny][nz] = sum([ _Lmk[prc]*_C[prc] for prc in range(NPRC)])
                    #de la fuente for dummies
                    Q[group][state][nx][ny][nz] = chi_g[group]*Precs[state,nx,ny,nz]
Q


#  NUMBER OF---COLUMNS, ROWS, PLANES, GROUPS, UPSCAT, DOWNSCAT, REGIONS, AND ZONES      44  22  32   1   0   0  29040   1999
# 
 > 0 Specifies that the whole space-energy fixed source is provided on an I/O device (logical 17) in the form:
 starting with the highest energy, the source (neutrons/sec) at each space point is given by,
     DO 1 K = 1,KMAX
     READ(17)(((S(J,I,KB),J=1,JMAX),I=1,IMAX),
         KB=1,KBMAX)
1    CONTINUE
 
 where J indexes along rows, I along columns, and KB through the planes.
 No more data cards are required.Q0={0: Q,\
    1: [[[[0.0 if Flux[state,nx,ny,nz,1] != 0 else 0.0\
                    for nz in range(Flux.shape[3])]\
                    for ny in range(Flux.shape[2])]\
                    for nx in range(Flux.shape[1])] \
                    for state in range(Flux.shape[0])]}

#se está suponiendo una fuente que aparece en un solo grupo
# In[231]:

Q0=Q


# In[232]:

with open('source.dat','w') as fod:
    for group in range(NOGBIB):
        for state in range(Flux.shape[0]):
            for nz in range(Flux.shape[3]):
                for ny in range(Flux.shape[2]):
                    for nx in range(Flux.shape[1]):
                        fod.write('{:12.5E}'.format(Q0[group][state][nx][ny][nz])) 
        fod.write('\n')


# In[233]:

get_ipython().system('cat source.dat')


# In[234]:

get_ipython().system('cp -f source.dat Source_To_LU17/')

