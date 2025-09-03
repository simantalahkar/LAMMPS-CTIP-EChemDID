import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
#import math
import time

Marr = np.array([1,0,1])
M = np.sum(Marr)

Marr2 = np.array([1,0,0,1])       #to verify function output with 2019sasikumar parameters
M2 = np.sum(Marr2)

Marr_Hf = np.array([1,0,0])
M_Hf = np.sum(Marr_Hf)

Marr_Ta = np.array([0,0,1])
M_Ta = np.sum(Marr_Ta)

Marr_HfTa = np.array([1,0,1])
M_HfTa = np.sum(Marr_HfTa)

Marr_HfTi = np.array([1,1,0])
M_HfTi = np.sum(Marr_HfTi)

#global re_pure
re_pure = np.array([3.036695,2.833527,2.791418])    #in angstroms for Hf Ti Ta as 0 1 2
#global fe_cross
fe_pure = np.array([2.279201,1.953310,3.091057])
#global rhoe_cross
rhoe_pure = np.array([39.447272,32.274144,31.431176])
rhos_pure = np.array([18.953257,10.692859,28.577385])
alpha_pure = np.array([7.340695,6.939830,8.259481])
beta_pure = np.array([2.986907,3.146787,4.143922])
A_pure = np.array([0.626761,0.627344,0.689443])     #in eV
B_pure = np.array([0.645988,0.672419,1.019847])
kappa_pure = np.array([0.440318,0.438699,0.176531])
lmda_pure = np.array([1.172855,1.163775,0.390894])
m_pure = np.array([38,21,17])
n_pure = np.array([21,23,21])
Fn0_pure = np.array([-4.369507,-3.093380,-5.326794])
Fn1_pure = np.array([-0.516980,-0.321162,-0.579329])
Fn2_pure = np.array([0.710178,0.496175,1.230484])
Fn3_pure = np.array([-3.142348,-2.276043,-3.516980])
F0_pure = np.array([-4.412630,-3.120445,-5.378197])
F1_pure = np.array([0,0,0])
F2_pure = np.array([1.695061,1.089671,2.310024])
F3minus_pure = np.array([-1.476944,-0.754620,0.169535])
F3plus_pure = np.array([1.534853,0.239451,2.618849])
eta_pure = np.array([0.356582,0.262594,0.805663])
Fe_pure = np.array([-4.646960,-3.315911,-5.427809])

#global re_cross

re_cross = np.zeros([4,4])                          #O as the 1st element, index 0 ???????????????
alpha_cross = np.zeros([4,4])
beta_cross = np.zeros([4,4])
A_cross = np.zeros([4,4])
B_cross = np.zeros([4,4])
kappa_cross = np.zeros([4,4])
lmda_cross = np.zeros([4,4])
m_cross = np.zeros([4,4])
n_cross = np.zeros([4,4])



i,j = 1, 2
re_cross[i,j] = 3.229927
re_cross[j,i] = re_cross[i,j]
i,j = 1, 3
re_cross[i,j] = 3.226335
re_cross[j,i] = re_cross[i,j]
i,j = 1, 0
re_cross[i,j] = 1.938218 
re_cross[j,i] = re_cross[i,j]
i,j = 2, 3
re_cross[i,j] = 5.910024 
re_cross[j,i] = re_cross[i,j]
i,j = 2, 0
re_cross[i,j] = 1.651221
re_cross[j,i] = re_cross[i,j]
i,j = 3, 0
re_cross[i,j] = 1.985520 
re_cross[j,i] = re_cross[i,j]
i,j = 0, 0
re_cross[i,j] = 2.411822 
re_cross[j,i] = re_cross[i,j]

i,j = 1,2
alpha_cross[i,j] = 7.929612
alpha_cross[j,i] = alpha_cross[i,j]
i,j = 1,3
alpha_cross[i,j] = 9.210666
alpha_cross[j,i] = alpha_cross[i,j]
i,j = 1,0
alpha_cross[i,j] = 7.899236 
alpha_cross[j,i] = alpha_cross[i,j]
i,j = 2,3
alpha_cross[i,j] = 16.068098
alpha_cross[j,i] = alpha_cross[i,j]
i,j = 2,0
alpha_cross[i,j] = 7.071813
alpha_cross[j,i] = alpha_cross[i,j]
i,j = 3,0
alpha_cross[i,j] = 6.569020
alpha_cross[j,i] = alpha_cross[i,j]
i,j = 0,0
alpha_cross[i,j] = 6.664614
alpha_cross[j,i] = alpha_cross[i,j]

i,j = 1,2
beta_cross[i,j] = 3.149661
beta_cross[j,i] = beta_cross[i,j]
i,j = 1,3
beta_cross[i,j] = 3.547917
beta_cross[j,i] = beta_cross[i,j]
i,j = 1,0
beta_cross[i,j] = 3.162444
beta_cross[j,i] = beta_cross[i,j]
i,j = 2,3
beta_cross[i,j] = 9.446161
beta_cross[j,i] = beta_cross[i,j]
i,j = 2,0
beta_cross[i,j] = 2.730021
beta_cross[j,i] = beta_cross[i,j]
i,j = 3,0
beta_cross[i,j] = 3.931190
beta_cross[j,i] = beta_cross[i,j]
i,j = 0,0
beta_cross[i,j] = 3.205781
beta_cross[j,i] = beta_cross[i,j]

i,j = 1,2
A_cross[i,j] = 0.289205
A_cross[j,i] = A_cross[i,j]
i,j = 1,3
A_cross[i,j] = 0.235810
A_cross[j,i] = A_cross[i,j]
i,j = 1,0
A_cross[i,j] = 0.983139
A_cross[j,i] = A_cross[i,j]
i,j = 2,3
A_cross[i,j] = 0.000177
A_cross[j,i] = A_cross[i,j]
i,j = 2,0
A_cross[i,j] = 1.065502
A_cross[j,i] = A_cross[i,j]
i,j = 3,0
A_cross[i,j] = 1.173747
A_cross[j,i] = A_cross[i,j]
i,j = 0,0
A_cross[i,j] = 1.787134
A_cross[j,i] = A_cross[i,j]

i,j = 1,2
B_cross[i,j] = 0.486683
B_cross[j,i] = B_cross[i,j]
i,j = 1,3
B_cross[i,j] = 0.503209
B_cross[j,i] = B_cross[i,j]
i,j = 1,0
B_cross[i,j] = 1.574250
B_cross[j,i] = B_cross[i,j]
i,j = 2,3
B_cross[i,j] = 0.006090
B_cross[j,i] = B_cross[i,j]
i,j = 2,0
B_cross[i,j] = 1.474285
B_cross[j,i] = B_cross[i,j]
i,j = 3,0
B_cross[i,j] = 1.842358
B_cross[j,i] = B_cross[i,j]
i,j = 0,0
B_cross[i,j] = 1.397276
B_cross[j,i] = B_cross[i,j]

i,j = 1,2
kappa_cross[i,j] = 0.803686
kappa_cross[j,i] = kappa_cross[i,j]
i,j = 1,3
kappa_cross[i,j] = 0.826245
kappa_cross[j,i] = kappa_cross[i,j]
i,j = 1,0
kappa_cross[i,j] = 0.464820
kappa_cross[j,i] = kappa_cross[i,j]
i,j = 2,3
kappa_cross[i,j] = -0.428474
kappa_cross[j,i] = kappa_cross[i,j]
i,j = 2,0
kappa_cross[i,j] = 0.359030
kappa_cross[j,i] = kappa_cross[i,j]
i,j = 3,0
kappa_cross[i,j] = 0.621463
kappa_cross[j,i] = kappa_cross[i,j]
i,j = 0,0
kappa_cross[i,j] = 0.438147
kappa_cross[j,i] = kappa_cross[i,j]

i,j = 1,2
lmda_cross[i,j] = 0.975948
lmda_cross[j,i] = lmda_cross[i,j]
i,j = 1,3
lmda_cross[i,j] = 1.0753
lmda_cross[j,i] = lmda_cross[i,j]
i,j = 1,0
lmda_cross[i,j] = 1.028556
lmda_cross[j,i] = lmda_cross[i,j]
i,j = 2,3
lmda_cross[i,j] = 1.009379
lmda_cross[j,i] = lmda_cross[i,j]
i,j = 2,0
lmda_cross[i,j] = 0.654934
lmda_cross[j,i] = lmda_cross[i,j]
i,j = 3,0
lmda_cross[i,j] = 0.868149
lmda_cross[j,i] = lmda_cross[i,j]
i,j = 0,0
lmda_cross[i,j] = 0.704935
lmda_cross[j,i] = lmda_cross[i,j]

i,j = 1,2
m_cross[i,j] = 20
m_cross[j,i] = m_cross[i,j]
i,j = 1,3
m_cross[i,j] = 18
m_cross[j,i] = m_cross[i,j]
i,j = 1,0
m_cross[i,j] = 20
m_cross[j,i] = m_cross[i,j]
i,j = 2,3
m_cross[i,j] = 20
m_cross[j,i] = m_cross[i,j]
i,j = 2,0
m_cross[i,j] = 20
m_cross[j,i] = m_cross[i,j]
i,j = 3,0
m_cross[i,j] = 20
m_cross[j,i] = m_cross[i,j]
i,j = 0,0
m_cross[i,j] = 20
m_cross[j,i] = m_cross[i,j]

i,j = 1,2
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]
i,j = 1,3
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]
i,j = 1,0
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]
i,j = 2,3
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]
i,j = 2,0
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]
i,j = 3,0
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]
i,j = 0,0
n_cross[i,j] = 20
n_cross[j,i] = n_cross[i,j]

F_rhoe_oxygen = np.array([[-2.09295,4.26048,6.4839,0.13047,59.777080],
                     [-1.13802,7.72000,9.02556,0,70.325976],
                     [0.81325,13.01013,8.35517,0,84.894593],
                     [2.40945,16.30812,7.11618,0,94.647314]])
                     
F_rhoe_oxygen2 = np.array([[-1.385712,-2.753538,-1.349940,0.017886,58.452061],
                     [-1.797200,-2.247050,3.308021,0,68.767131],
                     [-2.093336,-1.164633,5.082761,0,81.249210],
                     [0,0,0,0,1]])                                  ## Parameters from 2019sasikumar
                     
F_rhoe_oxygen3 = np.array([[-2.09295,4.26048,6.4839,0.13047,59.777080],
                     [-1.13802,7.72000,9.02556,0,70.325976],
                     [0.81325,13.01013,8.35517,0,84.894593],
                     [2.40945,16.30812,7.11618,0,94.647314],
                     [5.30016,20.90029,3.35381,0,110.459837],
                     [8.21773,25.16077,7.49615,0,125.491179]])

###################

def rho_ij(fe,beta,re,lmda,n,r):
  rho = fe*np.exp(-beta*(r/re - 1))/(1 + (r/re - lmda)**n)
  return rho


dr = 0.3789308713749051E-02
Nr = 2000
stopr = Nr*dr
r = np.arange(0,stopr,dr)
#print(r)
#print(len(r))

rho_metalarr = []
for i in range(3):
  rho_metalarr.append(rho_ij(fe_pure[i],beta_pure[i],re_pure[i],lmda_pure[i],n_pure[i],r))
  print('rho_metalarr shape and last element:',np.shape(rho_metalarr),rho_metalarr[-1][-1])

filename = 'rho_Hf.dat'
rho_Hf_formatted = np.array(rho_metalarr[0]).reshape(-1,5)
np.savetxt(filename, rho_Hf_formatted, delimiter='  ')

filename = 'rho_Ta.dat'
rho_Ta_formatted = np.array(rho_metalarr[2]).reshape(-1,5)
np.savetxt(filename, rho_Ta_formatted, delimiter='  ')

fe_oxygen = 1.383394   ## from 2019sasikumar

## the cross term parameters are likely not to be used for Oxygen
## electronic density calculation as even originally Zhou et al suggested different optimizable parameters

#rho_oxygenarr = []
#rho_oxygenarr.append(rho_ij(fe_oxygen,beta_cross[3,3],re_cross[3,3],lmda_cross[3,3],n_cross[3,3],r))
#print(np.shape(rho_oxygenarr),rho_oxygenarr[-1][-1])

## function output with 2019sasikumar parameters
gamma2 = 1.568833
nu2 = 1.013103
re2 = 3.275047 
rho_oxygenarr2 = rho_ij(fe_oxygen,gamma2,re2,nu2,20,r)
print('rho_oxygenarr',np.shape(rho_oxygenarr2),rho_oxygenarr2[-1])

re3 = 2.411822 ## tried with re from 2023wu
rho_oxygenarr3 = rho_ij(fe_oxygen,gamma2,re3,nu2,20,r)

####
filename = 'rho_O2.dat'
rho_oxygenarr2_formatted = np.array(rho_oxygenarr2).reshape(-1,5)
np.savetxt(filename, rho_oxygenarr2_formatted, delimiter='  ')

filename = 'rho_O3.dat'
rho_oxygenarr3_formatted = np.array(rho_oxygenarr3).reshape(-1,5)
np.savetxt(filename, rho_oxygenarr3_formatted, delimiter='  ')

## Mixing the parameters for pairwise terms in this calculation even if from same study also gives way off results 
#rho_oxygenarr3 = []
#beta3 = 2.557467
#lmda3 = 0.710976
#rho_oxygenarr3.append(rho_ij(fe_oxygen,beta3,re2,lmda3,20,r))
#print(np.shape(rho_oxygenarr3),rho_oxygenarr3[-1][-1])


################


def F_term(F,rho_const,index,rho):
  F_index = F*(rho/rho_const - 1)**index
  return F_index


def F_parta(rho_c,F0c,F1c,F2c,F3c,rho):
  F_p1 = F_term(F0c,rho_c,0,rho) + F_term(F1c,rho_c,1,rho) + F_term(F2c,rho_c,2,rho) + F_term(F3c,rho_c,3,rho)
  return F_p1
  
def F_partb(rhos,Fe,eta,rho):
  F_p4 = Fe*(1-eta*np.log(rho/rhos))*(rho/rhos)**eta
  return F_p4


def F_metal(rhoe,rhos,eta,Fe,Fn0,Fn1,Fn2,Fn3,F0,F1,F2,F3m,F3p,rho):
  rhon = 0.85*rhoe
  rho0 = 1.15*rhoe
  F_m = []
  for i in range(Nrho):
    if rho[i]<rhon:
      F_m.append(F_parta(rhon,Fn0,Fn1,Fn2,Fn3,rho[i]))
    elif (rho[i]>=rhon) and (rho[i]<rhoe):
      F_m.append(F_parta(rhoe,F0,F1,F2,F3m,rho[i]))
    elif (rho[i]>=rhoe) and (rho[i]<rho0):
      F_m.append(F_parta(rhoe,F0,F1,F2,F3p,rho[i]))
    elif (rho[i]>=rho0):
      F_m.append(F_partb(rhos,Fe,eta,rho[i]))
  return F_m
  

drho = 0.8128985017538071E-01
Nrho = 2000
stoprho = Nrho*drho
rho = np.arange(0,stoprho,drho)

F_metalarr = []
for i in range(3):
  F_metalarr.append(F_metal(rhoe_pure[i],rhos_pure[i],eta_pure[i],Fe_pure[i],Fn0_pure[i],Fn1_pure[i],Fn2_pure[i],
                    Fn3_pure[i],F0_pure[i],F1_pure[i],F2_pure[i],F3minus_pure[i],F3plus_pure[i],rho))
  print('F_metalarr:',np.shape(F_metalarr),F_metalarr[-1][-1])

filename = 'F_Hf.dat'
F_Hf_formatted = np.array(F_metalarr[0]).reshape(-1,5)
np.savetxt(filename, F_Hf_formatted, delimiter='  ')

filename = 'F_Ta.dat'
F_Ta_formatted = np.array(F_metalarr[2]).reshape(-1,5)
np.savetxt(filename, F_Ta_formatted, delimiter='  ')

####################

#def F_o_part(F_rhoe_i,rho):
#  F_part = F_term(F_rhoe_i[0],F_rhoe_i[4],0,rho) + \
#        F_term(F_rhoe_i[1],F_rhoe_i[4],1,rho) + \
#        F_term(F_rhoe_i[2],F_rhoe_i[4],2,rho) + \
#        F_term(F_rhoe_i[3],F_rhoe_i[4],3,rho)
#  return F_part

def F_o(F_rhoe,rhoarr,M):
  F_oarr = []
  for rho in rhoarr:
    if rho < F_rhoe[0,4]:
#    print('M=0')
      i = 0
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= F_rhoe[0,4]) and (((rho < 0.5*(F_rhoe[1,4]+F_rhoe[2,4])) and M>1) or M==1):
#      print('M=1')
      i = 1
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= 0.5*(F_rhoe[1,4]+F_rhoe[2,4])) and (((rho < 0.5*(F_rhoe[2,4]+F_rhoe[3,4])) and M>2) or M==2):
#     print('M=2')
      i = 2
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= 0.5*(F_rhoe[2,4]+F_rhoe[3,4])) and M==3:
#      print('M=3')
      i = 3
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
  return F_oarr

def F_o3(F_rhoe,rhoarr):
  F_oarr = []
  for rho in rhoarr:
    if rho < F_rhoe[0,4]:
      i = 0
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= F_rhoe[0,4]) and (rho < 0.5*(F_rhoe[1,4]+F_rhoe[2,4])):
      i = 1
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= 0.5*(F_rhoe[1,4]+F_rhoe[2,4])) and (rho < 0.5*(F_rhoe[2,4]+F_rhoe[3,4])):
      i = 2
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= 0.5*(F_rhoe[2,4]+F_rhoe[3,4])) and (rho < 0.5*(F_rhoe[3,4]+F_rhoe[4,4])):
      i = 3
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif (rho >= 0.5*(F_rhoe[3,4]+F_rhoe[4,4])) and (rho < 0.5*(F_rhoe[4,4]+F_rhoe[5,4])):
      i = 4
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
    elif rho >= 0.5*(F_rhoe[4,4]+F_rhoe[5,4]):
      i = 5
      F_oarr.append(F_parta(F_rhoe[i,4],F_rhoe[i,0],F_rhoe[i,1],F_rhoe[i,2],F_rhoe[i,3],rho))
  return F_oarr

F_oxygenarr = F_o(F_rhoe_oxygen,rho,M)
print('F_oxygenarr:',np.shape(F_oxygenarr),F_oxygenarr[-1])

print('now test 2019 sasikumar')
F_oxygenarr2 = F_o(F_rhoe_oxygen2,rho,2)
print('F_oxygenarr 19sasikumar:',np.shape(F_oxygenarr2),F_oxygenarr2[-1])

print('now test 2023wu with M=5')
F_oxygenarr3 = F_o3(F_rhoe_oxygen3,rho)
print('F_oxygenarr3:',np.shape(F_oxygenarr3),F_oxygenarr3[-1])


####
filename = 'F_O1.dat'

#O_M1 = np.append(F_oxygenarr_M1,rho_oxygen, axis = 0)
F_oxygenarr_formatted = np.array(F_oxygenarr).reshape(-1,5)
np.savetxt(filename, F_oxygenarr_formatted, delimiter='  ')

filename = 'F_O2-sasikumar.dat'

F_oxygenarr2_formatted = np.array(F_oxygenarr2).reshape(-1,5)
np.savetxt(filename, F_oxygenarr2_formatted, delimiter='  ')

filename = 'F_O3.dat'

F_oxygenarr3_formatted = np.array(F_oxygenarr3).reshape(-1,5)
np.savetxt(filename, F_oxygenarr3_formatted, delimiter='  ')

#
#exit()

##########################

## test 2019sasikumar
def cross_sasi(re,A,B,alpha,beta,kappa,lmda,m,n,r):
  Ha_argument = r/re-kappa
  Hb_argument = r/re-lmda
  #Ha_argument[Ha_argument >= 0] = 1
  
  phi = A*np.exp(-alpha*(r/re-1))/(1+(Ha_argument**m)) - \
        B*np.exp(-beta*(r/re-1))/(1+(Hb_argument**n))
        
  r_phi = r*phi
  
  return phi, r_phi


def cross_new(re,A,B,alpha,beta,kappa,lmda,m,n,r):
#  print('r phi calculation')
  Ha_argument = r/re-kappa
  Hb_argument = r/re-lmda
  #Ha_argument[Ha_argument >= 0] = 1
  Ha = np.array([1 if  i >= 0 else 0 for i in Ha_argument])
  Hb = np.array([1 if  i >= 0 else 0 for i in Hb_argument])
  
  phi = A*np.exp(-alpha*(r/re-1))/(1+(Ha_argument**m)*Ha) - \
        B*np.exp(-beta*(r/re-1))/(1+(Hb_argument**n)*Hb)
        
  r_phi = r*phi
  
  return r_phi
 
phi_sasi, r_phi_sasi = cross_sasi(3.275047,0.624804,0.755903,5.410041,2.557467,0.336468,0.710976,20,20,r)
print('phi, r_phi 19sasikumar:\n',np.shape(phi_sasi),np.shape(r_phi_sasi),phi_sasi[-1],r_phi_sasi[-1])

filename = 'rphi_OO_sasi.dat'

#O_M1 = np.append(F_oxygenarr_M1,rho_oxygen, axis = 0)
r_phi_sasi_formatted = np.array(r_phi_sasi).reshape(-1,5)
np.savetxt(filename, r_phi_sasi_formatted, delimiter='  ')


 
cross_potential = []
indices = [[0,0],[1,0],[1,1],[2,0],[2,1],[2,2],[3,0],[3,1],[3,2],[3,3]]
for i in indices:
    if i[0]!=i[1] or (i[0]==0 and i[1]==0):
        cross_potential.append(cross_new(re_cross[i[0],i[1]],A_cross[i[0],i[1]],
                            B_cross[i[0],i[1]],alpha_cross[i[0],i[1]],
                            beta_cross[i[0],i[1]],kappa_cross[i[0],i[1]],
                            lmda_cross[i[0],i[1]],m_cross[i[0],i[1]],
                            n_cross[i[0],i[1]],r))
                            
    elif i[0]==i[1] and i[0]!=0:
        cross_potential.append(cross_new(re_pure[i[0]-1],A_pure[i[0]-1],
                            B_pure[i[0]-1],alpha_pure[i[0]-1],
                            beta_pure[i[0]-1],kappa_pure[i[0]-1],
                            lmda_pure[i[0]-1],m_pure[i[0]-1],
                            n_pure[i[0]-1],r))

print('cross_potential:',np.shape(cross_potential),cross_potential[-1][-1])

filename = 'rphi_OO.dat'
r_phi_OO_formatted = np.array(cross_potential[0]).reshape(-1,5)
np.savetxt(filename, r_phi_OO_formatted, delimiter='  ')

filename = 'rphi_HfO.dat'
r_phi_HfO_formatted = np.array(cross_potential[1]).reshape(-1,5)
np.savetxt(filename, r_phi_HfO_formatted, delimiter='  ')

filename = 'rphi_HfHf.dat'
r_phi_HfHf_formatted = np.array(cross_potential[2]).reshape(-1,5)
np.savetxt(filename, r_phi_HfHf_formatted, delimiter='  ')

filename = 'rphi_TiO.dat'
r_phi_TiO_formatted = np.array(cross_potential[3]).reshape(-1,5)
np.savetxt(filename, r_phi_TiO_formatted, delimiter='  ')

filename = 'rphi_TiHf.dat'
r_phi_TiHf_formatted = np.array(cross_potential[4]).reshape(-1,5)
np.savetxt(filename, r_phi_TiHf_formatted, delimiter='  ')

filename = 'rphi_TiTi.dat'
r_phi_TiTi_formatted = np.array(cross_potential[5]).reshape(-1,5)
np.savetxt(filename, r_phi_TiTi_formatted, delimiter='  ')

filename = 'rphi_TaO.dat'
r_phi_TaO_formatted = np.array(cross_potential[6]).reshape(-1,5)
np.savetxt(filename, r_phi_TaO_formatted, delimiter='  ')

filename = 'rphi_TaHf.dat'
r_phi_TaHf_formatted = np.array(cross_potential[7]).reshape(-1,5)
np.savetxt(filename, r_phi_TaHf_formatted, delimiter='  ')

filename = 'rphi_TaTi.dat'
r_phi_TaTi_formatted = np.array(cross_potential[8]).reshape(-1,5)
np.savetxt(filename, r_phi_TaTi_formatted, delimiter='  ')

filename = 'rphi_TaTa.dat'
r_phi_TaTa_formatted = np.array(cross_potential[9]).reshape(-1,5)
np.savetxt(filename, r_phi_TaTa_formatted, delimiter='  ')

exit()

#########

F_oxygenarr_M1 = np.array(F_o(F_rhoe_oxygen,rho,1))
F_oxygenarr_M2 = np.array(F_o(F_rhoe_oxygen,rho,2))

rho_oxygen = np.array(rho_oxygenarr2)

filename = 'O_M1.dat'

print('\n',np.shape(F_oxygenarr_M1),np.shape(rho_oxygen),'\n')

O_M1 = np.append(F_oxygenarr_M1,rho_oxygen, axis = 0)

print('\n',np.shape(O_M1),'\n')
print(O_M1[0],F_oxygenarr_M1[0],O_M1[-1],rho_oxygen[-1],'\n')

O_M1_formatted = O_M1.reshape(-1,5)

print('\n',np.shape(O_M1_formatted),'\n')
print(O_M1_formatted[0:2],'\n',O_M1[0:10],'\n')

np.savetxt(filename, O_M1_formatted, delimiter='  ')

#

filename = 'O_M2.dat'
O_M2 = np.append(F_oxygenarr_M2,rho_oxygen, axis = 0)
O_M2_formatted = O_M2.reshape(-1,5)
np.savetxt(filename, O_M2_formatted, delimiter='  ')

#######

F_metalarr = np.array(F_metalarr)
rho_metalarr = np.array(rho_metalarr)

filename = 'Hf.dat'

Hf_embed = np.append(F_metalarr[0],rho_metalarr[0],axis = 0)
Hf_embed_formatted = Hf_embed.reshape(-1,5)
print('\n',np.shape(Hf_embed_formatted),'\n')
print(Hf_embed_formatted[0:2],'\n',Hf_embed[0:10],'\n')

np.savetxt(filename, Hf_embed_formatted, delimiter='  ')

#

filename = 'Ti.dat'

Ti_embed = np.append(F_metalarr[1],rho_metalarr[1],axis = 0)
Ti_embed_formatted = Ti_embed.reshape(-1,5)
print('\n',np.shape(Ti_embed_formatted),'\n')
print(Ti_embed_formatted[0:2],'\n',Ti_embed[0:10],'\n')

np.savetxt(filename, Ti_embed_formatted, delimiter='  ')

#

filename = 'Ta.dat'

Ta_embed = np.append(F_metalarr[2],rho_metalarr[2],axis = 0)
Ta_embed_formatted = Ta_embed.reshape(-1,5)

np.savetxt(filename, Ta_embed_formatted, delimiter='  ')

####

cross_potential = np.array(cross_potential)

filename = 'OHf_cross.dat'

OHf_cross = np.concatenate((cross_potential[0],cross_potential[1],cross_potential[2]),axis = 0)
OHf_cross_formatted = OHf_cross.reshape(-1,5)
print('\n',np.shape(OHf_cross_formatted),'\n')
print(OHf_cross_formatted[0:2],'\n',OHf_cross[0:10],'\n')

np.savetxt(filename, OHf_cross_formatted, delimiter='  ')

#

filename = 'OTa_cross.dat'

OTa_cross = np.concatenate((cross_potential[0],cross_potential[6],cross_potential[9]),axis = 0)
OTa_cross_formatted = OTa_cross.reshape(-1,5)
print('\n',np.shape(OTa_cross_formatted),'\n')

np.savetxt(filename, OTa_cross_formatted, delimiter='  ')

#

filename = 'OHfTi_cross.dat'

OHfTi_cross = np.concatenate((cross_potential[0],cross_potential[1],cross_potential[2],
                            cross_potential[3],cross_potential[4],cross_potential[5]),axis = 0)
OHfTi_cross_formatted = OHfTi_cross.reshape(-1,5)
print('\n',np.shape(OHfTi_cross_formatted),'\n')

np.savetxt(filename, OHfTi_cross_formatted, delimiter='  ')

#

filename = 'OHfTa_cross.dat'

OHfTa_cross = np.concatenate((cross_potential[0],cross_potential[1],cross_potential[2],
                            cross_potential[6],cross_potential[7],cross_potential[9]),axis = 0)
OHfTa_cross_formatted = OHfTa_cross.reshape(-1,5)
print('\n',np.shape(OHfTa_cross_formatted),'\n')

np.savetxt(filename, OHfTa_cross_formatted, delimiter='  ')

#

filename = 'Ti_cross.dat'

Ti_cross_formatted = cross_potential[5].reshape(-1,5)
print('\n',np.shape(Ti_cross_formatted),'\n')

np.savetxt(filename, Ti_cross_formatted, delimiter='  ')



exit()





