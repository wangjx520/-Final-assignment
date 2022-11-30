# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 21:42:06 2022

@author: 20160
"""
import sympy as sy
import numpy as np

v=sy.symbols('v') #, real=True)
R=8.314

def vdw(a,b,T,P):
	
	return R*T/(v-b)-a/v/v-P

def phi(a,b,v,P,T):

	ln_phi=b/(v-b)-np.log((v-b)*P/R/T)-2*a/R/T/v
	return np.exp(ln_phi)


Tc=370
Pc=42.44
a=27*(R*Tc)**2/64/Pc
b=R*Tc/8/Pc

phi_v=1
phi_l=0
P=10
T=273.15+60
i=0
while i<5000 or abs(phi_v-phi_l)>10**(-5):
	i+=1
	V_list=sy.solve(vdw(a,b,T,P=P))
	V_list=[float(str(V).split(" ")[0]) for V in V_list]
	print(V_list)

	V_v,V_l=np.array(V_list).max(),np.array(V_list).min()
	phi_v=phi(a,b,V_v,P,T)
	phi_l=phi(a,b,V_l,P,T)
	
	
# 	if (phi_v-phi_l)*P<1:
# 		delta_P=(phi_v-phi_l)*P
# 	else:
# 		delta_P=1
	P=P+0.01
print('迭代了:',phi_v)
print('phi_v=phi_l=',phi_v)
print('丙烷在60℃的饱和压强为：',P)
	
	
	
	
