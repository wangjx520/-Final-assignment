# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 22:51:25 2022

@author: 20160
"""
import numpy as np
import sympy as sy
v=sy.symbols('v') #,real=True) #, real=True)
def vdw(a,b,T,P):
	
	return R*T/(v-b)-a/v/v-P

def ab(R,Tc,Pc):
	return 27*(R*T)**2/64/Pc,R*Tc/8/Pc

[Tc3,Pc3,A3,B3,C3]=[370, 42.44, 9.1058, 1872.46, -25.16]
[Tc4,Pc4,A4,B4,C4]=[425.2,37.9, 9.0580, 2154.90, -34.42]

R=8.314
T=273.15+20

a3,b3=ab(R,Tc3,Pc3*10**5)
a4,b4=ab(R,Tc4,Pc4*10**5)

P3_sat=np.exp(A3-B3/(T+C3))*10**5
P4_sat=np.exp(A4-B4/(T+C4))*10**5

x3=0.49
x4=0.51
y3=0.2  #随机给的初始值
y4=0.2  #随机给的初始值

P0=x3*P3_sat+x4*P4_sat
y3=x3*P3_sat/P0
y4=x4*P4_sat/P0
P=P0

i=0
while i==0 or i<50000 or abs(y3+y4-1)>10**(-5):
	i+=1
    
	a_l=x3*x3*a3 + 2*x3*x4*np.sqrt(a3*a4) + x4*x4*a4
	a_v=y3*y3*a3 + 2*y3*y4*np.sqrt(a3*a4) + y4*y4*a4
	b_l=x3*b3+x4*b4	
	b_v=y3*b3+y4*b4
	
	V_v_list=sy.solve(vdw(a_v,b_v,T,P))
	print(V_v_list)
	V_v_list=[float(str(V).split(" ")[0]) for V in V_v_list]
	V_v=np.array(V_v_list).max()
	
	V_l_list=sy.solve(vdw(a_l,b_l,T,P))
	print(V_l_list)
	V_l_list=[float(str(V).split(" ")[0]) for V in V_l_list]
	V_l=np.array(V_l_list).min()
	
	ln_phi_v3=b3/(V_v-b_v)-np.log((V_v-b_v)*P/R/T)-2*(y3*a3+y4*a4)/R/T/V_v
	ln_phi_v4=b4/(V_v-b_v)-np.log((V_v-b_v)*P/R/T)-2*(y3*a3+y4*a4)/R/T/V_v
	ln_phi_l3=b3/(V_l-b_l)-np.log((V_l-b_l)*P/R/T)-2*(y3*a3+y4*a4)/R/T/V_l
	ln_phi_l4=b4/(V_l-b_l)-np.log((V_l-b_l)*P/R/T)-2*(y3*a3+y4*a4)/R/T/V_l

	phi_v3=np.exp(ln_phi_v3)
	phi_v4=np.exp(ln_phi_v4)
	phi_l3=np.exp(ln_phi_l3)
	phi_l4=np.exp(ln_phi_l4)
	
	y3=0.3*phi_l3/phi_v3
	y4=0.7*phi_l4/phi_v4

	P=P*(y3+y4)

    
print('迭代了：',i)
print('混合物的压强为：',P)
