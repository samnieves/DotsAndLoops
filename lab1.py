import numpy as np
from matplotlib import pyplot as plt



def PIB_Func(x, n, L):
        return np.sqrt(2/L)*np.sin(n*np.pi*x/L)

def Gauss_Packet(sig,x, x0,  k0):

        ci = 0 + 1j
        pre = 1/(sig*np.sqrt(2*np.pi))
        gx = np.exp(-0.5*((x-x0)/sig)**2)
        pw = np.exp(ci*k0*x)
        return pre*gx*pw

def FourierAnalysis(x, PsiX, n, L):
    cn = np.zeros(len(n),dtype=complex)
    dx = x[1]-x[0]
    for i in range (0,len(cn)):
        som = 0+0j
        psi_i = PIB_Func(x, n[i], L)

        for j in range (0, len(x)):
           som = som + psi_i[j]*PsiX[j]*dx

        cn[i] = som

    return cn

L = 500.
sig = 20.
k0 = 0.5
x0 = 200.
N = 500
x = np.linspace(0,L,500)
n = np.linspace(1, 100,100)

y = Gauss_Packet(sig, x, x0, k0)
P = np.real(np.conj(y)*y)
cn = FourierAnalysis(x, y, n, L)

psi_exp = np.zeros(len(x))
for i in range (0,len(cn)):
    psi_exp = psi_exp + cn[i]*PIB_Func(x, i+1, L)


plt.plot(x,np.real(psi_exp),'r--', x, np.real(y), 'blue')
plt.show()
