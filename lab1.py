import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 500), ylim=(-0.03, 0.03))
line, = ax.plot([], [], lw=2)


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

def PIB_En(n, L):
    En = (n*n * np.pi*np.pi)/(2*L*L)
    return En

### Give the quantum number and the current time, evaluate the time-dependent part of the wavefunction at current time
### and return its value
def PIB_Time(n, L, t):
    E = PIB_En(n, L)
    ci = 0.+1j
    phi_n_t = np.exp(-1*ci*E*t)
    ### Write code here to define phi_n_t
    return phi_n_t


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

def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially to generate the animation
def animate(i):
    
    ### Once PIB_Func and PIB_En are defined, the following
    ### code can be used to plot the time-evolution of an energy eigenfunction

    ### Define x-grid - this will be for a particle in a box of length L=30 atomic units (Bohr radii)
    ### We will represent the function with 1000 grid points (dx = 30/1000)
    L = 500.
    x = np.linspace(0, L, 2000)

    ### Imaginary unit i
    ci = 0.+1j
    fwhm = 7*np.pi/L
    k0 = 5*np.pi/L
    psi_t = np.zeros(len(x),dtype=complex)
    #psi = PIB_Func(x, 2, L)
    #ft = PIB_Time(2, L, i)
    #psi = PIB_Func(x, 10, L)
    #ft = PIB_Time(10, L, 4*i)
    #psi_t = psi*ft
    for j in range(0,len(cn)):
      psi = PIB_Func(x, nt[j], L) 
      ft  = PIB_Time(nt[j], L, 4*i)
      psi_t = psi_t + cn[j]*psi*ft
   
    psi_t_star = np.conj(psi_t)

    y = np.real(psi_t)
    z = np.imag(psi_t)
    p = np.real(psi_t_star * psi_t)
    line.set_data(x, y)
    return line,


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=10000, interval=20, blit=True)
### uncomment to save animation as mp4 
#anim.save('pib_wp.mp4', fps=20, extra_args=['-vcodec', 'libx264'])
plt.show()


#### Static plot
#plt.plot(x, np.real(Psi), 'b--', x, np.real(psi_exp), 'red', x, P, 'black')
#plt.show()



#plt.plot(x,np.real(psi_exp),'r--', x, np.real(y), 'blue')
#plt.show()
