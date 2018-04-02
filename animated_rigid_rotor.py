import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(xlim=(0, 2*np.pi), ylim=(-5, 5))
line, = ax.plot([], [], lw=2)

### array of angle values
x = np.linspace(0,np.pi*2,500)


### pass array of angles (phi), central angle (phi0), spread (sig), and mean
### angular momentum (k0)
def Gauss_Packet(phi, phi0, sig, k0):
    ci = 0.+1j
    pre = 1./(sig*np.sqrt(2.*np.pi))
    psi_x = pre*np.exp(-0.5*((phi-phi0)/sig)**2)*np.exp(ci*k0*phi)
    return psi_x

def rigid_eigenfunc (phi, ml):
   ci = 0. + 1j
   return (1/np.sqrt(np.pi*2))*np.exp(ci*ml*phi)


def rigid_energy (ml):
    return ml*ml/(2*mu*r*r)

def rigid_time (ml, t):
    ci = 0. + 1j
    return np.exp(-ci*rigid_energy(ml)*t)

def FourierAnalysis(phi, Psi, ml):
    cn = np.zeros(len(ml),dtype=complex)
    ### dphi can be gotten from any two elements of phi bc it is a linspace array
    dphi = phi[1] - phi[0]
    for i in range (0,len(ml)):

      som = 0+0j
      psi_i = rigid_eigenfunc (phi, ml[i])

      ### very important to take the complex conjugate of psi_i in the following integral!!!
      ### was not necessary in the PIB case because the energy eigenfunctions were real valued, but
      ### RR energy eigenfunctions are complex
      for j in range (0, len(phi)):
        som = som + np.conj(psi_i[j])*Psi[j]*dphi

      cn[i] = som

    return cn


### set up parameters for Gaussian Wavepacket in theta-space
sig = 0.1
### include quantum numbers from m=-100 to 100, probably more than necessary
nt=np.linspace(-100, 100, 201)
phi0 = np.pi/3
m1 = 1
m2 = 1
mu = (m1*m2)/(m1+m2)
r = 1.39839789
k0 = 5*np.pi
### compute wavepacket
Psi = Gauss_Packet(x, phi0, sig, k0)
### compute expansion coefficients for building wavepacket from RR energy eigenfunctions
cn = FourierAnalysis(x, Psi, nt)
### uncomment if you want to print out all expansion coefficients
#print(cn)

### build the initial expansion just to see if FourierAnalysis is working
psi_exp = np.zeros_like(Psi)

for i in range(0,len(cn)):
  psi_exp = psi_exp + cn[i]*rigid_eigenfunc(x, nt[i])

### Plot real part of Gaussian wavepacket and real part of expansion 
plt.plot(x, np.real(Psi), 'red', x, np.real(psi_exp), 'b--')
plt.show()

def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially to generate the animation
'''  #Delete the three quote marks to uncomment this animation section
def animate(i):
    
    ### Once PIB_Func and PIB_En are defined, the following
    ### code can be used to plot the time-evolution of an energy eigenfunction

    ### Define x-grid - this will be for a particle in a box of length L=30 atomic units (Bohr radii)
    ### We will represent the function with 1000 grid points (dx = 30/1000)

    ### nt is defined above, using nt instead of m
    #m = np.linspace(-5, 5, 11)
    ### Imaginary unit i
    ci = 0.+1j
    psi_t = np.zeros_like(Psi)
    for j in range(0,len(cn)):
      p1 = rigid_eigenfunc(x, nt[j])
      ft  = rigid_time(nt[j], i/100)
      psi_t = psi_t + cn[j]*p1*ft
   
    psi_t_star = np.conj(psi_t)

    y = np.real(psi_t)
    z = np.imag(psi_t)
    p = np.real(psi_t_star * psi_t)
    line.set_data(x, p)
    return line,


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=10000, interval=200, blit=True)
### uncomment to save animation as mp4 
#anim.save('pib_wp.mp4', fps=20, extra_args=['-vcodec', 'libx264'])
plt.show()

'''
