import matplotlib.pyplot as plt
import numpy as np


def PlotPhi(axis, x, phi, t):
    axis.plot(x, phi, label = '{0:05.2f}[s]'.format(t))
    axis.set_xlim([np.min(x), np.max(x)])
    axis.set_ylim([-0.1, 1.1])
    axis.set_xlabel('x [m]')
    axis.set_ylabel('φ [-]')
    return


## Input Parameters
c = 0.25 # Speed [m/s]
x_left = -3 # left side of x coordinate[m]
x_right= 10 # right side of x coordinate[m]
x_num  = 150 # Number of Lattice points[-]
t      = 0 #simulation time[s]
dt     = 0.2# Simulation time step[s]
t_max  = 40 # time to finish simulation[s]
num_loop = 0 # Loop counter for simulation[-]
num_out  = 25# plot ever num_out time step[-]


## Initial Condition
x = np.linspace(x_left, x_right, x_num, endpoint=True) # x coordinate system[m]
phi_init = np.exp(-x**2) # initial distribution of phi


## Simulation
phi_old = np.copy(phi_init) # phi for nth time step
phi_new = np.copy(phi_init) # phi for n+1th time step

# figure for plot
fig = plt.figure(figsize=(12, 7))
ax_exact = fig.add_subplot(211)
ax_num   = fig.add_subplot(212)
ax_exact.set_title('Exact solution')
ax_num.set_title('Numerical solution')
fig.subplots_adjust(hspace = 0.5)

while True:
    # exact solution
    phi_exact = np.exp(-(x-c*t)**2) 

    # numerical solution
    for i in range(1, x_num):
        phi_new[i] = phi_old[i]#-Write Your Code Here-#
  
  
    # update
    phi_old[:] = phi_new[:]
    num_loop = num_loop + 1
    t      = t + dt

    # judge time step and output
    if t >= t_max:
        break
    if num_loop%num_out == 0:
        print('{0:03d}th iteration, {1:05.2f}[s]'.format(num_loop, t))
        PlotPhi(ax_exact, x, phi_exact, t)
        PlotPhi(ax_num, x, phi_new, t)

ax_exact.legend(loc='upper left')
ax_num.legend(loc='upper left')
plt.show()