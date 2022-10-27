import matplotlib.pyplot as plt
import numpy as np

## Input Parameters
m    = 1.0 # mass of a small ball [kg]
g    = 9.8 # gravitational acceleration [m/s^2]
k    = 1.0 # coefficient of air resistance [N*s/m]
v0   = 0.0 # Initial velocity [m/s]
dt   = 0.1 # time step for simulation [s]
t_max =  10 # time which simulation is stopped [s]

## Calculation
t = 0 # time [s]
i = 0  # counter
v_hist = [v0] # list to hold velocity value
t_hist = [0]  # list to hold t value

v_old = v0
while True:
  v_new = v_old+dt*(g-k/m*v_old) # recurrence formula
  t += dt
  i += 1

  if t > t_max:
    break

  v_old = v_new  
  t_hist.append(t)
  v_hist.append(v_new)

## Exact solution
t_exact = np.array(t_hist)
v_exact = m*g/k*(1-np.exp(-k*t_exact/m))

## Visualization
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(t_hist, v_hist, label='Simulation')
ax.scatter(t_exact[0::3], v_exact[0::3], label='Exact', color = "#ff7f50")
ax.set_xlabel('t [s]')
ax.set_ylabel('v [m/s]')
ax.set_title('Velocity change')
ax.legend()
plt.show()