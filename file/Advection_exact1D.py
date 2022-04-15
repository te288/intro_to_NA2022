# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

## Input Parameters
c = 1
t_init = 0
t_max = 15
dt    = 0.2
x = np.linspace(-3, 10, 200, endpoint=True)
t = np.arange(t_init, t_max+dt, dt, )

## Plot Data
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for i in t:
  ax.plot(x, np.exp(-(x-c*i)**2), color = '#4169e1')
  ax.set_ylim(-0.1,1.1)
  ax.set_xlim(-3, 10)
  ax.set_title('t=%05.1f'%(i)+'[s]')
  # fig.savefig('Advection_exact%05.1f'%(i)+'.png') # Save figure
  # ax.cla()

plt.show()