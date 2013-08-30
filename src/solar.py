from scipy.integrate.ode import *
from RungeKutta45 import RungeKutta45
import math

from numpy import *
from pylab import *

from solarBase import *


# constants : 
G   = 6.67384e-11 / 1.496597870e11**3 * 3.1556926e7**2
M   = 4*math.pi**2 / G
# initial condition lists :
m = []
r0 = []

# Jupiter :
ma = 1.8986e27
x  = 5.204267
y  = 0.0
vx = 0.0
vy = sqrt(G * M / x)
TJ = 2*math.pi*x / vy
print 'Jupiter : a = %f; T = %f ' % (x, TJ)
m.append(ma)
r0.append([x, y, vx, vy])

# create a number of asteroids with random radii :
a = 2.0  # start of asteroid belt
b = 3.5  # end of asteroid belt
r = a    # radius incrementor variable
n = 150    # number of asteroids
for i in range(n):
  mult = random.random()*(1.02-0.98) + 0.98
  r += (b-a) / n
  ma = 1
  x  = r
  y  = 0.0
  vx = 0.0
  vy = mult * sqrt(G * M / x)
  T  = 2*math.pi*x / vy
  print 'Asteroid %i : a = %f; T = %f; mult = %f' % (i+1, x, T, mult)
  m.append(ma)
  r0.append([x, y, vx, vy])




# initial conditions :
t0  = 0      # start time in years
t1  = 5000   # end time (years)
dt  = 0.25   # time step to call ode.integrate()




# set up integrator :
i = ode(motion)
i.set_integrator('RungeKutta45', atol=1e-15, rtol=1e-9)
i.set_initial_value(r0, t0)
i.set_f_params(m, G, M)

# integrate :
# format is rf[time, planet, component]
rf   = []
time = []
# add the first step to the result and time arrays
rf.append(r0)
time.append(0)
while i.t < t1:
    i.integrate(i.t+dt)
    time.append(i.t)
    rf.append(i.y)
rf   = array(rf)
time = array(time)




# calculate the conservative forces :
etot, ltot, e, l = conservative(rf, m, G, M)

# find the percent change in conservative forces :
#eChange = etot / etot[0] - 1
#lChange = ltot / ltot[0] - 1
#diff    = eChange + 2*lChange

# find the semi-major axes of all the asteroids :
a = semi_major_axes(rf, m, e, G, M)

clf()
hist(a.flatten(), bins=50, histtype='stepfilled', color='#4d4d4d')
xlabel('Semi-major Axis (AU)')
ylabel('Number of Asteroids')
savefig('hist.png')
show()

'''
# Plotting :
ion()
# figsize arguement results in 1920 x 1080 pixel output image
fig = plt.figure(figsize=(19.2,10.8))
ax1 = fig.add_subplot(221, aspect='equal')
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
colors = ['r', 'g', 'b', 'c', 'm']
for t in [len(time)]:
  # clear the axes :
  ax1.cla()
  ax2.cla()
  ax3.cla()
  ax4.cla()
  ioff()
  # plot the results :
  for p in range(len(m)):
    if t > 100:
      ax1.plot(rf[:,p][t-100:t,0], rf[:,p][t-100:t,1], 
               'k', 
               label='Planet %i' % p)
      ax1.plot(rf[:,p][t-1:t,0], rf[:,p][t-1:t,1], 
               'ko-')
    else:
      ax1.plot(rf[:,p][:t,0], rf[:,p][:t,1], 
               'k', 
               label='Planet %i' % p)
      ax1.plot(rf[:,p][t-1:t,0], rf[:,p][t-1:t,1], 
               'ko-')
  ax2.plot(time[:t], eChange[:t], 'k')
  ax3.plot(time[:t], lChange[:t], 'k')
  ax4.plot(time[:t], diff[:t], 'k')
  
  # these have to be called every time to keep visible/set :
  ax1.set_xlabel('AU')
  ax1.set_ylabel('AU')
  ax1.set_title('Paths of n-body Solar System')
  ax2.set_xlabel('Time (years)')
  ax2.set_title('% Ch. in Total Energy')
  ax3.set_xlabel('Time (years)') 
  ax3.set_title('% Ch. in Total Ang. Mom.')
  ax4.set_xlabel('Time (years)')
  ax4.set_title('Sum of E and 2*L ch.')
  # axis format is [min_x, max_x, min_y, max_y]
  ax1.axis([-6, 6, -6, 6])
  ax2.axis([0,time[-1], min(eChange), max(eChange)])
  ax3.axis([0,time[-1], min(lChange), max(lChange)])
  ax4.axis([0,time[-1], min(diff), max(diff)])
  ax1.grid()
  ax2.grid()
  ax3.grid()
  ax4.grid()
  draw()
  # to save a sequence of image files to make a movie :
  #savefig('pics/%i.png' % t)

savefig('output.png')
show()
'''





