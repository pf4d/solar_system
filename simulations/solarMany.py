import sys
from scipy.integrate.ode import *
from RungeKutta45 import RungeKutta45
import math

from numpy import *
from pylab import *


def motion(t, r, m, G, M):
  '''
  PURPOSE:
    determine the acceleration of a system of n bodies which 
    exert gravitational forces on one another.
  INPUTS:
    t      - time
    m      - array of masses, length = n
    r[:,0] - x-component distances to the Sun
    r[:,1] - y-component distances to the Sun
    r[:,2] - x-component velocities
    r[:,3] - y-component velocities
    G      - gravitational constant
    M      - mass of the Sun
  OUTPUTS:
    a      - array of component velocities and accelerations in this format :
             [[v1x, v1y, a1x, a1y], ..., [vnx, vny, anx, any]]
  '''
  n = len(m)
  r = array(r)
  a = zeros(( n, 2 ))

  # fill the return array with the first term :
  for i in range(n):

    rmag = sqrt(r[i,0]**2 + r[i,1]**2)

    # accelerations :
    ax = -(G*M*r[i,0])/rmag**3
    ay = -(G*M*r[i,1])/rmag**3
    a[i] = [ax, ay]
    
  # update planet i with respect to every other planet j :
  for i in range(n):
    for j in range(i+1, n):

      # accelerations :
      rbtwn = r[j,:2] - r[i,:2]
      rbtwnMag = sqrt(rbtwn[0]**2 + rbtwn[1]**2)
      a[i] += (G*m[j]*rbtwn)/rbtwnMag**3
      a[j] -= (G*m[i]*rbtwn)/rbtwnMag**3
  
  # return the velocity and acceleration:
  return hstack( (r[:,2:4], array(a)) )



def conservative(r, m, G, M):
  '''
  PURPOSE:
    determine the energy and angular momentum of a system of n
    bodies which exert gravitational forces on one another.
  INPUTS:
    m           - array of masses, length = n
    r[:,i][:,0] - x-component distances to the Sun
    r[:,i][:,1] - y-component distances to the Sun
    r[:,i][:,2] - x-component velocities
    r[:,i][:,3] - y-component velocities
    G           - gravitational constant
    M           - mass of the Sun
  OUTPUTS:
    etot        - array of total energies between planets
    ltot        - array of total angular momentum between planets
    e           - array of individual energies for time step, length n
    l           - array of individual momentums for time step, length n
  '''
  n    = len(m)
  etot = array( [0.] * len(rf) )
  ltot = array( [0.] * len(rf) )
  e    = []
  l    = []
 
  # fill the return arrays with the first terms :
  for i in range(n):
    
    # individual planet kinetic energy and potential 
    # energy with respect to the sun :
    rmag = sqrt(r[:,i][:,0]**2 + r[:,i][:,1]**2)
    vmag = sqrt(r[:,i][:,2]**2 + r[:,i][:,3]**2)
    e.append( .5*m[i]*vmag**2 - (G*M*m[i])/rmag )
    
    # angular momentum array :
    x  = r[:,i][:,0]
    y  = r[:,i][:,1]
    vx = r[:,i][:,2]
    vy = r[:,i][:,3]
    l.append( m[i]*(x*vy - y*vx) )
  
  # update planet i with respect to every other planet j :
  for i in range(1):
    etot += e[i]
    ltot += l[i]
    for j in range(i+1, n):
      
      # vector r2 - r1 = distance vectors between planet i and j
      rbtwn = r[:,j][:,:2] - r[:,i][:,:2]
      # rbtwn[:,0] = x vectors,   rbtwn[:,1] = y vectors
      rbtwnMag = sqrt(rbtwn[:,0]**2 + rbtwn[:,1]**2)
      
      # energies :
      potBtwn = (G*m[i]*m[j]) / rbtwnMag   # potential between planets
      etot -= potBtwn
      e[i] -= potBtwn
      e[j] -= potBtwn
  
  # return the arrays :
  return array(etot), array(ltot), array(e), array(l)



# constants : 
G   = 6.67384e-11 / 1.496597870e11**3 * 3.1556926e7**2
M   = 4*math.pi**2 / G

# Jupiter :
m1  = 1.8986e27
x1  = 5.204267
y1  = 0.0
vx1 = 0.0
vy1 = sqrt(G * M / x1)
T1  = 2*math.pi*x1 / vy1
print T1, 'jup', x1

# asteroid with 1/2 period of Jupiter :
m2  = 5e10
x2  = ( ((T1/2)**2 * G * M)/(4*math.pi**2) )**(1/3.)
y2  = 0.0
vx2 = 0.0
vy2 = sqrt(G * M / x2)
T2  = 2*math.pi*x2 / vy2
print T2/T1, 1/2., x2

# asteroid with 3/7 period of Jupiter :
m3  = 5e10
x3  = ( ((3*T1/7)**2 * G * M)/(4*math.pi**2) )**(1/3.)
y3  = 0.0
vx3 = 0.0
vy3 = sqrt(G * M / x3)
T3  = 2*math.pi*x3 / vy3
print T3/T1, 3/7., x3

# asteroid with 2/5 period of Jupiter :
m4  = 5e10
x4  = ( ((2*T1/5)**2 * G * M)/(4*math.pi**2) )**(1/3.)
y4  = 0.0
vx4 = 0.0
vy4 = sqrt(G * M / x4)
T4  = 2*math.pi*x4 / vy4
print T4/T1, 2/5., x4

# asteroid with 2/3 period of Jupiter :
m5  = 5e10
x5  = ( ((2*T1/3)**2 * G * M)/(4*math.pi**2) )**(1/3.)
y5  = 0.0
vx5 = 0.0
vy5 = sqrt(G * M / x5)
T5  = 2*math.pi*x5 / vy5
print T5/T1, 2/3., x5

# initial conditions :
m2=m3=m4=m5=1.0
m   = [m1, m2, m3, m4, m5]
r0  = [[x1, y1, vx1, vy1],
       [x2, y2, vx2, vy2],
       [x3, y3, vx3, vy3],
       [x4, y4, vx4, vy4],
       [x5, y5, vx5, vy5]]
t0  = 0     # start time in years
t1  = 500   # end time (years)
dt  = 0.25  # time step to call ode.integrate()

# set up integrator :
i = ode(motion)
i.set_integrator('RungeKutta45')
i.set_initial_value(r0, t0)
i.set_f_params(m, G, M)

# integrate :
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
eChange = etot / etot[0] - 1
lChange = (ltot - ltot[0]) / ltot[0]
diff    = eChange + 2*lChange


# Plotting :
ion()
# figsize arguement results in 1920 x 1080 pixel output image
fig = plt.figure(figsize=(19.2,10.8))
ax1 = fig.add_subplot(111, aspect='equal')
#ax2 = fig.add_subplot(222)
#ax3 = fig.add_subplot(223)
#ax4 = fig.add_subplot(224)
colors = ['r', 'g', 'b', 'c', 'm']
for t in range(len(time)):
  # clear the axes :
  ax1.cla()
#  ax2.cla()
#  ax3.cla()
#  ax4.cla()
  ioff()
  # plot the results :
  for p in range(len(m)):
    if t > 100:
      ax1.plot(rf[:,p][t-100:t,0], rf[:,p][t-100:t,1], 
               colors[p], 
               label='Planet %i' % p)
      ax1.plot(rf[:,p][t-1:t,0], rf[:,p][t-1:t,1], 
               colors[p] + 'o-')
    else:
      ax1.plot(rf[:,p][:t,0], rf[:,p][:t,1], 
               colors[p], 
               label='Planet %i' % p)
      ax1.plot(rf[:,p][t-1:t,0], rf[:,p][t-1:t,1], 
               colors[p] + 'o-')
#  ax2.plot(time[:t], eChange[:t], 'k')
#  ax3.plot(time[:t], lChange[:t], 'k')
#  ax4.plot(time[:t], diff[:t], 'k')
  
  # these have to be called every time to keep visible/set :
  ax1.set_xlabel('AU')
  ax1.set_ylabel('AU')
  ax1.set_title('Paths of n-body Solar System')
#  ax2.set_xlabel('Time (years)')
#  ax2.set_title('% Ch. in Total Energy')
#  ax3.set_xlabel('Time (years)') 
#  ax3.set_title('% Ch. in Total Ang. Mom.')
#  ax4.set_xlabel('Time (years)')
#  ax4.set_title('Sum of E and 2*L ch.')
  # axis format is [min_x, max_x, min_y, max_y]
  ax1.axis([-6, 6, -6, 6])
#  ax2.axis([0,time[-1], min(eChange), max(eChange)])
#  ax3.axis([0,time[-1], min(lChange), max(lChange)])
#  ax4.axis([0,time[-1], min(diff), max(diff)])
  ax1.grid()
#  ax2.grid()
#  ax3.grid()
#  ax4.grid()
  draw()
  # to save a sequence of image files to make a movie :
  #savefig('pics/%i.png' % t)

#ax1.title('Paths of Two-Planet Solar System')
savefig('output.png')
show()






