#    Copyright (C) <2012>  <cummings.evan@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from numpy import *


def motion(t, r, m, G, M):
  '''
  PURPOSE:
    determine the acceleration of a system of n bodies which 
    exert gravitational forces on one another with a fixed Sun.
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
  for i in range(1):
    for j in range(i+1, n):

      # accelerations :
      rbtwn = r[j,:2] - r[i,:2]
      rbtwnMag = sqrt(rbtwn[0]**2 + rbtwn[1]**2)
      a[i] += (G*m[j]*rbtwn)/rbtwnMag**3
      a[j] -= (G*m[i]*rbtwn)/rbtwnMag**3
  
  # return the velocity and acceleration:
  return hstack( (r[:,2:4], array(a)) )


def motion_free_sun(t, r, m, G):
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
  OUTPUTS:
    a      - array of component velocities and accelerations in this format :
             [[v1x, v1y, a1x, a1y], ..., [vnx, vny, anx, any]]
  '''
  n = len(m)
  r = array(r)
  a = zeros(( n, 2 ))

  # fill the return array with the first term :
  for i in range(n):
    a[i] = [0, 0]
    
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
    bodies which exert gravitational forces on one another with 
    a fixed Sun.
  INPUTS:
    m        - array of masses, length = n
    r[:,i,0] - x-component distances to the Sun
    r[:,i,1] - y-component distances to the Sun
    r[:,i,2] - x-component velocities
    r[:,i,3] - y-component velocities
    G        - gravitational constant
    M        - mass of the Sun
  OUTPUTS:
    etot     - array of total energies between planets
    ltot     - array of total angular momentum between planets
    
    e        - array of individual energies.  
               format is e[time, planet] 
               where len(planet) = n and len(time) = dt
    
    l        - array of individual angular momentums. 
               format is l[time, planet] 
               where len(planet) = n and len(time) = dt
  '''
  n    = len(m)
  etot = array( [0.] * len(r) )
  ltot = array( [0.] * len(r) )
  e    = []
  l    = []
 
  # fill the return arrays with the first terms :
  for i in range(n):
    
    # individual planet kinetic energy and potential 
    # energy with respect to the sun :
    # format is r[time, planet, component]
    rmag = sqrt(r[:,i,0]**2 + r[:,i,1]**2)
    vmag = sqrt(r[:,i,2]**2 + r[:,i,3]**2)
    e.append( .5*m[i]*vmag**2 - (G*M*m[i])/rmag )
    
    # angular momentum array :
    x  = r[:,i,0]
    y  = r[:,i,1]
    vx = r[:,i,2]
    vy = r[:,i,3]
    l.append( m[i]*(x*vy - y*vx) )
  
  # update planet i with respect to every other planet j :
  for i in range(0):
    etot += e[i]
    ltot += l[i]
    for j in range(i+1, n):
      
      # vector r2 - r1 = distance vectors between planet i and j
      rbtwn = r[:,j,:2] - r[:,i,:2]
      # rbtwn[:,0] = x vectors,   rbtwn[:,1] = y vectors
      rbtwnMag = sqrt(rbtwn[:,0]**2 + rbtwn[:,1]**2)
      
      # energies :
      potBtwn = (G*m[i]*m[j]) / rbtwnMag   # potential between planets
      etot -= potBtwn
      e[i] -= potBtwn
      e[j] -= potBtwn
  
  # return the arrays :
  return array(etot), array(ltot), array(e).T, array(l).T


def conservative_free_sun(r, m, G):
  '''
  PURPOSE:
    determine the energy and angular momentum of a system of n
    bodies which exert gravitational forces on one another.
  INPUTS:
    m        - array of masses, length = n
    r[:,i,0] - x-component distances to the Sun
    r[:,i,1] - y-component distances to the Sun
    r[:,i,2] - x-component velocities
    r[:,i,3] - y-component velocities
    G        - gravitational constant
  OUTPUTS:
    etot     - array of total energies between planets
    ltot     - array of total angular momentum between planets

    
    e        - array of individual energies.  
               format is e[time, planet] 
               where len(planet) = n and len(time) = dt

    
    l        - array of individual angular momentums. 
               format is l[time, planet] 
               where len(planet) = n and len(time) = dt
  '''
  n    = len(m)
  etot = array( [0.] * len(r) )
  ltot = array( [0.] * len(r) )
  e    = []
  l    = []
 
  # fill the return arrays with the first terms :
  for i in range(n):
    
    # individual planet kinetic energy :
    vmag = sqrt(r[:,i,2]**2 + r[:,i,3]**2)
    e.append( .5*m[i]*vmag**2)
    
    # angular momentum array :
    x  = r[:,i,0]
    y  = r[:,i,1]
    vx = r[:,i,2]
    vy = r[:,i,3]
    l.append( m[i]*(x*vy - y*vx) )
  
  # update planet i with respect to every other planet j :
  for i in range(n):
    etot += e[i]
    ltot += l[i]
    for j in range(i+1, n):
      
      # vector r2 - r1 = distance vectors between planet i and j
      rbtwn = r[:,j,:2] - r[:,i,:2]
      # rbtwn[:,0] = x vectors,   rbtwn[:,1] = y vectors
      rbtwnMag = sqrt(rbtwn[:,0]**2 + rbtwn[:,1]**2)
      
      # potential energies :
      potBtwn = (G*m[i]*m[j]) / rbtwnMag   # potential between planets
      etot -= potBtwn
      e[i] -= potBtwn
      e[j] -= potBtwn
  
  # return the arrays :
  return array(etot), array(ltot), array(e).T, array(l).T


def semi_major_axes(r, m, e, G, M):
  '''
  PURPOSE:
    Calculate the semi-major axes for an array of planets.
  INPUTS:
    m - array of planet masses, length = n
    r - planet poistion and velocity array format is 
        r[time, planet, component]
        where len(time) = dt and len(planet) = n.
    e - array of planet individual energies for a given dt.  format is
        e[time, planet] where len(planet) = n and len(time) = dt
    G - gravitational constant
    M - mass of the sun
  OUTPUT:
    a - semi-major axes of each planet at each time-step
  '''
  # operate on every planet except Jupiter and every time but the last time :
  e  = e[:,1:].copy()
  m  = array(m[1:])
  x  = r[:,1:,0].copy()
  y  = r[:,1:,1].copy()
  vx = r[:,1:,2].copy()
  vy = r[:,1:,3].copy()
  
  rmag     = sqrt(x**2 + y**2)
  vmag     = sqrt(vx**2 + vy**2)
  mu       = G*(M + m)
  episilon = vmag**2/2 - mu/rmag
  a = -mu / (2*episilon)
  
  return array(a)



