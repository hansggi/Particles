import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button, RadioButtons
"""
Same as particles_LJ but updates positions and animates real time
"""

def init_pos(radius,numberOfParticles):
    to_close = True
    while(to_close):
        to_close = False
        r = np.random.uniform(0,radius, numberOfParticles)
        theta = np.random.uniform(0,2*np.pi, numberOfParticles)
        pos_x = r * np.cos(theta)
        pos_y = r * np.sin(theta)
            
        for i in range(numberOfParticles):
            for j in range(i+1,numberOfParticles):
                d = np.sqrt((pos_x[i] - pos_x[j])**2 + (pos_y[i] - pos_y[j])**2)
                if d < 1:
                    to_close = True
                    break
            if to_close:
                break
    return pos_x, pos_y

def random_velocities(numberOfParticles):
    angles =  np.random.uniform(-np.pi, np.pi,numberOfParticles)
    vx = np.cos(angles) / np.sqrt(numberOfParticles / 2)  
    vy = np.sin(angles) / np.sqrt(numberOfParticles / 2)  
    return vx, vy

def vel_one_particle(numberOfParticles):
    vx = np.zeros(numberOfParticles)
    vy = np.zeros(numberOfParticles)
    angle = np.random.uniform(-np.pi,np.pi, 1)
    vx[0] = np.cos(angle) * np.sqrt(2)
    vy[0] = np.sin(angle) * np.sqrt(2)
    return vx, vy


def one_step_billiard(numberOfParticles,radius, dt, KK, pos, v ):
    
    distance = np.sqrt(pos[0,:]**2 + pos[1,:]**2)
    acceleration = np.zeros((2, numberOfParticles))
    
    #test if we are outside the circle and calculate force from wall
    acceleration[0, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[0, distance > radius] / distance[distance > radius]
    acceleration[1, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[1, distance > radius] / distance[distance > radius]
    
    for j in range(numberOfParticles):
            for k in range(j+1, numberOfParticles):
                d_jk = np.sqrt((pos[0,j] - pos[0,k])**2 + (pos[1,j] - pos[1,k])**2)
                if d_jk < 1:
                    F = 12 * epsilon * (1/d_jk**7 - 1/d_jk**13)
                    
                    acceleration[0,j] += - F * (pos[0,j] - pos[0,k]) / d_jk 
                    acceleration[1,j] += - F * (pos[1,j] - pos[1,k]) / d_jk 
                    acceleration[0,k] += F * (pos[0,j] - pos[0,k]) / d_jk 
                    acceleration[1,k] += F * (pos[1,j] - pos[1,k]) / d_jk
                    
    #update positions
    pos[0, :] =  pos[0, :] + v[0, :] * dt + (acceleration[0, :] * dt**2)/2
    pos[1, :] =  pos[1, :] + v[1, :] * dt + (acceleration[1, :] * dt**2)/2
    
    distance = np.sqrt(pos[0,:]**2 + pos[1,:]**2)
    acceleration2 = np.zeros((2, numberOfParticles))

    #test if we are outside the circle and calculate force from wall
    acceleration2[0, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[0, distance > radius] / distance[distance > radius]
    acceleration2[1, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[1, distance > radius] / distance[distance > radius]
    for j in range(numberOfParticles):
            for k in range(j+1, numberOfParticles):
                d_jk = np.sqrt((pos[0,j] - pos[0,k])**2 + (pos[1,j] - pos[1,k])**2)
                if d_jk < 1:
                    F = 12 * epsilon * (1/d_jk**7 - 1/d_jk**13)
                    
                    acceleration2[0,j] += - F * (pos[0,j] - pos[0,k]) / d_jk 
                    acceleration2[1,j] += - F * (pos[1,j] - pos[1,k]) / d_jk 
                    acceleration2[0,k] += F * (pos[0,j] - pos[0,k]) / d_jk 
                    acceleration2[1,k] += F * (pos[1,j] - pos[1,k]) / d_jk
                    
    #update velocities
    v[0, :] +=  (acceleration[0,:] + acceleration2[0,:])/2 * dt
    v[1, :] += (acceleration[1,:] + acceleration2[1,:])/2 * dt
    return pos, v

epsilon = 1
num_part = 10
num_it = 50000
r = 5
dt = 0.02
KK = 5
a = 1

pos_0 = np.zeros((2, num_part))
v_0 = np.zeros((2, num_part))
pos_0[0,:], pos_0[1,:] = init_pos(r,num_part)

v_0[0, :], v_0[1, :] = random_velocities(num_part)

pos = np.zeros((2,num_part, num_it + 1))
v = np.zeros((2,  num_part, num_it + 1))
pos[:,:, 0] = pos_0
v[:,:, 0] = v_0

def init(): #must use name init for funcAnimation. Initiates lines and circles.
    
    for j in range(num_part):
        circles[j].center = (pos[0,j,0],pos[1,j,0])
        ax.add_patch(circles[j])
    


    ang = np.linspace(0, 2 * np.pi, 100)
    line1.set_data(np.cos(ang) * r, np.sin(ang) * r)
    
    return circles


def animate(i): # amination-function. i = frame number. 
    pos[:,:,  i+1], v[:,:, i+1] = one_step_billiard(num_part, r, dt, KK, pos[:,:, i], v[:,:,i])
    if i %a ==0:
        
        for j in range(num_part):
            x1, y1 = pos[0,j, i], pos[ 1,j, i]
            circles[j].center = (x1, y1)

    return circles

circles = []
for i in range(num_part):
    circles.append(plt.Circle((0, 0), 0.5, color='r'))
        

plt.style.use('dark_background')
fig = plt.figure( figsize = (8, 8))
 
ax = plt.axes(xlim = (-1.5 * r, 1.5*r), ylim = (-1.5*r , 1.5*r))
line1, = ax.plot([], [], lw = 2, color = "c")

frames =  int( num_it / a ) 
plt.axis("off") 

anim = animation.FuncAnimation(fig, animate,  init_func= init, frames = frames, interval = 1, blit = True) #interval given in ms
plt.draw()
plt.show()


