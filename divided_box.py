import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
"""
animation of particles in divided box. higher potential on right side of box
"""

def init_pos(numberOfParticles):
    """
    initial positions are taken from uniform distribution
    """
    pos_x, pos_y = np.random.uniform(-Lx, Lx, numberOfParticles), np.random.uniform(-Ly, Ly, numberOfParticles)
    
    return pos_x, pos_y

def random_velocities(numberOfParticles):
    """
    initial velocities are taken from uniform distribution and normalized
    """
    vx_old = np.random.uniform(-1,1, numberOfParticles) / np.sqrt(numberOfParticles / 2)
    vy_old = np.random.uniform(-1,1, numberOfParticles) / np.sqrt(numberOfParticles / 2)
    vx = vx_old / np.sqrt(vx_old**2 + vy_old**2) / 10
    vy = vy_old / np.sqrt(vx_old**2 + vy_old**2) / 10
    return vx, vy


def part_in_box(numberOfParticles, numberOfIterations, Lx, Ly, dt, KK, epsilon_step = 0.1, vel_dist = random_velocities):
    pos = np.zeros((2, numberOfParticles, numberOfIterations+1))
    pos[0,:,0], pos[1,:,0] = init_pos(numberOfParticles)
    vx, vy = vel_dist(numberOfParticles)
    vx_list,vy_list = np.zeros((numberOfIterations,numberOfParticles,)), np.zeros((numberOfIterations,numberOfParticles,))

    for i in range(numberOfIterations):
        vx_list[i] = vx
        vy_list[i] = vy

        acceleration = np.zeros((2, numberOfParticles))
        #test if we are outside the square and calculate force from wall
        acceleration[0, np.absolute( pos[0,:,i]) > Lx] += -KK *pos[0,:,i][np.absolute( pos[0,:,i]) > Lx]
        acceleration[1, np.absolute( pos[1,:,i]) > Ly] += -KK *pos[1,:,i][np.absolute( pos[1,:,i]) > Ly]

        #constant force to the left if part is in small region in the middle ( "step potential"):
        acceleration[0, ( Lx / 10 > pos[0,:,i]) > 0] += -epsilon_step 
                    
        #update positions
        pos[0, :, i+1] =  pos[0, :, i] + vx * dt + (acceleration[0, :] * dt**2) / 2
        pos[1, :, i+1] =  pos[1, :, i] + vy * dt + (acceleration[1, :] * dt**2) / 2
        
        #do it again
        acceleration2 = np.zeros((2, numberOfParticles))
        acceleration2[0, np.absolute( pos[0,:,i]) > Lx] += -KK *pos[0,:,i][np.absolute( pos[0,:,i]) > Lx]
        acceleration2[1, np.absolute(pos[1,:,i]) > Ly] += -KK *pos[1,:,i][np.absolute( pos[1,:,i]) > Ly]
        acceleration2[0, (Lx / 10 > np.absolute(pos[0,:,i]) ) > 0] += -epsilon_step 
              
        #update velocities
        vx += (acceleration[0,:] + acceleration2[0,:])/2 * dt
        vy += (acceleration[1,:] + acceleration2[1,:])/2 * dt

    return pos

num_part = 30
num_it = 30000
dt = 0.01
Lx, Ly = 10, 5
pos = part_in_box(num_part, num_it, Lx, Ly, dt = 0.01, KK = 0.2, epsilon_step = 0.05)

#plot every a'th value
a = 10
pos = pos[:, :, ::a]

def init(): 
    for j in range(num_part):
        circles[j].center = (pos[0,j,0],pos[1,j,0])
        ax.add_patch(circles[j])
    
    ax.add_patch(rect)
    ax.add_patch(rect2)
    return circles

def animate(i): # amination-function. i = frame number. 
    for j in range(num_part):
        x1, y1 = pos[0,j, i], pos[ 1,j, i]
        circles[j].center = (x1, y1)
    return circles

circles = []
for i in range(num_part):
    circles.append(plt.Circle((0, 0), 0.2, color='r'))

rect = patches.Rectangle((0,-Ly), Lx, 2* Ly, edgecolor = "blue")
rect2 = patches.Rectangle((-Lx,-Ly), Lx, 2* Ly, edgecolor = "blue", fill = False)
plt.style.use('dark_background') 
#plt.style.use('dark_background') 
fig = plt.figure( figsize = (6, 6))
ax = plt.axes(xlim = (-12, 12), ylim = (-12 , 12))
line1, = ax.plot([], [], lw = 2, color = "c")

frames =  int(num_it / a) 
plt.axis("off") 

#animate
anim = animation.FuncAnimation(fig, animate, init_func= init, frames = frames, interval = 10, blit = True) #interval given in ms
plt.draw()
plt.show()