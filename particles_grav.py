import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def init_pos(radius,numberOfParticles):
    to_close = True
    itr = 0
    while(to_close):
        to_close = False
        r = np.random.uniform(0,radius, numberOfParticles)
        itr += 1
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
    return pos_x, pos_y, itr

def random_velocities(numberOfParticles):
    angles =  np.random.uniform(-np.pi, np.pi,numberOfParticles)
    vx = np.cos(angles) / np.sqrt(numberOfParticles / 2)  /4
    vy = np.sin(angles) / np.sqrt(numberOfParticles / 2)  /4
    return vx, vy

def vel_one_particle(numberOfParticles):
    vx = np.zeros(numberOfParticles)
    vy = np.zeros(numberOfParticles)
    angle = np.random.uniform(-np.pi,np.pi, 1)
    vx[0] = np.cos(angle) * np.sqrt(2)
    vy[0] = np.sin(angle) * np.sqrt(2)
    return vx, vy

def gravity(numberOfParticles, numberOfIterations, radius, dt, KK, vel_dist = random_velocities):
    pos = np.zeros((2, numberOfParticles, numberOfIterations+1))
    pos[0,:,0], pos[1,:,0], itr = init_pos(radius,numberOfParticles)
    print("Placing the particles non-overlapping was a success. It took {} attempts.".format(itr))
    vx, vy = vel_dist(numberOfParticles)
    vx_list = np.zeros((numberOfIterations,numberOfParticles,))
    vy_list = np.zeros((numberOfIterations,numberOfParticles,))
    
    
    for i in range(numberOfIterations):
        vx_list[i] = vx
        vy_list[i] = vy
        distance = np.sqrt(pos[0,:,i]**2 + pos[1,:,i]**2)
        acceleration = np.zeros((2, numberOfParticles))
        
        #test if we are outside the circle and calculate force from wall
        acceleration[0, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[0, distance > radius, i] / distance[distance > radius]
        acceleration[1, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[1, distance > radius, i] / distance[distance > radius]
        for j in range(numberOfParticles):
                for k in range(j+1, numberOfParticles):
                    F = 0
                    d_jk = np.sqrt((pos[0,j,i] - pos[0,k,i])**2 + (pos[1,j,i] - pos[1,k,i])**2)
                    if d_jk < 1:
                        F += 12 * epsilon * (1/d_jk**7 - 1/d_jk**13)
                    F +=  epsilon * (1/d_jk**2) #Here: postive force = attraction
                    
                    
                    acceleration[0,j] += - F * (pos[0,j,i] - pos[0,k,i]) / d_jk 
                    acceleration[1,j] += - F * (pos[1,j,i] - pos[1,k,i]) / d_jk 
                    acceleration[0,k] += F * (pos[0,j,i] - pos[0,k,i]) / d_jk 
                    acceleration[1,k] += F * (pos[1,j,i] - pos[1,k,i]) / d_jk
                    
        #update positions
        pos[0, :, i+1] =  pos[0, :, i] + vx * dt + (acceleration[0, :] * dt**2)/2
        pos[1, :, i+1] =  pos[1, :, i] + vy * dt + (acceleration[1, :] * dt**2)/2
        
        distance = np.sqrt(pos[0,:,i+1]**2 + pos[1,:,i+1]**2)
        acceleration2 = np.zeros((2, numberOfParticles))

        #test if we are outside the circle and calculate force from wall
        acceleration2[0, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[0, distance > radius, i+1] / distance[distance > radius]
        acceleration2[1, distance > radius] = -KK * (distance[distance > radius] - radius) * pos[1, distance > radius, i+1] / distance[distance > radius]
        for j in range(numberOfParticles):
                for k in range(j+1, numberOfParticles):
                    F = 0
                    d_jk = np.sqrt((pos[0,j,i+1] - pos[0,k,i+1])**2 + (pos[1,j,i+1] - pos[1,k,i+1])**2)
                    if d_jk < 1:
                        F += 12 * epsilon * (1/d_jk**7 - 1/d_jk**13)
                    F +=  epsilon * (1/d_jk**2) #Here: postive force = attraction
                    
                    acceleration2[0,j] += - F * (pos[0,j,i+1] - pos[0,k,i+1]) / d_jk 
                    acceleration2[1,j] += - F * (pos[1,j,i+1] - pos[1,k,i+1]) / d_jk 
                    acceleration2[0,k] += F * (pos[0,j,i+1] - pos[0,k,i+1]) / d_jk 
                    acceleration2[1,k] += F * (pos[1,j,i+1] - pos[1,k,i+1]) / d_jk
                        
        #update velocities
        vx += (acceleration[0,:] + acceleration2[0,:])/2 * dt
        vy += (acceleration[1,:] + acceleration2[1,:])/2 * dt
    

    return pos

num_part = 10
num_it = 10000
r = 10
dt = 0.01
KK = 0.01 #hardness of wall
epsilon = 0.1 #strength of gravity and LJ
pos = gravity(num_part, num_it, r, dt, KK )
a = 10 #plot every a'th value
pos = pos[:, :, ::a]

def init(): 
    for j in range(num_part):
        circles[j].center = (pos[0,j,0],pos[1,j,0])
        ax.add_patch(circles[j])
    
    if not KK  < 0.3:
        ang = np.linspace(0, 2 * np.pi, 100)
        line1.set_data(np.cos(ang) * r, np.sin(ang) * r)
    return circles

def animate(i): # amination-function. i = frame number. 
    for j in range(num_part):
        x1, y1 = pos[0,j, i], pos[ 1,j, i]
        circles[j].center = (x1, y1)
    return circles

circles = []
for i in range(num_part):
    circles.append(plt.Circle((0, 0), 0.5, color='r'))

plt.style.use('dark_background') 
#plt.style.use('dark_background') 
fig = plt.figure( figsize = (6, 6))
ax = plt.axes(xlim = (-1.5 * r, 1.5*r), ylim = (-1.5*r , 1.5*r))
line1, = ax.plot([], [], lw = 2, color = "c")

frames =  int(num_it / a) 
plt.axis("off") 

#animate
anim = animation.FuncAnimation(fig, animate, init_func= init, frames = frames, interval = 10, blit = True) #interval given in ms
plt.draw()
plt.show()