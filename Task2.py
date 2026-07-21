# modelling the brownian motion of a small but visible particle using the zero-momentum frame
import matplotlib.pyplot as plt
import numpy as np
import random
import matplotlib.animation as animation
from matplotlib.patches import Circle
from math import cos
from math import sin
from math import pi
from math import sqrt


# variables: time in ps and length in nm

N=250
R=10
r=R/10
m=5*10**-26
M=50*m
sim_length = 15*R
Kn=15
C=0.8
mean_free_path = r*Kn
small_speed= 1000 # converting m/s into nm/ps
dir_change_timestep = mean_free_path/small_speed
dt=0.001
steps_per_frame=2

physics_time=0

def dot(vec1,vec2):
    return vec1[0]*vec2[0] + vec1[1]*vec2[1]

def distance(x,y):
    return sqrt(x**2 + y**2)




# initialize all small particles
small_positions = np.array([[random.uniform(r,sim_length-r),random.uniform(r,sim_length-r)] for _ in range(N)])
# ensure none spawn touching the large circle
for s in range(len(small_positions)):
    modify=small_positions[s]
    while distance(modify[0]-sim_length/2,modify[1]-sim_length/2)<r+R:
        small_positions[s] = [random.uniform(r,sim_length-r),random.uniform(r,sim_length-r)]

small_dirs = np.array([2*pi*random.random() for _ in range(N)])
small_vels = np.array([[small_speed*cos(small_dirs[i]),small_speed*sin(small_dirs[i])]for i in range(N)])
small_next_dir_changes = np.array([random.uniform(0,dir_change_timestep) for _ in range(N)])

large_position = np.array([sim_length/2,sim_length/2])
large_vel = np.array([0,0])

fig,ax = plt.subplots()

ax.set_xlim(0,sim_length)
ax.set_ylim(0,sim_length)
ax.set_aspect("equal")
ax.set_xlabel("x / nm")
ax.set_ylabel("y / nm")
ax.set_title("Brownian motion simulation")

small_scatter = ax.scatter(
    small_positions[:,0],
    small_positions[:,1],
    s=5
)

large_circle = Circle(
    large_position,
    R,
    fill=True,
    color='red',
    linewidth=2
)
ax.add_patch(large_circle)

trail_x = [large_position[0]]
trail_y = [large_position[1]]

trail_line, = ax.plot(trail_x,trail_y,linewidth=1)

start_marker, = ax.plot(
    large_position[0],
    large_position[1],
    marker="*",
    markersize=10,
    linestyle="None"
)
curr_marker, = ax.plot(
    large_position[0],
    large_position[1],
    marker="o",
    markersize=6,
    linestyle="None"
)

def physics_step():
    # one physics update run
    global physics_time
    global large_position
    global large_vel

    for x in range(N):
        if small_next_dir_changes[x]<physics_time:
            # choose and apply a new random direction
            small_dirs[x] = 2*pi*random.random()
            small_vels[x]= [small_speed*cos(small_dirs[x]),small_speed*sin(small_dirs[x])]
            small_next_dir_changes[x]+=random.uniform(0,dir_change_timestep)

        small_pos = small_positions[x]
        small_vel = small_vels[x]

        separation_vector = small_pos-large_position
        magnitude=distance(separation_vector[0],separation_vector[1])
        magnitude = 0.001 if magnitude==0 else magnitude
        relative_vel=small_vel-large_vel
        normal = separation_vector/magnitude



        # check against conditions for handling a collision
        if magnitude < R+r:
            if dot(relative_vel,normal)<0:
                # handle collision using the ZMF
                small_normal_vel = dot(normal,small_vel)*normal
                small_tangent_vel = small_vel-small_normal_vel

                large_normal_vel = dot(normal,large_vel)*normal
                large_tangent_vel = large_vel-large_normal_vel

                V = (m*small_normal_vel + M*large_normal_vel)/(M+m)

                small_vels[x] = C*(V-small_normal_vel)+V + small_tangent_vel
                large_vel = C*(V-large_normal_vel)+V + large_tangent_vel

                small_positions[x] = large_position + normal*(R+r)
                # make sure small particle no longer inside large one
        
        small_positions[x] += small_vels[x]*dt
    
    large_position += large_vel*dt

    physics_time+=dt

    trail_x.append(large_position[0])
    trail_y.append(large_position[1])

def update_frame(frame):
    for _ in range(steps_per_frame):
        physics_step()

    small_scatter.set_offsets(small_positions)
    large_circle.center = large_position
    trail_line.set_data(trail_x,trail_y)

    curr_marker.set_data(
        [large_position[0]],
        [large_position[1]]
    )

ani = animation.FuncAnimation(fig=fig,func=update_frame,frames=100,interval=300)
plt.show()

    

    














