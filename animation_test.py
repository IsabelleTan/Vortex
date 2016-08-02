import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

#fig = plt.figure()

''' A function to create a lamb oseen vortex with core radius CR'''
def LambOseen(CR):
    # Set x and y gridsize
    nx = 100
    ny = 100
    circ = 10
    cx = 50
    cy = 50

    data = np.zeros((nx,ny))

    for i in range(ny):
        for j in range(nx):
            rad = np.sqrt((i-cy)**2 + (j-cx)**2)
            data[i,j] = circ/(np.pi * CR**2) * np.exp(- (rad**2)/(CR**2))

    return data

# Prepare some stuff
CR_0=10
fig2 = plt.figure()
ims = []

# Append the image list with each frame for the animation
for i in range(1,10):
    # TODO Make the image use the same colormap
    ims.append([plt.imshow(LambOseen(CR_0+i*10), cmap=plt.get_cmap('viridis'), animated=True)])

# Transform the list into an animation
ani2 = anim.ArtistAnimation(fig2, ims, interval=50, blit=False, repeat=False)
plt.show()

# Save the animation
# TODO Need to download FFMPEG!
'''
Writer = anim.writers['ffmpeg']
writer = Writer(fps=5, metadata=dict(artist = 'Me'), bitrate=1800)

ani2.save('test.mp4', writer=writer)'''
