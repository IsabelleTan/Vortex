import matplotlib.pyplot as plt
import numpy as np
import os.path
import matplotlib.animation as anim

# Plot array with matplotlib imshow()

''' A function to plot a vector. The data is assumed to be ordered lexicographically, in binary file format with double precision.
'''
def plotArray(array):
    hf = plt.figure()
    ha = hf.add_subplot(111)
    plt.imshow(array)
    plt.show()

    return

''' A function that reads the data from an binary file according to the format
xxx
'''
def readFile(filename):
    # Open the file
    f = open(filename, "rb")

    # Read and print the first decimal (!) number which represents the number of particles
    line = f.readline()
    nParticles = int(line)
    print("Number of particles: ", nParticles)

    # Read in the rest of the data
    data = np.fromfile(f, np.float64, nParticles, "")

    # Compute the size of the rows and columns and transform the vector into an array
    dimension = np.sqrt(nParticles)
    array = np.reshape(data, (dimension, dimension))

    # Return the number of particles and the data array
    return nParticles, data

''' A function to create a lamb oseen vortex with core radius CR'''
def LambOseen(CR):
    # Set x and y gridsize
    nx = 400
    ny = 400
    circ = 10
    cx = 200
    cy = 200


    data = np.zeros((nx,ny))

    for i in range(ny):
        for j in range(nx):
            rad = np.sqrt((i-cy)**2 + (j-cx)**2)
            data[i,j] = circ/(np.pi * CR**2) * np.exp(- (rad**2)/(CR**2))

    return data


##########################################

# Read in a vector from a txt file
# Go to the location of the script, depends on your personal system.
curDir = os.getcwd()
os.chdir(curDir)

# Prepare animation figure
fig = plt.figure()

# Set animation variables
t0 = 0
t_end = 100
dt = 1
times = np.arange(t0, t_end, dt)

# Read initial data
file_0 = "vortex_0"
nParticles, data_0 = readFile(file_0)

# Set initial image
im = plt.imshow(data_0, cmap=plt.get_cmap('viridis'), animated=True)

# Define the update function for the animation
def updatefig(*args):


    # Set array of image to the array specified by the filename
    im.set_array(filename)
    return im,

##########
# Try 2, now with for-loop to append image list with frames

# Prepare some stuff
CR_0=10
fig2 = plt.figure()
ims = []

# Append the image list with each frame for the animation
for t in times:
    # Define new filename
    filename = "vortex_" + str(t)
    # TODO Make the image use the same colormap
    ims.append([plt.imshow(readFile(filename), cmap=plt.get_cmap('viridis'), animated=True)])

# Transform the list into an animation
ani2 = anim.ArtistAnimation(fig2, ims, interval=50, blit=False, repeat=False)
plt.show()

