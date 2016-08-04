import os.path

import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np

# Plot array with matplotlib imshow()

''' A function to plot "array". The array is assumed to be a matrix.
'''
def plotArray(array):
    hf = plt.figure()
    ha = hf.add_subplot(111)
    plt.imshow(array)
    plt.show()

    return

''' A function that reads a vector from a binary file according to the format

int (in decimal, number of following floats) \n
float/double (in binary) float/double (in binary) ... float/double (in binary)

The function returns the read number of particles as an integer and the data as a square array.
If the number of particles is not a perfect square the function prints an error message and exits.
'''
def readFile(filename):
    # Go to the location of the script, depends on your personal system.
    curDir = os.getcwd()
    os.chdir(curDir)

    # Open the file
    f = open(filename, "rb")

    # Read and print the first decimal (!) number which represents the number of particles
    line = f.readline()
    nParticles = int(line)
    print("Number of particles: ", nParticles)

    # Read in the rest of the data (binary)
    data = np.fromfile(f, np.float32, nParticles, "")

    # Compute the size of the rows and columns and transform the vector into an array
    dimension = np.sqrt(nParticles)

    # Check if the nParticles is a perfect square
    if (int(dimension) * int(dimension)!= nParticles):
        print("The read number of particles is not a perfect square so the vector cannot be reshaped into a square matrix.")
        return nParticles, 0
    else:
        # Reshape the data into a square matrix
        array = np.reshape(data, (dimension, dimension))

        # Return the number of particles and the data array
        return nParticles, array

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

''' A function that animates 2D data where each timesteps data is contained in a file in the folder
input:
t0:         initial frame                                       DEF = 0
t_end:      final frame                                         DEF = 4
v_min:      minimum value for color scale                       DEF = 0
v_max:      maximum value for color scale                       DEF = 4
folder:     location of files w.r.t. this script's directory    DEF = "Test_output_files/5t/"
colormap:   type of colormap to use                             DEF = plt.get_cmap("BrBG")
'''
def animate(t0 = 0, t_end = 4, v_min = 0, v_max = 40, folder = "Test_output_files/5t/", colormap = plt.get_cmap("BrBG")):
    # Go to the script's directory
    curDir = os.getcwd()
    os.chdir(curDir)

    # Create array with time values (t's are inclusive)
    times = np.arange(t0+1, t_end + 1, 1)

    # Prepare image and frame array
    fig = plt.figure()
    fig.add_subplot(111)
    frames = []

    # Read initial data
    file_0 = folder + "test_0"
    nParticles, data_0 = readFile(file_0)
    frames.append([plt.imshow(data_0, cmap = colormap, animated=True, vmin = v_min, vmax = v_max)])

    # Append the image list with each frame for the animation
    for t in times:
        # Define new filename
        filename = folder + "test_" + str(t)
        frames.append([plt.imshow(readFile(filename)[1], cmap = colormap, animated=True, vmin = v_min, vmax = v_max)])

    # Transform the list into an animation
    ani = anim.ArtistAnimation(fig, frames, interval=500, blit=False, repeat=False)
    plt.colorbar()
    plt.show()

    Writer = anim.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    ani.save('test.mp4', writer=writer)

    return

##########################################
animate()
