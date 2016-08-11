import os.path

import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# Plot array with matplotlib imshow()

''' A function to plot "array". The array is assumed to be a matrix.
'''
def plotArray(array):
    hf = plt.figure()
    ha = hf.add_subplot(111)
    plt.imshow(array, cmap = plt.get_cmap("viridis"))
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

    # Read in the rest of the data (binary)
    data = np.fromfile(f, np.float64, nParticles, "")

    # Compute the size of the rows and columns and transform the vector into an array
    dimension = np.sqrt(nParticles)

    # Check if the nParticles is a perfect square
    if (int(dimension) * int(dimension)!= nParticles):
        print("The read number of particles is not a perfect square so the vector cannot be reshaped into a square matrix.")
        return nParticles, 0
    else:
        # Reshape the data into a square matrix
        #data = np.reshape(data, (dimension, dimension))

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



''' A function that animates 2D data where each timesteps data is contained in a file in the folder
input:
t0:         initial frame                                       DEF = 0
t_end:      final frame                                         DEF = 4
v_min:      minimum value for color scale                       DEF = 0
v_max:      maximum value for color scale                       DEF = 4
folder:     location of files w.r.t. this script's directory    DEF = "Test_output_files/5t/"
colormap:   type of colormap to use                             DEF = plt.get_cmap("BrBG")
'''
def animateParticles(t0, t_end, writeFreq, folder, colormap = plt.get_cmap("viridis")):
    # Go to the script's directory
    curDir = os.getcwd()
    os.chdir(curDir)

    # Create array with time values (t's are inclusive)
    times = np.arange(t0+1, t_end + 1, writeFreq)

    # Prepare image and frame array
    fig = plt.figure()
    fig.add_subplot(111)
    frames = []

    # Read initial data
    file_0_Q = folder + "/" + "0_Q.txt"
    file_0_X = folder + "/" + "0_X.txt"
    file_0_Y = folder + "/" + "0_Y.txt"

    nParticlesQ, data_0_Q = readFile(file_0_Q)
    nParticlesX, data_0_X = readFile(file_0_X)
    nParticlesY, data_0_Y = readFile(file_0_Y)

    # Check if the number of particles in the files is equal
    if(nParticlesQ != nParticlesX or nParticlesX != nParticlesY or nParticlesY != nParticlesQ):
        print("The number of particles in the input files are not equal.")

    # Compute the extent of the initial grid and data values
    q_min = np.min(data_0_Q)
    q_max = np.max(data_0_Q)
    xmin = np.min(data_0_X)
    ymin = np.min(data_0_Y)
    xmax = np.max(data_0_X)
    ymax = np.max(data_0_Y)

    #Create a new grid
    grid_x, grid_y = np.mgrid[xmin:xmax:100j, ymin:ymax:200j]

    # Interpolate data onto grid
    coordinates = np.column_stack((data_0_X, data_0_Y))
    data_0_grid = griddata(coordinates, data_0_Q.squeeze(), (grid_x, grid_y), method='cubic')

    # Add frame to the list of frames
    frames.append([plt.imshow(data_0_grid.T, cmap=colormap, extent=(xmin,xmax,ymin,ymax), animated=True, vmin=10*q_min, vmax=10*q_max)])


    # Append the image list with each frame for the animation
    for t in times:
        print(t)
        # Define new filenames
        filenameQ = folder + "/" + str(t) + "_Q.txt"
        filenameX = folder + "/" + str(t) + "_X.txt"
        filenameY = folder + "/" + str(t) + "_Y.txt"

        # Read the files
        dataQ = readFile(filenameQ)[1]
        dataX = readFile(filenameX)[1]
        dataY = readFile(filenameY)[1]

        # Interpolate the data
        coordinates = np.column_stack((dataX, dataY))
        data_grid = griddata(coordinates, dataQ.squeeze(), (grid_x, grid_y), method='linear')

        frames.append([plt.imshow(data_grid.T, cmap = colormap, extent=(xmin,xmax,ymin,ymax), animated=True, vmin = 10*q_min, vmax = 10*q_max)])

    # Transform the list into an animation
    ani = anim.ArtistAnimation(fig, frames, interval=100, blit=False, repeat=False)
    plt.colorbar()
    plt.show()

    Writer = anim.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)

    ani.save('animation.mp4', writer=writer)

    print(q_min, q_max)

    return


def scatterPlot(folder):
    # Go to the script's directory
    curDir = os.getcwd()
    os.chdir(curDir)

    # Prepare image and frame array
    fig = plt.figure()
    fig.add_subplot(111)

    # Read data
    filenameX = folder + "/1_X.txt"
    filenameY = folder + "/1_Y.txt"
    nParticles, dataX = readFile(filenameX)
    nParticles, dataY = readFile(filenameY)

    plt.scatter(dataX, dataY, s=5)
    plt.show()

    return


##########################################
# SET THESE PARAMETERS FOR YOUR ANIMATION
# Files are assumed to be named:
# [ 0_X.txt, 0_Y.txt, 0_Q.txt, 1_X.txt,    ...,    t_end_X.txt, t_end_Y.txt, t_end_Q.txt]

# Animate from t_0 to t_end (inclusive)
t_0 = 0
t_end = 15
writeFreq = 1
foldername = "/home/shoshijak/Documents/ETH-FS16/HPC/p-shared"

# Make an animation
#animateParticles(t_0, t_end, writeFreq, foldername)

# Make a scatter pot of the x and y data
scatterPlot(foldername)
