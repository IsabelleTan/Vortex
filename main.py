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

''' A function to interpolate data onto a grid'''
def interpolate(dataX, dataY, dataQ):
    # check if the number of particles in the three files is equal
    if(len(dataQ) != len(dataX) or len(dataX) != len(dataY) or len(dataY) != len(dataQ)):
        print("The number of particles in the input files are not equal.")


    # Compute the extent of the coordinates
    xmin = np.min(dataX)
    ymin = np.min(dataY)
    xmax = np.max(dataX)
    ymax = np.max(dataY)


    #Create a new grid
    grid_x, grid_y = np.mgrid[xmin:xmax:200j, ymin:ymax:200j]

    # Interpolate data onto grid
    coordinates = np.column_stack((dataX, dataY))
    dataGrid = griddata(coordinates, dataQ.squeeze(), (grid_x, grid_y), method='cubic')

    return dataGrid

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
    times = np.arange(t_0+writeFreq, t_end + 1, writeFreq)

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


    # Compute the minimum and maximum q, x and y values
    q_min = np.min(data_0_Q)
    q_max = np.max(data_0_Q)
    xmin = np.min(data_0_X)
    ymin = np.min(data_0_Y)
    xmax = np.max(data_0_X)
    ymax = np.max(data_0_Y)

    print(q_min, q_max)

    # Interpolate the data
    data_0_grid = interpolate(data_0_X, data_0_Y, data_0_Q)

    # Add frame to the list of frames
    frames.append([plt.imshow(data_0_grid.T, cmap=colormap, extent=(xmin,xmax,ymin,ymax), animated=True, vmin=q_min, vmax=q_max)])


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
        data_grid = interpolate(dataX, dataY, dataQ)

        frames.append([plt.imshow(data_grid.T, cmap = colormap, extent=(xmin,xmax,ymin,ymax), animated=True, vmin = q_min, vmax =q_max)])

    # Transform the list into an animation
    ani = anim.ArtistAnimation(fig, frames, interval=100, blit=False, repeat=False)
    plt.colorbar()
    plt.show()

    Writer = anim.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)

    ani.save('animation.mp4', writer=writer)

    print(q_min, q_max)

    return


def scatterPlot(folder, t_0, t_end, writeFreq):
    # Go to the script's directory
    curDir = os.getcwd()
    os.chdir(curDir)

    # Prepare image and frame array
    fig = plt.figure()
    fig.add_subplot(111)

    # Loop over timesteps
    times = np.arange(t_0, t_end + 1, writeFreq)

    for t in times:
        # Create filename string
#        filenameX = folder + "/" + str(t) + "_X.txt"
#        filenameY = folder + "/" + str(t) + "_Y.txt"
        nParticles, dataX = readFile(filenameX)
        nParticles, dataY = readFile(filenameY)
        plt.scatter(dataX, dataY, s=5)
        title = "(x,y): step " + str(t)
        plt.title(title)
        plt.show()

    return

def valuesPlot(folder, t_0, t_end, writeFreq, colormap = plt.get_cmap("viridis")):
    # Go to the script's directory
    curDir = os.getcwd()
    os.chdir(curDir)

    # Prepare image and frame array
    fig = plt.figure()
    fig.add_subplot(111)
    fig2 = plt.figure()
    fig2.add_subplot(111)

    times = np.arange(t_0, t_end + 1, 1)

    for t in times:
        # Read data
        filenameX = folder + "/" + str(t) + "_X.txt"
        filenameY = folder + "/" + str(t) + "_Y.txt"
        filenameQ = folder + "/" + str(t) + "_Q.txt"
        filenameV = folder + "/" + str(t) + "_V.txt"
        filenameP = folder + "/" + str(t) + "_P.txt"
        nParticles, dataX = readFile(filenameX)
        nParticles, dataY = readFile(filenameY)
        nParticles, dataQ = readFile(filenameQ)
        nParticles, dataV = readFile(filenameV)
        nParticles, dataP = readFile(filenameP)

        # Compute the minimum and maximum q, x and y values
        v_min = np.min(dataV)
        v_max = np.max(dataV)
        q_min = np.min(dataQ)
        q_max = np.max(dataQ)
        p_min = np.min(dataP)
        p_max = np.max(dataP)

        xmin = np.min(dataX)
        ymin = np.min(dataY)
        xmax = np.max(dataX)
        ymax = np.max(dataY)

        print("(x,y) min/max values : ", xmin, xmax, ymin, ymax)

        # Interpolate the data
        data_gridV = interpolate(dataX, dataY, dataV)
        data_gridQ = interpolate(dataX, dataY, dataQ)
        data_gridP = interpolate(dataX, dataY, dataP)

        plt.imshow(data_gridV.T, cmap=colormap, extent=(xmin,xmax,ymin,ymax), vmin = v_min, vmax = v_max)
        plt.title("velocity at step " + str(t))
        plt.colorbar()
        plt.show()
        
        print("velocity :    ", np.min(dataV), np.max(dataV))

        plt.imshow(data_gridQ.T, cmap=colormap, extent=(xmin,xmax,ymin,ymax), vmin = q_min, vmax = q_max)
        plt.title("vorticity at step " + str(t))
        plt.colorbar()
        plt.show()
        
        print("vorticity : ", np.min(dataQ), np.max(dataQ))

        plt.imshow(data_gridP.T, cmap=colormap, extent=(xmin,xmax,ymin,ymax), vmin = p_min, vmax = p_max)
        plt.title("potential at step " + str(t))
        plt.colorbar()
        plt.show()
        
        print("potential : ", np.min(dataP), np.max(dataP))

    return

##########################################
# SET THESE PARAMETERS FOR YOUR ANIMATION
# Files are assumed to be named:
# [ 0_X.txt, 0_Y.txt, 0_Q.txt, 1_X.txt,    ...,    t_end_X.txt, t_end_Y.txt, t_end_Q.txt]

# Animate from t_0 to t_end (inclusive)
t_0 = 0
t_end = 99
writeFreq = 1
#foldername = "/home/shoshijak/Documents/ETH-FS16/HPC/p-shared"
foldername = "/Users/Isabelle/Documents/Studie/Master/Vakken/SS16/HPCSE2/Vortex/Test_output_files/simulation"

# Make an animation
#animateParticles(t_0, t_end, writeFreq, foldername)

# Make a scatter pot of the x and y data
#scatterPlot(foldername, 0, 6, 1)

# Plot the velocity
#velocityPlot(foldername)

# Plot the voriticty
#vorticityPlot(foldername)

#valuesPlot(foldername, t_0, t_end, writeFreq)

nParticles, data = readFile("/Users/Isabelle/Documents/Studie/Master/Vakken/SS16/HPCSE2/Vortex/cut_off.txt")
plt.plot(data)
plt.show()
