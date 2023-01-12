# Slade Brooks
# brooksl@mail.uc.edu
# 01.11.2023
# AeroComputing Assignment 1

# time for python

# This code plots line or contour plots of a function when given parameters.

# import libraries
import numpy as np
import matplotlib.pyplot as plt


# function to get and check inputs
def inputs(A, B, n, m, O, Y):

    # ask if inputs are correct
    # make flag for correct nums
    correct = False

    # set up loop to check if values are correct
    while correct == False:
        # display inputs
        print("You input: A =",A," B =",B," n =",n," m =",m," O =",O," Y =",Y)

        # query if inputs are right
        ans = input("Is that correct? (Y or N) ")

        # if ans=N, ask for inputs again
        if ans == "N":
            A = input("A =")
            B = input("B =")
            n = input("n =")
            m = input("m =")
            O = input("O =")
            Y = input("Y =")
        # if ans=Y, set flag to True
        elif ans == "Y":
            correct = True

    # return correct values
    return [A, B, n, m, O, Y]

# function to perform plotting
def plot(inputs):

    # set fxn stop flag
    stop = False

    # loop to ask for graph type or to stop
    while stop == False:
        
        # ask for type of plot or to end
        ans = input("Would you like a contour plot (1), a line plot (2), or to end (3)?")

        # if ans=1, generate contour plot
        if ans == "1":
            # ask for x and y ranges
            xmin = float(input("Input the minimum x val: "))
            xmax = float(input("Input the maximum x val: "))
            ymin = float(input("Input the minimum y val: "))
            ymax = float(input("Input the maximum y val: "))

            # initialize arrays
            xs = np.empty([100, 100])
            ys = np.empty([100, 100])
            fs = np.empty([100, 100])

            # find x and y values
            for i in range(100):
                xs[i] = np.linspace(xmin, xmax, num=100)
                ys[i] = np.linspace(ymin, ymax, num=100)

            # calculate f values
            for k, j in range(xs.size[0],xs.shape[1]):
                fs[k, j] = inputs[0]*np.sin(inputs[2]*np.pi*xs[k, j] + inputs[4]) + inputs[1]*np.sin(inputs[3]*np.pi*ys[k, j] + inputs[5])

            # plot contour
            plt.contourf(xs, ys, fs, 20, cmap="RdGy")
            plt.colorbar()
            # set plot range
            plt.axis([xmin, xmax, ymin, ymax])
            plt.title("Contour Plot")
            plt.xlabel("X Axis")
            plt.ylabel("Y Axis")

            plt.show()

        # if ans=2, generate line plot

        # if ans=3, end loop and output filenames of plots
        stop = True
    

# testing
if __name__ == "__main__":

    testinputs = inputs(1, 3, 1, 0.5, (np.pi)/4, (-np.pi)/3)
    testplot = plot(testinputs)