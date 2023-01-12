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
        ans = input("Would you like a contour plot (1), a line plot (2), or to end (3)? ")

        # if ans=1, generate contour plot
        if ans == "1":
            # ask for x and y ranges
            vmin = float(input("Input the minimum val of the range: "))
            vmax = float(input("Input the maximum val of the range: "))

            # create list of points between min and max vals
            xs = np.arange(vmin, vmax, 0.1)
            ys = np.arange(vmin, vmax, 0.1)
            fs = np.zeros((ys.size,xs.size))
            
            # create a grid of x and y values
            [X, Y] = np.meshgrid(xs, ys)

            # calculate f values
            fs = inputs[0]*np.sin(inputs[2]*np.pi*X + inputs[4])
            + inputs[1]*np.sin(inputs[3]*np.pi*Y + inputs[5])

            # plot contour
            plt.contourf(X, Y, fs, 20, cmap="RdGy")

            # set plot range
            plt.axis([vmin, vmax, vmin, vmax])

            # add plot formatting
            plt.colorbar()
            plt.title("Contour Plot")
            plt.xlabel("X Axis")
            plt.ylabel("Y Axis")

            # display plot
            plt.show()

        # if ans=2, generate line plot
        if ans == "2":
            # ask which direction line plot
            ans = input("Would you like a vertical (V) or horizontal (H) line plot? ")

            # if V, do this
            if ans == "V":
                # 
                fixedv = float(input("Input the value of x: "))
                vmin = float(input("Input the minimum value of y: "))
                vmax = float(input("Input the maximum value of y: "))

            # if H, do this
            if ans == "H":
                #
                fixedv = float(input("Input the value of y: "))
                vmin = float(input("Input the minimum value of x: "))
                vmax = float(input("Input the maximum value of x: "))

        # if ans=3, end loop and output filenames of plots
        stop = True
    

# testing
if __name__ == "__main__":

    testinputs = inputs(1, 3, 1, 0.5, (np.pi)/4, -(np.pi)/3)
    testplot = plot(testinputs)