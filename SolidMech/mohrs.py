# Slade Brooks
# brooksl@mail.uc.edu
# 02.14.2023

# Mohr's Circles Plotting Function


import matplotlib.pyplot as plt
import numpy as np


def mohrs(point1, point2, point3, type:str):
    y=0
    theta = np.linspace(0, 360, 10000)

    r1 = (np.abs(point1-point3))/2
    r2 = (np.abs(point1-point2))/2
    r3 = (np.abs(point2-point3))/2

    c1 = (point1+point3)/2
    c2 = (point1+point2)/2
    c3 = (point2+point3)/2

    mohr1x = r1*np.cos(theta) + c1
    mohr1y = r1*np.sin(theta)
    mohr2x = r2*np.cos(theta) + c2
    mohr2y = r2*np.sin(theta)
    mohr3x = r3*np.cos(theta) + c3
    mohr3y = r3*np.sin(theta)

    plt.scatter(mohr1x,mohr1y,color="black",s=0.5)
    plt.scatter(mohr2x,mohr2y,color="black",s=0.5)
    plt.scatter(mohr3x,mohr3y,color="black",s=0.5)

    if type == "stress":
        plt.plot(point3, y, "o", markersize=10,label="Sigma 3")
        plt.plot(point2, y, "o", markersize=10,label="Sigma 2")
        plt.plot(point1, y, "o", markersize=10,label="Sigma 1")

        plt.xlabel("Normal Stress (MPa)")
        plt.ylabel("Shear Stress (MPa)")
        plt.title("Mohr's Stress Circles")

    elif type == "strain":
        plt.plot(point3, y, "o", markersize=10,label="Epsilon 3")
        plt.plot(point2, y, "o", markersize=10,label="Epsilon 2")
        plt.plot(point1, y, "o", markersize=10,label="Epsilon 1")

        plt.xlabel("$\epsilon_{n}$")
        plt.ylabel("$(1/2)\gamma_{n}$")
        plt.title("Mohr's Strain Circles")
            
    plt.axis("scaled")
    plt.legend(loc="upper left")
    plt.show()


if __name__ == "__main__":
    mohrs(-.0004794, .0001345, .0006449, "strain")
    mohrs(-120.834, -6.6697, 53.5037, "stress")