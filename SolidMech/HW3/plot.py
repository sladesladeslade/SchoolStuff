import matplotlib.pyplot as plt
import math
import numpy as np

theta = np.linspace(0, 360, 400)

P = 22.5
w = .01
h = .015

sigman = (P*(np.cos(theta*math.pi/180))**2)/(w*h)

taun = np.array([[(P*np.cos(theta*math.pi/180) - P*(np.cos(theta*math.pi/180))**3)/(w*h)],
                [0],
                [(-P*((np.cos(theta*math.pi/180))**2)*np.sin(theta*math.pi/180))/(w*h)]])

plt.plot(theta, sigman, label="Normal Stress", linewidth=2.5)
plt.plot(theta, np.linalg.norm(taun), label="Shear Stress", linewidth=2.5)
plt.xlabel("Theta (degrees)")
plt.ylabel("Stress (kPa)")
plt.title("Normal and Shear Stress")
plt.legend(loc="upper left")
plt.xlim(0, 360)
plt.ylim(0)
plt.show()