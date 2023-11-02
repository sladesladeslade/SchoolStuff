import numpy as np
import matplotlib.pyplot as plt
import control as ctl
from control.matlab import *
import lib.mssimparams as P


# tf num and denom
kD = 7.2
kP = 3.05
num = [1]
den = [P.m1, (P.b + kD), (P.k + kP), 0.]

# run tf
sys=ctl.tf(num, den)
print(sys)
ctl.rlocus(sys)
plt.show()