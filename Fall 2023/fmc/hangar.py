# Hangar for Anim Aircraft
# Slade Brooks
# brooksl@mail.uc.edu
# definitely gonna spend too much time on this

import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

"""
Vehicle List:
-------------
- sample UAV
- vision jet (Cirrus SF50)
"""

# sample UAV
fuse_l1 = 2.8
fuse_l2 = 0.1
fuse_l3 = 5
fuse_h = 1
fuse_w = 1
wing_l = 1.5
wing_w = 7
vtail_h = 1.5
vtail_l = 1
htail_l = 0.75
htail_w = 3.3
sampleUAV_verts = np.array([[fuse_l1, 0, 0],
                            [fuse_l2, fuse_w/2, fuse_h/2],
                            [fuse_l2, -fuse_w/2, fuse_h/2],
                            [fuse_l2, -fuse_w/2, -fuse_h/2],
                            [fuse_l2, fuse_w/2, -fuse_h/2],
                            [-fuse_l3, 0, 0],
                            [0, wing_w/2, 0],
                            [-wing_l, wing_w/2, 0],
                            [-wing_l, -wing_w/2, 0],
                            [0, -wing_w/2, 0],
                            [-(fuse_l3 - htail_l), htail_w/2, 0],
                            [-fuse_l3, htail_w/2, 0],
                            [-fuse_l3, -htail_w/2, 0],
                            [-(fuse_l3 - htail_l), -htail_w/2, 0],
                            [-(fuse_l3 - vtail_l), 0, 0],
                            [-fuse_l3, 0, vtail_h]
                            ])*[1, 1, -1]
def sampleUAV_obj(v):
    return np.array([[v[0], v[1], v[2]],
                     [v[0], v[1], v[4]],
                     [v[0], v[2], v[3]],
                     [v[1], v[3], v[4]],
                     [v[2], v[3], v[5]],
                     [v[1], v[2], v[5]],
                     [v[1], v[4], v[5]],
                     [v[3], v[4], v[5]],
                     [v[6], v[7], v[8]],
                     [v[6], v[8], v[9]],
                     [v[10], v[11], v[12]],
                     [v[10], v[12], v[13]],
                     [v[5], v[14], v[15]]
                     ])


# Vision Jet
visionJet_mesh = mesh.Mesh.from_file("C:\\Users\\spbro\\SchoolStuff\\Fall 2023\\fmc\\stls\\visionjet.stl")
visionJet_verts = visionJet_mesh.vectors*np.array([-1, 1, -1])
def visionJet_obj(): return