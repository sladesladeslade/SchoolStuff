# Sliders for MPL Anims
# Slade Brooks
# brooksl@mail.uc.edu
# also stolen

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np


class sliders():
    
    def __init__(self):
        self.flag_init = True
        self.fig = plt.figure(1)
        self.axroll = plt.axes([0.25, 0.01, 0.65, 0.03])
        self.roll=0
        self.roll_slider = Slider(
            ax=self.axroll,
            label='roll',
            valmin=-np.pi,
            valmax=np.pi,
            valstep=0.1,
            valinit=0,
        )
        self.roll_slider.on_changed(self.update)
        self.axyaw = plt.axes([0.01, 0.25, 0.0225, 0.63])
        self.yaw_slider = Slider(
            ax=self.axyaw,
            label="yaw",
            valmin=-np.pi,
            valmax=np.pi,
            valinit=0,
            valstep=0.01,
            orientation="vertical"
        )
        self.yaw_slider.on_changed(self.update)

        self.axptich = plt.axes([0.05, 0.25, 0.0225, 0.63])# location of the slider on figure
        self.pitch_slider = Slider(
            ax=self.axptich,
            label="pitch",
            valmin=-np.pi/2,
            valmax=np.pi/2,
            valinit=0,
            valstep=0.1,
            orientation="vertical"
        )
        self.pitch_slider.on_changed(self.update)
        
        
    def update(self,val):
        self.roll=self.roll_slider.val
        self.yaw=self.yaw_slider.val
        self.pitch=self.pitch_slider.val
        #plt.pause(0.001)
        self.fig.canvas.draw_idle()