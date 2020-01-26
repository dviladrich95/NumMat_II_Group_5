import time

import numpy as np
from scipy import optimize
import scipy.sparse as sp
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation


### PLOT PARAMTER
### SET THIS TO DETERMINE THE NUMER
### OF TIMES SKIPPED IN THE ANIMATION
skip = 10
###
###

# end time
t_end =0.5



# paramter
a = 0.25
alpha = 0
beta = 0
L = 10

N = 100
h = 1.0/(N+1)
x = np.linspace(0,L, N)
l = 0.00015
g = a*l/h**2


# data
u0 = np.sin(np.pi*x/L) + 2*np.sin(3*np.pi*x/L) + np.sin(5*np.pi*x/L)
    
U = np.load("heatsol.npy")



#### Plot

fr = list(range(0, U.shape[1]))
fr = [skip*t for t in range(0, U.shape[1])]
T = np.arange(0,t_end, l)


fig, ax = plt.subplots()
ln, = plt.plot(x, u0)

def init():
    return ln,

def update(idx):
    ln.set_data(x, U[:,idx])
    ax.set_title("t="+ str(T[idx]))
    return ln,

ani = FuncAnimation(fig, update, frames=fr,
                    init_func=init, blit=False)
plt.show()
















