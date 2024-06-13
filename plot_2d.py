import math

from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np

f = open("pde_output_2d.txt", "r")
d = f.readlines()

timespan = len(d)
dataSet = np.zeros((timespan, 102, 102), dtype=np.float32)

for t in range(timespan):
    s1 = d[t].split(":")
    for x in range(len(s1)):
        s2 = np.fromstring(str(s1[x][1:-1]), sep=",", dtype=np.float32)
        dataSet[t, x, :] = s2

dw = dataSet.shape[1]
maximum = np.max(dataSet)
minimum = np.min(dataSet)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

x = np.arange(dw)
y = np.arange(dw)

(X, Y) = np.meshgrid(x, y)
line = ax.plot_surface(X, Y, dataSet[0, :, :], cmap='viridis', antialiased=True)

ax.set_zlim(minimum+minimum/10, maximum+maximum/10)

def update(frame):
    ax.clear()
    ax.set_zlim(minimum+minimum/10, maximum+maximum/10)
    return ax.plot_surface(X, Y, dataSet[frame, :, :], cmap='viridis', antialiased=True)

ani = animation.FuncAnimation(fig=fig, func=update, frames=timespan, interval=5)
ani.save("output2.gif", writer="pillow")
plt.show()