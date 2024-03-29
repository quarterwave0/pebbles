from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np

data = np.loadtxt("pde_output.txt", delimiter=":")
dw = data.shape[1]
data = data[:, 0:dw]

maximum = np.max(data)
minimum = np.min(data)

fig, ax = plt.subplots()
line, = ax.plot(data[0])

ax.set(xlim=[0, dw-1], ylim=[minimum+minimum/10, maximum+maximum/10])

print(data.shape)

def update(frame):
    line.set_ydata(data[frame])
    return line


ani = animation.FuncAnimation(fig=fig, func=update, frames=data.shape[0], interval=1)

#ani.save("output.gif", writer="pillow")
plt.show()
