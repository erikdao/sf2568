from matplotlib import pyplot as plt
import numpy as np
import math

file_name = "possion.txt"

F = open(file_name,"r")

for line in F:
    lin = line.split(" ")[:-1]
    lin = [float(i) for i in lin]

#manufactured solution:

u = lambda x: x*(x-0.5)*(x-1)
x = np.linspace(0,1,1000)
y = u(x)

fig, ax = plt.subplots()
approx, = plt.plot(x,lin, label="Approximated u", linewidth=2., alpha=0.8)
true, = plt.plot(x,y, label="True u", linewidth=2., alpha=0.8)
plt.xlabel("x", fontsize=14)
plt.ylabel(r"$u(x)$", fontsize=14)
ax.set_xlim(xmin=0,xmax=1)
plt.legend(handles=[approx, true])
fig.savefig('u_plot_8_processors_10e6_3.png')
plt.show()