from matplotlib import pyplot as plt
import numpy as np

file_name = "possion.txt"

F = open(file_name,"r")

for line in F:
    lin = line.split(" ")[:-1]
    lin = [float(i) for i in lin]

#manufactured solution:

u = lambda x: x*(x-1)

x = np.linspace(0,1,1000)
y = u(x)

fig, ax = plt.subplots()
approx, = plt.plot(x,lin, label="Approximated u")
true, = plt.plot(x,y, label="True u")
plt.xlabel("x")


plt.ylabel("u(x)")
ax.set_xlim(xmin=0,xmax=1)
plt.legend(handles=[approx, true])
fig.savefig('u_plot_8_processors_10e6.png')
plt.show()