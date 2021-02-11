import numpy as np
import matplotlib.pyplot as plt
import sys


file_name = sys.argv[1]
F = open(file_name,"r")
A = []
for line in F:
    lin = line.split(" ")[:-1]
    lin = [int(i) for i in lin]
    A.append(lin)

A = np.array(A)
print(np.max(np.max(A)))
#plt.figure()
#plt.imshow(A.T,cmap='RdBu')
plt.imsave('Im2048_511_zoom_redgy.png',A.T,cmap='RdGy')
#plt.show()
