import numpy as np
import matplotlib.pyplot as plt
import sys


file_name = sys.argv[1]
F = open(file_name,"r")
A = []
el_count = 2048
for line in F:
    lin = line.split(" ")[:-1]
    if len(lin) != el_count:
        break
    lin = [int(i) for i in lin]
    A.append(lin)

A = np.array(A)
print(np.max(np.max(A)))
plt.figure()
# # plt.imshow(A.T,cmap='RdGy')
plt.imshow(A.T, cmap='RdGy')
# # plt.imsave('Im2048_511_zoom_redgy.png',A.T,cmap='RdGy')
plt.show()
