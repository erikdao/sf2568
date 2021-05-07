"""
Visualize edges
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    fname = None
    try:
        fname = sys.argv[1]
    except IndexError:
        sys.exit(-1)
    
    array = []
    with open(fname, "r") as f:
        for line in f:
            line_array = [int(n) for n in line.strip().split(" ")]
            array.append(line_array)
    
    data = np.array(array)
    print("data", data.shape, np.max(data), np.min(data))

    image = np.flip(data)

    fig, ax = plt.subplots()
    plt.axis('off')
    plt.imshow(image, cmap='gray')
    plt.savefig(sys.argv[2], bbox_inches='tight')

if __name__ == '__main__':
    main()
