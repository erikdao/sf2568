"""
Visualize edges
"""
import os
import sys
import numpy as np
# import png
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

    fig, ax = plt.subplots()
    plt.imshow(data)
    plt.savefig(f"{fname}.png")
    # png.from_array(data, 'L').save(f"{fname}.png")

if __name__ == '__main__':
    main()
