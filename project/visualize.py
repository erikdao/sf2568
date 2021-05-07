"""
Visualize edges
"""
import sys
import numpy as np
from PIL import Image

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

    image = Image.fromarray(np.flip(data, axis=0).astype(np.uint8))
    image.save(sys.argv[2])

if __name__ == '__main__':
    main()
