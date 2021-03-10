"""
Validate if an array is sorted.
The array is read from a text file
"""
import sys
import numpy as np


def main():
    arr = []

    with open(sys.argv[1], "r") as f:
        for line in f:
            arr.append(float(line.strip()))
    
    print("Array read! %d elements" % len(arr))

    ori_arr = np.array(arr)
    sorted_arr = np.sort(ori_arr)

    if (ori_arr == sorted_arr).all():
        print("Array is properly sorted!")
    else:
        print("Array might have not been sorted properly!")


if __name__ == '__main__':
    main()

