from glob import glob
import numpy as np

def split(filename):

    data = np.loadtxt(filename, skiprows=1)
    
    # Put on Angstrom scale.
    data[:,0] *= 10.

    output_filename_base = ".".join(filename.split(".")[:-1])

    # Split at points:
    split_points = (5500, 6250, 7000)
    color = ("blue", "green", "red", "ir")
    previous_point = 0
    for i, point in enumerate(split_points):
        index = data[:,0].searchsorted(point)
        
        this_channel = data[previous_point:index]
        sampling_rate = int(np.floor(len(this_channel)/2000.))
        np.savetxt("{0}_{1}.txt".format(output_filename_base, color[i]), this_channel[::sampling_rate])
        previous_point = index
    last_channel = data[previous_point:]
    sampling_rate = int(np.floor(len(last_channel)/2000.))
    np.savetxt("{0}_{1}.txt".format(output_filename_base, color[-1]), last_channel[::sampling_rate])

    print("done")


if __name__ == "__main__":
    files = glob("*narval.txt")
    map(split, files)
