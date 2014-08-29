import numpy as np

def split(filename):

    data = np.loadtxt(filename, skiprows=1)
    
    # Put on Angstrom scale.
    data[:,0] *= 10.

    output_filename_base = ".".join(filename.split(".")[:-1])

    # Split at points:
    split_points = (5500, 6250, 7000)
    previous_point = 0
    for i, point in enumerate(split_points):
        index = data[:,0].searchsorted(point)
        
        this_channel = data[previous_point:index]
        np.savetxt("{0}_{1}.txt".format(output_filename_base, i + 1), this_channel[::10])
        previous_point = index
    last_channel = data[previous_point:]
    np.savetxt("{0}_{1}.txt".format(output_filename_base, i + 2), last_channel[::10])

    print("done")

