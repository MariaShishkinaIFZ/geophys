#!/usr/bin/env python3

import numpy as np
import argparse

# Stores the spatial data in predefined format:
# ZYX, meaning that data[k][j][i] returns a point
# with coordinates (x == i, y == j, z == k).
class VolumeData:

    # self.data contains the ZYX fload target data
    def __init__(self, file_in, format_in):

        self.data = np.zeros()


        a = array(floats,'float32')
        output_file = open('file', 'wb')
        a.tofile(output_file)
        output_file.close()


    # Reads the given binary file, which contents are interpreted as
    # according to the parameters into internal float ZYX numpy array.
    #
    # @file_path the path to the file to open
    # @iformat tuple of incoming data axis order:
    #   the order of axis scan (last axis increased
    #   with every next point, second with every next
    #   stride, first with every next layer)
    #   Example: ('x', 'z', 'y') means that first point
    #       in array has coords (ZYX order): (0; 0; 0)
    #       second (0; 1; 0), third (0; 2; 0), ....
    # @zyx_shape a tuple of 3 sizes (for Z, Y, X axis correspondingly)
    #       NOTE: must fit to the actual binary data size
    #   Example: (30, 60, 120)
    #       means 30 z=const layers, 60 y=const strides,
    #       120 x=const points in every stride in incoming data
    def read_binary(self, file_path, iformat, zyx_shape):

        if len(iformat) != 3: raise ValueError("iformat must be 3 items in size")
        if len(zyx_shape) != 3: raise ValueError("zyx_shape must be 3 items in size")

        read_shape = []
        zyx_shape_map = {'x' : 2, 'y' : 1, 'z' : 0}
        for i in range(len(iformat)):
            read_shape.append(zyx_shape[zyx_shape_map[iformat[i]]])

        # raw 1D array
        raw_float_data = np.fromfile(file_path, dtype='float32')
        # interpreted accoding to input shape
        raw_float_data = np.reshape(raw_float_data, tuple(read_shape))
        # reorder the axis to fit target (zyx) axis order
        self.data = np.transpose(raw_float_data, tuple(read_shape))

    # Reads the given ASCII file, which contents are interpreted as
    # according to the parameters into internal float ZYX numpy array.
    #
    # @icolumns is a tuple of 4 chars which define how to interpret each
    #   input file column. Possible char values: 'x', 'y', 'z', 'd', '-'
    #   (where xyz mean corresponding x value, and d means data value
    #   at given point, and - means ignore the column)
    def read_ascii(self, file_path, icolumns):
        # raw 2D array (y runs over file lines, x runs over columns)
        raw_data = np.loadtxt(file_path, dtype='float32')
        
        if (raw_data.shape[1] > len(icolumns)):
            raise ValueError("icolumns must describe all existing file columns")

        data_column_idx = 0
        for i in range(len(icolumns)):
            if (icolumns[i] == 'd'):
                data_column_idx = i
                break

        # TODO: one can add a continuity check for coordinate columns here

        self.data = raw_data(

        

        


# Does the convertion from input file to output file
# using given input and output formats
def convert(file_in, file_out, format_in, format_out):




parser = argparse.ArgumentParser(description='Converts from/to given ASCII format
        to/from raw float format to fit computational framework.')
parser.add_argument('--input', dest='input_file_path', type=str, nargs=1,
                    help='an input file to convert')
parser.add_argument('--output', dest='output_file_path', type=str, nargs=1,
                    help='an output file to write the result to')
parser.add_argument('--iformat', dest='input_format', type=str, nargs=1,
        help='an input file format; is understood as follows:\n'
             '* char prefix "b" for binary format\n'
             '  * order char triplet: "ccc" where every of "c" can be \n'
             '                  / |\\                                 \n'
             '                 C1 C2 C3                               \n'
             '    either "x", "y", "z", and the triplet order states  \n' 
             '    in which order the binary format is written:        \n' 
             '    first read/write swipe runs along the last (C3)     \n'
             '    axis, then C2 axis is incremented and runs          \n' 
             '    the second swipe along C3 axis, etc..               \n' 
             '* char prefix "a" for ascii format;                     \n'
             '  * order char triplet: "ccc" where every of "c" can be \n'
             '                  / |\\                                 \n'
             '                 C1 C2 C3                               \n'
             '    either "x", "y", "z". Defines the points order      \n' 
             '    while reading/writing the file, has the same        \n' 
             '    meaning as binary format order triplet.             \n' 
             '  * columns char triplet: "ccc" where every of "c" can  \n'
             '                  / |\\                                 \n'
             '                 C1 C2 C3                               \n'
             '    be either "x", "y", "z". It defines the ASCII format\n' 
             '    column order.                                       \n' 
             )

args = parser.parse_args()

print(args)



a =


a = array(floats,'float32')
output_file = open('file', 'wb')
a.tofile(output_file)
output_file.close()
