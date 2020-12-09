#!/usr/bin/env python3

import argparse
import unittest
import tempfile as tf
import random
import io
import os

current_source_dir = os.path.dirname(os.path.abspath(__file__))
test_data_dir_a = os.path.join(current_source_dir, "testdata")

try:
    import numpy as np
except Exception as e:
    raise Exception("Could not import NumPy, please install it:"
                    "* sudo apt-get install python3-numpy # on Linux"
                   ) from e
try:
    from scipy import stats
except Exception as e:
    raise Exception("Could not import SciPy, please install it:"
                    "* sudo apt-get install python3-scipy # on Linux"
                   ) from e

# Stores the spatial data in predefined format:
# ZYX, meaning that data[k][j][i] returns a point
# with coordinates (x == i, y == j, z == k). And
# the stride along X leads to reading out the
# neighboring points in memory.
class VolumeData:

    # the 3D numpy array containing the actual spatial data
    # NOTE: spatial coordinates of the data can be converted to the
    #   data indexes using the linear transform:
    #       r_index = r_spatial * grid2idx_dict[c][0] + grid2idx_dict[c][1]
    data = None
    # dict axis symbol -> tuple (a,b), where the a,b are the
    # linear coefficient regression from grid to index
    # NOTE: see also @data
    grid2idx_dict = {'x': None, 'y': None, 'z': None}

    # RETURNS:
    #   True if the object is file-like for us
    #   False else
    def __is_file_like(self, obj):
        return (isinstance(obj, io.IOBase)
                or (hasattr(obj, 'read') and hasattr(obj, 'write')))

    # initializes the VolumeData, see read(...)
    # description for load_config info
    def __init__(self, load_config = {}):
        self.read(load_config)

    # self.data contains the ZYX fload target data
    # @format_in the first char defines if th
    # @load_config the dictionary which contains the loading
    #   data parameters - values depend on the type of loaded data
    #   * "file"
    #     -> the file object or a path to the file to read
    #   * "file_type"
    #       -> "ascii" read the column-organized data from ASCII file
    #       -> "binary" read the raw continuos fload data from binary file
    #   * * binary parameters see @read_binary(...)
    #   * * ascii parameters see @read_ascii(...)
    #   OR the numpy array
    #   OR any entity which can be passed to constructor of numpy array
    def read(self, load_config = {}):
        if not isinstance(load_config, dict):
            self.data = np.array(load_config)
            print("Created empty volume data of shape: %s"
                  % (str(self.data.shape)))
            return

        if not "file_type" in load_config:
            self.data = np.array([[[]]])
            print("Created empty volume data.")
            return

        if load_config["file_type"] == "ascii":
            if self.__is_file_like(load_config["file"]):
                self.__read_ascii(load_config)
            else:
                with open(load_config["file"], mode="r") as f:
                    config = load_config.copy()
                    config["file"] = f
                    self.__read_ascii(config)
        elif load_config["file_type"] == "binary":
            if self.__is_file_like(load_config["file"]):
                self.__read_binary(load_config)
            else:
                with open(load_config["file"], mode="rb") as f:
                    config = load_config.copy()
                    config["file"] = f
                    self.__read_binary(config)

    # Saves the self.data in the appropriate format
    # defined by save_config
    # @save_config the dictionary which contains the saving
    #   data parameters - values depend on the type of loaded data
    #   * "file"
    #     -> the file object or a path to the file to write
    #   * "file_type"
    #       -> "ascii" read the column-organized data from ASCII file
    #       -> "binary" read the raw continuos fload data from binary file
    #   * * binary parameters see @write_binary(...)
    #   * * ascii parameters see @write_ascii(...)
    def write(self, save_config):
        if save_config["file_type"] == "ascii":
            if self.__is_file_like(save_config["file"]):
                self.__write_ascii(save_config)
            else:
                with open(save_config["file"], mode="w") as f:
                    config = save_config.copy()
                    config["file"] = f
                    self.__write_ascii(config)
        elif save_config["file_type"] == "binary":
            if self.__is_file_like(save_config["file"]):
                self.__write_binary(save_config)
            else:
                with open(save_config["file"], mode="wb") as f:
                    config = save_config.copy()
                    config["file"] = f
                    self.__write_binary(config)

    # Reads the given binary file, which contents are interpreted as
    # according to the parameters into internal float ZYX numpy array.
    #
    # @load_config the dictionary which contains the loading
    #   data parameters:
    #   * "file"
    #     -> the opened binary file object to read
    #   * "idx_change_order"
    #     -> tuple of incoming data axis order:
    #        the order of axis scan (last axis increased
    #        with every next point, second with every next
    #        stride, first with every next layer)
    #          Example: ('x', 'z', 'y') means that first point
    #            in array has coords (ZYX order): (0; 0; 0)
    #            second (0; 1; 0), third (0; 2; 0), ....
    #   * "volume_shape_zyx"
    #     -> a tuple of 3 sizes (for Z, Y, X axis correspondingly)
    #          NOTE: must fit to the actual binary data size
    #          Example: (30, 60, 120)
    #            means 30 z=const layers, 60 y=const strides,
    #            120 x=const points in every stride in incoming data
    def __read_binary(self, load_config):
        file_p = load_config["file"]
        iformat = load_config["idx_change_order"]
        zyx_shape_expected = load_config["volume_shape_zyx"]
        if len(iformat) != 3:
            raise ValueError("\"format\" must be 3 chars in size")
        if len(zyx_shape_expected) != 3:
            raise ValueError("\"volume_shape_zyx\" must be 3 items in size")

        print("Loading binary file: %s, format: %s, expected shape (ZYX): %s"
              % (str(file_p), str(iformat), str(zyx_shape_expected)))

        zyx_shape_map = {'z' : 0, 'y' : 1, 'x' : 2}
        read_shape = [zyx_shape_expected[zyx_shape_map[c]] for c in iformat]

        # raw 1D array
        raw_float_data = np.fromfile(file_p)

        input_shape_map = { iformat[0] : 0, iformat[1] : 1, iformat[2] : 2 }
        transpose_map = [input_shape_map[c] for c in "zyx"]

        # interpreted accoding to input shape
        raw_float_data = np.reshape(raw_float_data, tuple(read_shape))
        # reorder the axis to fit target (zyx) axis order
        self.data = np.transpose(raw_float_data, tuple(transpose_map))
        self.original_file_path = file_p

        shape = np.shape(self.data)
        print("Load done, points loaded: %d, layers: %d, layer size (YX): (%d, %d)"
              % (np.size(self.data), shape[0], shape[1], shape[2]))

    # Writes the binary file in our standard format:
    # float32, slowest axis Z, fastest axis X
    # @save_config the dictionary which contains the saving data parameters
    #   * "file"
    #     -> the opened binary file object to write
    #   * "idx_change_order"
    #     -> tuple of incoming data axis order:
    #        the order of axis scan (last axis increased
    #        with every next point, second with every next
    #        stride, first with every next layer)
    #          Example: ('x', 'z', 'y') means that first point
    #            in array has coords (ZYX order): (0; 0; 0)
    #            second (0; 1; 0), third (0; 2; 0), ....
    def __write_binary(self, config):
        file_p = config["file"]
        oformat = config["idx_change_order"]
        print("The volume data of shape (ZYX): %s, saving to the "
              "file: %s, in binary format (last index is fastest): %s"
              % (str(self.data.shape), str(file_p), str(oformat)))

        initial_map = { "z": 0, "y": 1, "x": 2 }
        transponent_map = tuple([initial_map[c] for c in oformat])
        out_data = np.transpose(self.data, transponent_map).ravel(order='C')

        out_data.tofile(file_p, sep="")
        print("The volume data of shape (ZYX): %s, was written to the "
              "file: %s, in binary format (last index is fastest): %s"
              % (str(self.data.shape), str(file_p), str(oformat)))

    # Computes the transformation from the given spatial uniform
    # grid to the indexes [0; len(grid_1d) - 1]:
    #           index = a * x + b
    # where x is the spatial coordinate on the grid, and also
    # checks if the grid is really uniform and nice.
    #
    # RETURNS:
    #   tuple(a,b) if the grid is uniform and OK
    #   raises ValueError if grid is bad (non-linear/non-uniform)
    def grid2idx_transform(self, grid_1d):
        if len(grid_1d.shape) != 1:
            raise ValueError("only 1d grid is supported for now, "
                             "provided shape: %s" % (grid_1d.shape))
        if grid_1d.shape[0] == 0:
            raise ValueError("no points in grid!")
        if grid_1d.shape[0] == 1:
            # linear regression complains if only one point
            # is given, so running it customly here
            return(1.0, float(grid_1d[0]))

        N = grid_1d.shape[0]
        res = stats.linregress(grid_1d, np.linspace(0, N - 1, N))
        a,b,r_val = res[0], res[1], res[2]

        EPSILON = 0.0001
        if abs(1.0 - r_val**2) > EPSILON:
            raise ValueError("Grid values are not uniform R^2 == %f, grid: "
                             "%s\n" % (r_val**2, str(grid_1d)))
        return (a,b)

    # @axis one of 'x','y','z' defines the axis to work with
    # RETURNS: the grid index for given spatial coordinate @spatial_coord, for
    #   given @axis
    def grid2idx(self, axis, spatial_coord):
        if self.grid2idx_dict[axis] is None:
            raise ValueError("the grid %c transform was not yet computed")
        a = self.grid2idx_dict[axis][0]
        b = self.grid2idx_dict[axis][1]
        return int(round(a * spatial_coord + b))

    # Reads the given ASCII file, which contents are interpreted as
    # according to the parameters into internal float ZYX numpy array.
    #
    # @load_config the dictionary which contains the loading
    #   data parameters:
    #   * "input_columns"
    #     -> is a tuple chars which define how to interpret each
    #        input file column. Possible char values: 'x', 'y', 'z', 'd', '-'
    #        (where xyz mean corresponding coordinate value, and d means data value
    #        at given point, and - means ignore the column)
    #        Example: ('x', 'z', '-', '-', 'd', 'y')
    #   * "file"
    #     -> the opened text file object to read
    #
    # Input data example:
    #       0.0     0.0     0.0     239.2
    #       1.0     0.0     0.0     124.4
    #       2.0     0.0     0.0     243.7
    #       0.0     1.0     0.0     384.1
    def __read_ascii(self, load_config):
        icolumns = load_config["input_columns"]
        file_p = load_config["file"]

        print("Loading data volume from ascii file: %s" % (file_p))
        print("Columns config: %s" % (str(icolumns)))

        # raw 2D array (y runs over file lines, x runs over columns)
        raw_data = np.loadtxt(file_p, dtype='float32')

        if (len(raw_data.shape) <= 0):
            raise ValueError("Nothing was loaded from file: %s "
                             "(resulting data shape: %s)"
                             % (str(file_p), str(raw_data.shape)))
        if (len(raw_data.shape) == 1):
            raw_data = np.reshape(raw_data, (1, raw_data.shape[0]))
        print("Loaded raw ascii data shape: %s" % (str(raw_data.shape)))

        if (raw_data.shape[1] > len(icolumns)):
            raise ValueError("\"input_columns\" must describe "
                             "all existing file columns")
        points_count = raw_data.shape[0]
        columns_count = raw_data.shape[1]
        if (points_count <= 0 or columns_count <= 0):
            raise ValueError("empty input data comes from file: %s, "
                             "incoming data shape: %s"
                             % (str(file_p), str(raw_data.shape)))

        # the indexes of target columns
        zyxd_column_indexes = []
        for c in 'zyxd':
            idx = [index for index,val in enumerate(icolumns) if val == c]
            if (len(idx) != 1):
                raise ValueError("\"input_columns\" must contain "
                                 "exactly one flag '%c' now has: %d"
                                 % (c, len(data_column_idx)))
            zyxd_column_indexes.append(idx[0])

        # only relevant columns
        needed_data = raw_data[:, tuple(zyxd_column_indexes)]

        # TODO: a continuity check of the given data
        print("WARNING: continuity and uniformity checks are NOT performed "
              "on the input data!")

        # NOTE: main assertion for now: the data points grid is
        #   not rotated w.r.t. coordinate frame (stride along
        #   any single index changes spatial coordinate only of a
        #   single axis).

        # Example of unique grid vals:  [-1.3, 0.7, 2.7, 4.7]
        unique = {'z': np.unique(needed_data[:, 0])
                  , 'y': np.unique(needed_data[:, 1])
                  , 'x': np.unique(needed_data[:, 2])}

        print("Extracted data coordinate planes count (Z;Y;X): (%d;%d;%d)"
              % tuple([unique[c].shape[0] for c in 'zyx']))

        # now we need to compute the appropriate tranformation
        # from the spatial grid coordinates to array indexes
        for c in 'xyz':
            try:
                self.grid2idx_dict[c] = self.grid2idx_transform(unique[c])
                print("Grid %c space to idx: %c_idx = %f * %c + %f"
                      % (c, c, self.grid2idx_dict[c][0], c
                         , self.grid2idx_dict[c][1]))
            except Exception as e:
                raise ValueError("(the grid: %s); please fix it before "
                                 "continuing. Abort." % (c)) from e

        data_shape = np.array([ unique[c].shape[0] for c in 'zyx' ], dtype=int)

        print("Data volume index shape (ZYX): %s" % (str(data_shape)))
        print("Data points total count: %d" % (points_count))

        self.data = np.zeros(data_shape)

        for i in range(points_count):
            rec = needed_data[i]
            z = rec[0]
            y = rec[1]
            x = rec[2]
            val = rec[3]
            z_idx = self.grid2idx('z', z)
            y_idx = self.grid2idx('y', y)
            x_idx = self.grid2idx('x', x)
            self.data[z_idx][y_idx][x_idx] = val
        self.original_file_path = file_p

        shape = np.shape(self.data)
        print("Load done, points loaded: %d, layers: %d, layer size (XY): (%d, %d)"
              % (np.size(self.data), shape[0], shape[1], shape[2]))

    # Writes the volume data according to the given config.
    # @config the dictionary which contains the saving configuration
    #   parameters - values depend on the type of loaded
    #   * "output_columns"
    #     -> is a tuple chars which define how to interpret each
    #        input file column. Possible char values: 'x', 'y', 'z', 'd'
    #        (where xyz mean corresponding coordinate value, and d means data value
    #        at given point)
    #        Example: ('x', 'z', 'd', 'y')
    #   * "file"
    #     -> the opened text file object to write
    def __write_ascii(self, config):
        col_order = config["output_columns"]
        file_p = config["file"]
        print("The volume data of shape (ZYX): %s, saving to the "
              "file: %s, in ascii format (columns order): %s"
              % (str(self.data.shape), str(file_p), str(col_order)))

        NX = self.data.shape[0]
        NY = self.data.shape[1]
        NZ = self.data.shape[2]
        points_count = NX * NY * NZ
        col_count = len(col_order)

        points = np.empty((points_count, col_count), dtype=float)

        zz, yy, xx = np.meshgrid(np.linspace(0, NZ - 1, NZ)
                                 , np.linspace(0, NY - 1, NY)
                                 , np.linspace(0, NX - 1, NX)
                                 , indexing='ij')
        coords = np.empty((points_count, 3))
        coords[:,0] = np.reshape(xx, (points_count, ))
        coords[:,1] = np.reshape(yy, (points_count, ))
        coords[:,2] = np.reshape(zz, (points_count, ))

        for raw_idx in range(coords.shape[0]):
            k = int(coords[raw_idx][0])
            j = int(coords[raw_idx][1])
            i = int(coords[raw_idx][2])
            raw = np.empty((1, col_count), dtype=float)
            for col_idx in range(col_count):
                col_type = col_order[col_idx]
                if col_type == 'x': raw[0][col_idx] = i
                elif col_type == 'y': raw[0][col_idx] = j
                elif col_type == 'z': raw[0][col_idx] = k
                elif col_type == 'd':
                    raw[0][col_idx] = self.data[k][j][i]
                else:
                    raw[0][col_idx] = rn.randint(MIN_MOCK_DATA
                                                 , MAX_MOCK_DATA)
            points[raw_idx, :] = raw

        for idx in range(points_count):
            points[idx].tofile(file_p, sep="\t")
            file_p.write("\n")

        print("The volume data of shape (ZYX): %s, saved to the "
              "file: %s, in ascii format (columns order): %s"
              % (str(self.data.shape), str(file_p), str(col_order)))

#########################
# TESTS
#########################

# Tests the functionality of the VolumeData class IF
class VolumeDataTester(unittest.TestCase):

    def test_construction_ascii_simple(self):
        with tf.SpooledTemporaryFile(max_size=1000, mode='r+') as in_file:
            in_file.write("0 0 0 100\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": "xyzd" })
            # data geometry: ZYX (1 1 1)
            self.assertEqual(len(np.shape(vd.data)), 3)
            self.assertEqual(np.shape(vd.data)[0], 1)
            self.assertEqual(np.shape(vd.data)[1], 1)
            self.assertEqual(np.shape(vd.data)[2], 1)

            # data value
            self.assertEqual(vd.data[0][0][0], 100.)

    def test_construction_ascii_3points(self):
        with tf.SpooledTemporaryFile(max_size=1000, mode='r+') as in_file:
            in_file.write("0 0 0 100\n")
            in_file.write("0 0 1 200\n")
            in_file.write("0 0 2 300\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             ,"input_columns": "xyzd" })
            # data geometry: ZYX (3 1 1)
            self.assertEqual(len(np.shape(vd.data)), 3)
            self.assertEqual(np.shape(vd.data)[0], 3)
            self.assertEqual(np.shape(vd.data)[1], 1)
            self.assertEqual(np.shape(vd.data)[2], 1)

            # data value
            self.assertEqual(vd.data[0][0][0], 100.)
            self.assertEqual(vd.data[1][0][0], 200.)
            self.assertEqual(vd.data[2][0][0], 300.)

    def test_construction_ascii_3points_zyxd(self):
        with tf.SpooledTemporaryFile(max_size=1000, mode='r+') as in_file:
            in_file.write("0 0 0 100\n")
            in_file.write("0 0 1 200\n")
            in_file.write("0 0 2 300\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             ,"input_columns": "zyxd" })
            # data geometry: ZYX (3 1 1)
            self.assertEqual(len(np.shape(vd.data)), 3)
            self.assertEqual(np.shape(vd.data)[0], 1)
            self.assertEqual(np.shape(vd.data)[1], 1)
            self.assertEqual(np.shape(vd.data)[2], 3)

            # data value
            self.assertEqual(vd.data[0][0][0], 100.)
            self.assertEqual(vd.data[0][0][1], 200.)
            self.assertEqual(vd.data[0][0][2], 300.)

    def test_construction_ascii_12points__yd_zx_(self):
        with tf.SpooledTemporaryFile(max_size=1000, mode='r+') as in_file:
            # ZYX (3 2 2)
            #              -  y   d   -  z x  -
            in_file.write("80 0 1000 321 0 0 101\n")
            in_file.write("81 1 1010 322 0 0 102\n")
            in_file.write("82 1 1011 324 0 1 104\n")
            in_file.write("83 0 1100 325 1 0 105\n")
            in_file.write("84 0 1101 327 1 1 107\n")
            in_file.write("85 0 1001 323 0 1 103\n")
            in_file.write("86 1 1111 328 1 1 108\n")
            in_file.write("87 0 1200 329 2 0 109\n")
            in_file.write("88 1 1210 3210 2 0 1010\n")
            in_file.write("89 0 1201 3211 2 1 1011\n")
            in_file.write("90 1 1110 326 1 0 106\n")
            in_file.write("91 1 1211 3212 2 1 1012\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             ,"input_columns": "-yd-zx-" })
            # data geometry: ZYX (3 2 2)
            self.assertEqual(len(np.shape(vd.data)), 3)
            self.assertEqual(np.shape(vd.data)[0], 3)
            self.assertEqual(np.shape(vd.data)[1], 2)
            self.assertEqual(np.shape(vd.data)[2], 2)

            # data value
            self.assertEqual(vd.data[0][0][0], 1000.)
            self.assertEqual(vd.data[0][0][1], 1001.)
            self.assertEqual(vd.data[0][1][0], 1010.)
            self.assertEqual(vd.data[0][1][1], 1011.)
            self.assertEqual(vd.data[1][0][0], 1100.)
            self.assertEqual(vd.data[1][0][1], 1101.)
            self.assertEqual(vd.data[1][1][0], 1110.)
            self.assertEqual(vd.data[1][1][1], 1111.)
            self.assertEqual(vd.data[2][0][0], 1200.)
            self.assertEqual(vd.data[2][0][1], 1201.)
            self.assertEqual(vd.data[2][1][0], 1210.)
            self.assertEqual(vd.data[2][1][1], 1211.)

    def test_construction_ascii_simple_to_binary(self):
        with tf.SpooledTemporaryFile(max_size=1000
                                     , mode='r+') as in_file:
            in_file.write("0 0 0 42\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": "xyzd" })
            with tf.SpooledTemporaryFile(max_size=1000
                                         , mode='r+b') as out_file:
                vd.write({"file": out_file
                          , "file_type": 'binary'
                          , "idx_change_order": "zyx"})
                out_file.seek(0)
                vd2 = VolumeData({"file": out_file
                                  , "file_type": "binary"
                                  , "idx_change_order": "zyx"
                                  , "volume_shape_zyx": (1, 1, 1)})

                # data geometry: ZYX (1 1 1)
                self.assertEqual(len(np.shape(vd2.data)), 3)
                self.assertEqual(np.shape(vd2.data)[0], 1)
                self.assertEqual(np.shape(vd2.data)[1], 1)
                self.assertEqual(np.shape(vd2.data)[2], 1)

                # data value
                self.assertEqual(vd2.data[0][0][0], 42.)

    def test_construction_ascii_3linepoints_to_binary(self):
        with tf.SpooledTemporaryFile(max_size=1000
                                     , mode='r+') as in_file:
            #              Y Z X
            in_file.write("0 0 0 15\n")
            in_file.write("0 0 1 25\n")
            in_file.write("0 0 2 35\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": "yzxd" })
            with tf.SpooledTemporaryFile(max_size=1000
                                         , mode='r+b') as out_file:
                vd.write({"file": out_file
                          , "file_type": 'binary'
                          , "idx_change_order": "xzy"})
                out_file.seek(0)
                vd2 = VolumeData({"file": out_file
                                  , "file_type": "binary"
                                  , "idx_change_order": "xzy"
                                  , "volume_shape_zyx": (1, 1, 3)})

                # data geometry: ZYX (1 1 3)
                self.assertEqual(len(np.shape(vd2.data)), 3)
                self.assertEqual(np.shape(vd2.data)[0], 1)
                self.assertEqual(np.shape(vd2.data)[1], 1)
                self.assertEqual(np.shape(vd2.data)[2], 3)

                # data value
                self.assertEqual(vd2.data[0][0][0], 15.)
                self.assertEqual(vd2.data[0][0][1], 25.)
                self.assertEqual(vd2.data[0][0][2], 35.)

    def test_construction_ascii_1x2x3_to_binary(self):
        with tf.SpooledTemporaryFile(max_size=1000
                                     , mode='r+') as in_file:
            #              Y Z X
            in_file.write("0 0 0 1\n")
            in_file.write("0 0 1 2\n")
            in_file.write("0 0 2 3\n")
            in_file.write("1 0 0 4\n")
            in_file.write("1 0 1 5\n")
            in_file.write("1 0 2 6\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": "yzxd" })
            with tf.SpooledTemporaryFile(max_size=1000
                                         , mode='r+b') as out_file:
                vd.write({"file": out_file
                          , "file_type": 'binary'
                          , "idx_change_order": "zxy"})
                out_file.seek(0)
                expected_binary_file_data = (b'\x00\x00\x00\x00\x00\x00\xf0?'
                                             b'\x00\x00\x00\x00\x00\x00\x10@'
                                             b'\x00\x00\x00\x00\x00\x00\x00@'
                                             b'\x00\x00\x00\x00\x00\x00\x14@'
                                             b'\x00\x00\x00\x00\x00\x00\x08@'
                                             b'\x00\x00\x00\x00\x00\x00\x18@')
                binary_data = out_file.read()
                self.assertEqual(len(binary_data), len(expected_binary_file_data))
                self.assertEqual(binary_data, expected_binary_file_data)

    def test_construction_1x2x3_from_zxy_binary(self):
        with tf.SpooledTemporaryFile(max_size=1000
                                     , mode='r+b') as in_file:
            # effectively (in format "zxy"):
            # Z Y X D
            # 0 0 0 1.
            # 0 1 0 4.
            # 0 0 1 2.
            # 0 1 1 5.
            # 0 0 2 3.
            # 0 1 2 6.
            binary_file_data = (b'\x00\x00\x00\x00\x00\x00\xf0?'
                                b'\x00\x00\x00\x00\x00\x00\x10@'
                                b'\x00\x00\x00\x00\x00\x00\x00@'
                                b'\x00\x00\x00\x00\x00\x00\x14@'
                                b'\x00\x00\x00\x00\x00\x00\x08@'
                                b'\x00\x00\x00\x00\x00\x00\x18@')
            in_file.write(binary_file_data)
            in_file.seek(0)

            vd = VolumeData({"file": in_file
                             , "file_type": "binary"
                             , "idx_change_order": "zxy"
                             , "volume_shape_zyx": (1, 2, 3)})

            # data geometry: ZYX (1 2 3)
            self.assertEqual(len(np.shape(vd.data)), 3)
            self.assertEqual(np.shape(vd.data)[0], 1)
            self.assertEqual(np.shape(vd.data)[1], 2)
            self.assertEqual(np.shape(vd.data)[2], 3)

            # data value
            self.assertEqual(vd.data[0][0][0], 1.)
            self.assertEqual(vd.data[0][0][1], 2.)
            self.assertEqual(vd.data[0][0][2], 3.)
            self.assertEqual(vd.data[0][1][0], 4.)
            self.assertEqual(vd.data[0][1][1], 5.)
            self.assertEqual(vd.data[0][1][2], 6.)

    def test_construction_ascii_4points_to_from_binary(self):
        with tf.SpooledTemporaryFile(max_size=1000
                                     , mode='r+') as in_file:
            #              Y Z X
            in_file.write("0 0 0 42\n")
            in_file.write("0 0 1 52\n")
            in_file.write("1 0 0 62\n")
            in_file.write("1 0 1 72\n")
            in_file.seek(0)
            print(in_file.read())
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": "yzxd" })
            with tf.SpooledTemporaryFile(max_size=1000
                                         , mode='r+b') as out_file:
                vd.write({"file": out_file
                          , "file_type": 'binary'
                          , "idx_change_order": "zxy"})
                out_file.seek(0)
                vd2 = VolumeData({"file": out_file
                                  , "file_type": "binary"
                                  , "idx_change_order": "zxy"
                                  , "volume_shape_zyx": (1, 2, 2)})

                # data geometry: ZYX (1 2 2)
                self.assertEqual(len(np.shape(vd2.data)), 3)
                self.assertEqual(np.shape(vd2.data)[0], 1)
                self.assertEqual(np.shape(vd2.data)[1], 2)
                self.assertEqual(np.shape(vd2.data)[2], 2)

                # data value
                self.assertEqual(vd2.data[0][0][0], 42.)
                self.assertEqual(vd2.data[0][0][1], 52.)
                self.assertEqual(vd2.data[0][1][0], 62.)
                self.assertEqual(vd2.data[0][1][1], 72.)


    # generates the random shaped data, and random ascii
    # order and checks if loaded data matched
    #
    # RETURNS:
    #   (written volume data in ZYX format, columns order tuple)
    #
    # TODO: looks relatively ugly - refactor
    def helper_gen_data_ascii(self, max_shape_zyx, min_shape_zyx
                              , min_extra_columns, max_extra_columns
                              , file_out):
        MIN_MOCK_DATA = -1000
        MAX_MOCK_DATA = 1000

        import random as rn
        data_zyx_shape = (rn.randint(min_shape_zyx[0], max_shape_zyx[0])
                          ,rn.randint(min_shape_zyx[1],max_shape_zyx[1])
                          ,rn.randint(min_shape_zyx[2],max_shape_zyx[2]))
        NX = data_zyx_shape[0]
        NY = data_zyx_shape[1]
        NZ = data_zyx_shape[2]
        points_count = NX * NY * NZ

        print("Test data shape (ZYX): %s" % (str(data_zyx_shape)))
        test_data_zyx = np.random.random((data_zyx_shape)) * 100
        print("Test data: %s" % (str(test_data_zyx)))

        extra_columns_count = rn.randint(min_extra_columns
                                         , max_extra_columns)

        orig_col_order = "xyzd" + extra_columns_count * "-"
        col_order = ''.join(random.sample(orig_col_order, len(orig_col_order)))
        print("ASCII extra colums count: %d" % (extra_columns_count))
        print("ASCII columns: %s" % (str(col_order)))

        # major data to write with 4 columns (3 coords + value
        # + extra_ignore_columns)
        col_count = len(col_order)
        points = np.empty((points_count, col_count), dtype=float)

        zz, yy, xx = np.meshgrid(np.linspace(0, NZ - 1, NZ)
                                 , np.linspace(0, NY - 1, NY)
                                 , np.linspace(0, NX - 1, NX)
                                 , indexing='ij')
        coords = np.empty((points_count, 3))
        coords[:,0] = np.reshape(xx, (points_count, ))
        coords[:,1] = np.reshape(yy, (points_count, ))
        coords[:,2] = np.reshape(zz, (points_count, ))

        np.random.shuffle(coords)

        for raw_idx in range(coords.shape[0]):
            k = int(coords[raw_idx][0])
            j = int(coords[raw_idx][1])
            i = int(coords[raw_idx][2])
            raw = np.empty((1, col_count), dtype=float)
            for col_idx in range(col_count):
                col_type = col_order[col_idx]
                if col_type == 'x': raw[0][col_idx] = i
                elif col_type == 'y': raw[0][col_idx] = j
                elif col_type == 'z': raw[0][col_idx] = k
                elif col_type == 'd':
                    raw[0][col_idx] = test_data_zyx[k][j][i]
                else:
                    raw[0][col_idx] = rn.randint(MIN_MOCK_DATA
                                                 , MAX_MOCK_DATA)
            points[raw_idx, :] = raw
        print("ASCII data: shape %s" % (str(points.shape)))
        for idx in range(points_count):
            points[idx].tofile(file_out, sep="\t")
            file_out.write("\n")

        return (test_data_zyx, col_order)

    # helper to test if the intentionally generated ASCII file
    # is loaded properly into the parsing class
    def helper_check_ascii_read_correct(
            self, max_shape_zyx, min_shape_zyx, min_extra_columns
            , max_extra_columns):
        with tf.SpooledTemporaryFile(max_size=10000
                                     , mode='r+') as in_file:

            expected_data, col_oder = self.helper_gen_data_ascii(
                                            max_shape_zyx
                                            , min_shape_zyx
                                            , min_extra_columns
                                            , max_extra_columns
                                            , in_file)
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": col_oder })

            self.assertTrue(np.allclose(vd.data, expected_data)
                            , msg=("\nLoaded data:\n%s\nExpected data:\n%s\n"
                                   %(vd.data, expected_data)))

    # helper to test if the intentionally generated ASCII file
    # is loaded properly into the parsing class, then saved into
    # a binary, and then loaded back from the binary correctly
    def helper_check_binary_write_read_correct(
            self, max_shape_zyx, min_shape_zyx):
        with tf.SpooledTemporaryFile(max_size=10000
                                     , mode='r+') as in_file:

            expected_data, col_oder = self.helper_gen_data_ascii(
                                            max_shape_zyx
                                            , min_shape_zyx
                                            , 0, 0, in_file)
            in_file.seek(0)
            vd = VolumeData({"file": in_file
                             , "file_type": "ascii"
                             , "input_columns": col_oder })

            # ASCII was correctly loaded
            self.assertTrue(np.allclose(vd.data, expected_data)
                            , msg=("\nLoaded data:\n%s\nExpected data:\n%s\n"
                                   %(vd.data, expected_data)))

            with tf.SpooledTemporaryFile(max_size=1000
                                         , mode='r+b') as out_file:
                change_order = "xyz"
                change_order = ''.join(random.sample(change_order
                                                     , len(change_order)))
                vd.write({"file": out_file
                          , "file_type": 'binary'
                          , "idx_change_order": change_order})

                out_file.seek(0)
                vd2 = VolumeData({"file": out_file
                                  , "file_type": "binary"
                                  , "idx_change_order": change_order
                                  , "volume_shape_zyx": expected_data.shape})

                self.assertTrue(np.allclose(vd2.data, expected_data)
                                , msg=("\nLoaded from binary data:\n"
                                       "%s\nExpected data:\n%s\n"
                                       %(vd2.data, expected_data)))


    # Quite exhaustive test for ascii parsing
    def test_ascii_read_random_data_100iterations(self):
        for i in range(100):
            self.helper_check_ascii_read_correct(
                       (10,10,10), (1,1,1), 0, 3)

    # Quite exhaustive test for ascii parsing
    def test_ascii_binary_rw_random_data_100iterations(self):
        for i in range(100):
            self.helper_check_binary_write_read_correct(
                                    (10,10,10), (1,1,1))

    # tests direct data creation
    def test_direct_creation(self):
        expected_data = np.array([[[42.3, 52.4], [62.3, 72.8]]])

        vd = VolumeData()
        vd.data = expected_data
        self.assertTrue(np.allclose(vd.data, expected_data)
                        , msg=("\nVolumeData assigned:\n"
                               "%s\ndiffers from expected data:\n%s\n"
                               %(vd.data, expected_data)))

        vd = VolumeData(np.array([[[42.3, 52.4], [62.3, 72.8]]]))
        self.assertTrue(np.allclose(vd.data, expected_data)
                        , msg=("\nDirect constructor data (from ndarray):\n"
                               "%s\ndiffers from expected data:\n%s\n"
                               %(vd.data, expected_data)))

        vd = VolumeData([[[42.3, 52.4], [62.3, 72.8]]])
        self.assertTrue(np.allclose(vd.data, expected_data)
                        , msg=("\nDirect constructor data (from list):\n"
                               "%s\ndiffers from expected data:\n%s\n"
                               %(vd.data, expected_data)))

    # tests how we write ascii file:
    #     VolumeData -> automatic ascii ->  VolumeData
    def test_ascii_write(self):
        with tf.SpooledTemporaryFile(max_size=1000
                                     , mode='r+') as in_file:
            expected_data = np.array([[[42.3, 52.4], [62.3, 72.8]]])

            vd = VolumeData()
            vd.data = expected_data

            with tf.SpooledTemporaryFile(max_size=1000
                                         , mode='r+') as saved_file:
                col_order = "ydzx"
                vd.write({"file": saved_file
                          , "file_type": "ascii"
                          , "output_columns": col_order })
                saved_file.seek(0)

                vd2 = VolumeData({"file": saved_file
                                  , "file_type": "ascii"
                                  , "input_columns": col_order })

                self.assertTrue(np.allclose(vd2.data, expected_data)
                                    , msg=("\nLoaded from ascii data:\n"
                                           "%s\nExpected data:\n%s\n"
                                           %(vd.data, expected_data)))

    def test_construction_sparse_ascii_to_from_binary(self):
        in_ascii_file_path = os.path.join(test_data_dir_a
                                          , "volumedata_test_ascii_1.txt")
        vd = VolumeData({"file": in_ascii_file_path
                         , "file_type": "ascii"
                         , "input_columns": "xyzd"})
        self.assertTrue(vd is not None)

        # see test data file (this is part of the data to compare)
        expected_data = np.array([2551.81, 2552.82, 2553.83, 2554.84
                                  , 2555.85 , 2556.86 , 2557.87 , 2558.88
                                  , 2559.89 , 2550.90 , 2551.91])

        real_data = vd.data[0,0:11,0].flatten()
        self.assertTrue(np.allclose(real_data, expected_data)
                        , msg=("\nLoaded from ascii data (ZYX) [0,0:11,0]:\n"
                               "%s\nExpected data:\n%s\n"
                               %(real_data, expected_data)))


if __name__ == "__main__":
    unittest.main()

