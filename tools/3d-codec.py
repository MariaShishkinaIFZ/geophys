#!/usr/bin/env python3

import argparse
import traceback as tr
import os
import re
import importlib

current_source_dir = os.path.dirname(os.path.abspath(__file__))
components_dir_a = os.path.join(current_source_dir, "../components")

def load_geo_component(module_name):
    file_path = os.path.join(components_dir_a, "./" + module_name + ".py")
    try:
        spec = importlib.util.spec_from_file_location(module_name, file_path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
    except:
        import imp
        global vd
        vd = imp.load_source(module_name, file_path)

# local imports
load_geo_component("volumedata")

# Does the convertion from input file to output file
# using given input and output formats
# @config_in dictionary describing the input
#   * "file"
#     -> the file object or a path to the file to read
#   * "file_type"
#       -> "ascii" read the column-organized data from ASCII file
#       -> "binary" read the raw continuos fload data from binary file
#   * * binary parameters see @VolumeData.read_binary(...)
#   * * ascii parameters see @VolumeData.read_ascii(...)
# RETURNS:
#   True on success
#   str on failure (string error description)
def convert(config_in, config_out):
    vdata = None
    try:
        vdata = vd.VolumeData(config_in)
    except Exception as e:
        print("Failed to load data from file: %s, "
              "exception:\n"
              % (str(config_in["file"])))
        tr.print_exc()
        return "Input file read failure"

    try:
        vdata.write(config_out)
    except Exception as e:
        print("Failed to save data to file: %s, "
              "exception:\n"
              % (str(config_in["file"])))
        tr.print_exc()
        return "Output file write failure"

    return True


# Parses the arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Converts from/to given ASCII '
                                    'format to/from raw float format to fit '
                                    'computational framework.')
    parser.add_argument('--in', dest='input', type=str
                        , required=True, action='store'
                        , help='a path to an input file to convert')
    parser.add_argument('--out', dest='output', type=str
                        , required=True, action='store'
                        , help='a path to an output file to write the result to')
    parser.add_argument('--itype', dest='input_type', type=str
                        , choices=['binary', 'ascii'], required=True
                        , help='an input file type: can be one of "binary" or'
                               ' "ascii"')
    parser.add_argument('--otype', dest='output_type', type=str
                        , choices=['binary', 'ascii'], required=True
                        , help='an output file type: can be one of "binary" or'
                               ' "ascii"')
    parser.add_argument('--iformat', dest='input_format', type=str
                        , required=True
          , help='an input file format; is understood as follows:         \n'
                 '* When input file type (see "--itype") is binary:       \n'
                 '  * string of following format:                         \n'
                 '                   c(L)c(M)c(N)                         \n'
                 '    where "c" is one of "x","y","z" axis, and any axis, \n'
                 '    can be mentioned there only once. The meaning is    \n'
                 '    the coordinates change order in binary file: first  \n'
                 '    axis mentioned has lowest change rate, while the    \n'
                 '    last axis mentioned has the fastest change rate:    \n'
                 '    EXAMPLE: "x(L)z(M)y(N)" (ignore for now numbers     \n'
                 '      in braces) means that "y" index runs faster of    \n'
                 '      all when iterating over binary file, and the      \n'
                 '      "x" index runs slower of all: first point         \n'
                 '      in array has coords (ZYX order): (0; 0; 0)        \n'
                 '      second (0; 1; 0), third (0; 2; 0), ....           \n'
                 '      end of Y stride: (0; N-1; 0), next: (1; 0; 0),    \n'
                 '      ....                                              \n'
                 '      So jumping to next point in binary file changes   \n'
                 '      the last mentioned axis index, and then if it     \n'
                 '      reached its axis size limit, then intermediate    \n'
                 '      axis indes is changed, etc.                       \n'
                 '    Numbers L, M, N in braces define the corresponding  \n'
                 '    axis size.                                          \n'
                 '    EXAMPLE: "z(10)x(30)y(20)" says that binary file    \n'
                 '      Y axis iterates faster of all and has a size of   \n'
                 '      20 points, X axis is of intermediate change       \n'
                 '      rate and has a size of 30 points, and             \n'
                 '      Z axis iterates slower of all and has a size of   \n'
                 '      10 points.                                        \n'
                 '* When input file type (see "--itype") is ascii:        \n'
                 '  * column order string: "ccc...c" where every of "c"   \n'
                 '    can be either "x", "y", "z", "d" or "_".            \n'
                 '    this string defines how to interpret each input file\n'
                 '    column. "xyz" mean corresponding coordinate value,  \n'
                 '    "d" means data value, and "_" means "ignore this    \n'
                 '    column". --iformat string must describe **all**     \n'
                 '    ascii file columns and have only one occurence of   \n'
                 '    all chars except "_".                               \n'
                 '    EXAMPLE:                                            \n'
                 '          --iformat "_xy__zd"                           \n'
                 '                                                        \n'
                 '    NOTE: main assertion for ascii importing: the data  \n'
                 '      points grid is not rotated w.r.t. coordinate      \n'
                 '      frame (stride along any single index changes      \n'
                 '      coordinate only of a single axis).                \n'
                 '      Otherwise the automatic data grid parameters      \n'
                 '      computation will not work as it is.               \n'
                 )
    parser.add_argument('--oformat', dest='output_format', type=str
                        , required=True
          , help='an output file format; is understood as follows:        \n'
                 '* When output file type (see "--otype") is a binary:    \n'
                 '  * string of following format:                         \n'
                 '                       ccc                              \n'
                 '    where "c" is one of "x","y","z" axis, and any axis, \n'
                 '    can be mentioned there only once. The meaning is    \n'
                 '    the coordinates change order in binary file: first  \n'
                 '    axis mentioned has lowest change rate, while the    \n'
                 '    last axis mentioned has the fastest change rate:    \n'
                 '    EXAMPLE: "xzy" means that "y" index runs faster     \n'
                 '      of all when iterating over binary file, and the   \n'
                 '      "x" index runs slower of all: first point         \n'
                 '      in array has coords (ZYX order): (0; 0; 0)        \n'
                 '      second (0; 1; 0), third (0; 2; 0), ....           \n'
                 '      end of Y stride: (0; N-1; 0), next: (1; 0; 0),    \n'
                 '      ....                                              \n'
                 '      So jumping to next point in binary file changes   \n'
                 '      the last mentioned axis index, and then if it     \n'
                 '      reached its axis size limit, then intermediate    \n'
                 '      axis indes is changed, etc.                       \n'
                 '* When output file type (see "--otype") is ascii:       \n'
                 '  * column order string: "cccc" where every of "c"      \n'
                 '    can be either "x", "y", "z", "d".                   \n'
                 '    This string defines the ascii file column order.    \n'
                 '    "xyz" mean corresponding coordinate value,          \n'
                 '    "d" means data value. All 4 column types must be    \n'
                 '    mentioned exactly once.                             \n'
                 '    EXAMPLE:                                            \n'
                 '      --oformat "xyzd"                                  \n'
                 )

    return parser.parse_args()

# Arguments checks and conversions to fail early and clearly
# instead of long stack traced deep fails inside.
# @args the argparse arguments (parser.parse_args() result)
#
# RETURNS:
#   tuple (load configuration dict, save configuration dict)
#       on success
#   False on failure
#
# TODO: too big: split
def check_and_prepare_args(args):
    load_config = {}
    save_config = {}

    print("================= PARAMETERS ====================")
    print(args)
    print("=============== PARAMETERS END ==================")

    if not os.path.isfile(args.input):
        print("Input file should exist and be a regular file. "
              "Given input file: %s" % (str(args.input)))
        return False

    if not os.access(args.input, os.R_OK):
        print("Input file should be readable (now not the case) "
              "Current file: %s" % (str(args.input)))
        return False

    load_config["file"] = args.input
    save_config["file"] = args.output

    if (args.input_type == "binary"):
        load_config["file_type"] = "binary"

        # iformat: c(L)c(M)c(N)
        m = re.match(r"\A"
                     r"([xyz]{1})\((\d+)\)"
                     r"([xyz]{1})\((\d+)\)"
                     r"([xyz]{1})\((\d+)\)"
                     r"\Z"
                     , args.input_format)
        if not m:
            print("Input format (--iformat): %s doesn't match "
                  "expected pattern c(L)c(M)c(N). See help for "
                  "more information." % (str(args.input_format)))
            return False
        
        c1, c2, c3 = m.group(1), m.group(3), m.group(5)
        s1, s2, s3 = int(m.group(2)), int(m.group(4)), int(m.group(6))
        tr_map = {c1: s1 , c2: s2 , c3: s3}

        load_config["idx_change_order"] = c1 + c2 + c3
        load_config["volume_shape_zyx"] = (tr_map['z'], tr_map['y']
                                           , tr_map['x'])

    elif (args.input_type == "ascii"):
        load_config["file_type"] = "ascii"
        m = re.match(r"\A[xyzd_]{4,}\Z"
                     , args.input_format)
        if not m:
            print("Input format (--iformat): %s doesn't match "
                  "expected pattern for ascii input: only 'xyzd_' "
                  "symbols are allowed and all of 'xyzd' must be "
                  " mentioned. See help for more information."
                  % (str(args.input_format)))
            return False
        for c in "xyzd":
            m = re.match(r"\A[^%c]*%c{1}[^%c]*\Z" % (c,c,c)
                         , args.input_format)
            if not m:
                print("Input format (--iformat): %s doesn't match "
                      "expected pattern for ascii input: Only single "
                      "occurance of every char from 'xyzd' is allowed "
                      "and all must be mentioned. "
                      "See help for more information."
                      % (str(args.input_format)))
                return False
        load_config["input_columns"] = args.input_format.replace("_","-")
    else:
        print("Unknown input type format: %s" % (str(args.input_type)))
        return False

    if (args.output_type == "binary"):
        save_config["file_type"] = "binary"
        for c in "xyz":
            m = re.match(r"\A[^%c]*%c{1}[^%c]*\Z" % (c,c,c), args.output_format)
            if not m:
                print("Output format (--oformat): %s doesn't match "
                      "expected pattern for binary output: every column "
                      "('xyz') must be mentioned exactly once. See help "
                      "for more information." % (str(args.output_format)))
                return False
        save_config["idx_change_order"] = args.output_format
    elif (args.output_type == "ascii"):
        save_config["file_type"] = "ascii"

        for c in "xyzd":
            m = re.match(r"\A[^%c]*%c{1}[^%c]*\Z" % (c,c,c), args.output_format)
            if not m:
                print("Output format (--oformat): %s doesn't match "
                      "expected pattern for ascii input: every column ('xyzd') "
                      "must be mentioned exactly once. See help for more "
                      "information." % (str(args.output_format)))
                return False

        save_config["output_columns"] = args.output_format
    else:
        print("Unknown output type format: %s" % (str(args.output_type)))
        return False

    print("================= LOAD CONFIG ===================")
    print(load_config)
    print("=============== LOAD CONFIG END =================")
    print("================= SAVE CONFIG ===================")
    print(save_config)
    print("=============== SAVE CONFIG END =================")

    return (load_config, save_config)


configs = check_and_prepare_args(parse_args())

if not configs:
    print("Abort due to arguments processing failure.")
    exit(1)

result = convert(configs[0], configs[1])
if isinstance(result, str):
    print("Convertion failed.")
    exit(2)

exit(0)

