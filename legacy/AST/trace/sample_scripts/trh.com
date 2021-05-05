#!/bin/tcsh

#source coordinates
set xs = 1.0
set ys = 2.0
set zs = 0.2

trace_hw_only_1 <<$END
receivers..xyz
0
hw1-1.tpk
0. -1.0
traveltimes.3d
$xs $ys $zs
modified_interface.2d
hw_raytraces.dat
hw_refraction_points.dat
velmod.3d
velmod.3d
0
$END
