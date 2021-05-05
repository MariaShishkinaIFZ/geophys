#!/bin/tcsh

set xs = 1.0
set ys = 2.0
set zs = -0.2

trace_drw_only <<$END
receivers.xyz
0
drw1-1.tpk
0. -1.0
traveltimes.3d
$xs $ys $zs
modified_interface.2d
drw_raytraces.dat
drw_refraction_points.dat
modified_velmod.3d
modified_velmod.3d
0
$END
