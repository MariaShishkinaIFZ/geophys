#!/bin/tcsh

#source coordinates
set xs = 1.0
set ys = 2.0
set zs = 0.2

trace_dw_only <<$END
receivers.xyz
0
dw1-1.tpk
0. -1.0
traveltimes.3d
$xs $ys $zs
modified_interface.2d
dw_raytraces.dat
$start_number
$END
