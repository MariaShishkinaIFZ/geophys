#!/bin/tcsh

findheadw << END
rec.xyz
0
hw_s1.tp
0. 0.01
th_s1.3d
14.5 5.0 0.0
plane15.2d
hw-s1-sel.tp
dw-s1-sel.tp
vel_upl.3d
velcr.3d
0
END
