#include <stdio.h>

typedef
  struct
    {
      int nx;
      int ny;
      int node_offset;
      double x_min;
      double x_max;
      double y_min;
      double y_max;
      double z_min;
      double z_max;
      double x_inc;
      double y_inc;
      double z_scale_factor;
      double z_add_offset;
      char x_units[80];
      char y_units[80];
      char z_units[80];
      char title[80];
      char command[320];
      char remark[160];
    } 
  GMT_native_header;


