#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/

// See: FILE_COLUMNS_MAP parameter
//static const char DEFAULT_FILE_COLUMNS_MAP[] = "xyzd";

void help()
{
    printf("USAGE:\n");
    printf("    text23d <velocyty_model.txt> <velocyty_model.3d> [FILE_COLUMNS_MAP]\n");
    printf("\n");
    printf("PARAMETERS:\n");
    printf("    FILE_COLUMNS_MAP:  the string consisting of combination\n");
    printf("                       of 4 characters:  x y z d -         \n");
    printf("                       It defines the semantic of column   \n");
    printf("                       values of the file.                 \n");
    printf("                       VALUES:                             \n");
    printf("                           x - the X coordinate.           \n");
    printf("                           y - the Y coordinate.           \n");
    printf("                           z - the Z coordinate.           \n");
    printf("                           d - the D coordinate.           \n");
    printf("                           - - the ignore column.          \n");
    printf("                       DEFAULT:                            \n");
    printf("                                'xyzd'                     \n");
    printf("\n");
    printf("RETURNS:\n");
    printf("    0 : on success\n");
    printf("    1 : if arguments are wrong\n");
    printf("    2 : other error\n");
}

int main(int narg, char **argv)
{
  long int Nx, Ny, Nz;
  int i, j, k;
  float x0, y0, z0, h;
  float ***Vu;
  FILE *fp;

  // Maps the file column number (starting from 0) to the data (x, y, z, DATA):
  //    column_index -> data
  //
  // data: {'x', 'y', 'z', 'd', any_other}
  //
  // NOTE: any_other character means "ignore column"
  // NOTE: is set to "xyzd" by default.
//  char *file_columns_map;

  if (narg < 2 || narg > 3) {
    printf("Weired number of arguments. See usage!");
    help();
    return 1;
  }
  char *output_file_name = argv[2];

//  if (narg == 3) {
//    file_columns_map = argv[3];
//  } else {
//    file_columns_map = DEFAULT_FILE_COLUMNS_MAP;
//  }

  //char format_str[] = "";
  
  fp = file_open(GEOMFILE, "rt");
  fscanf(fp, "%ld %ld %ld", &Nx, &Ny, &Nz);
  fscanf(fp, "%g %g %g", &x0, &y0, &z0);
  fscanf(fp, "%g", &h);
  fclose(fp);
 
  // According to general agreements:
  //    * Z axes has 0 on the ground surface, and direction toward Earth center
  alloc_float_3d(&Vu, Nz, Ny, Nx);
 
  fp = file_open(argv[1], "rt");
  
  for (k = 0; k < Nz; k++) {
    for (j = 0; j < Nx; j++)
      for (i = 0; i < Ny; i++)
        fscanf(fp, "%g %g %g %g \n", &x0, &y0, &z0, &Vu[k][i][j]);
      }

  write_3D_volume(Vu, Nz, Ny, Nx, output_file_name);
  
  free_float_3d(Vu, Nz, Ny, Nx);
  exit(1);
}
