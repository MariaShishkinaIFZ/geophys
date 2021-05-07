#ifndef H_ISDGIO
#define H_ISDGIO

void read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *filename);

void write_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *dir, const char *filename);

void read_2D_surface(float **Buf, long Ny, long Nx, const char *filename);

void write_2D_surface(float **Buf, long Ny, long Nx, const char *dir, const char *filename);

#endif