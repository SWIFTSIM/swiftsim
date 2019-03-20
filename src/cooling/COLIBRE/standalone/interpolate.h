
int row_major_index_2d(const int x, const int y, const int Nx, const int Ny);
int row_major_index_3d(const int x, const int y, const int z, const int Nx, const int Ny, const int Nz);
int row_major_index_4d(const int x, const int y, const int z, const int w, 
                       const int Nx, const int Ny, const int Nz, const int Nw);
int row_major_index_5d(const int  x, const int  y, const int  z, const int  w, const int  v,
                       const int Nx, const int Ny, const int Nz, const int Nw, const int Nv);
void get_index_1d(const float *restrict table, const int size, const float x, int *i, float *restrict dx);

double interpolation4d(const float *table,
    const int xi, const int yi, const int zi, const int wi,
    const float dx, const float dy, const float dz, const float dw,
    const int Nx, const int Ny, const int Nz, const int Nw) ;

double interpolation4d_plus_summation(const float *table, const float *weights,
    const int istart, const int iend,
    const int xi, const int yi, const int zi, const int wi,
    const float dx, const float dy, const float dz, const float dw,
    const int Nx, const int Ny, const int Nz, const int Nw, const int Nv);

