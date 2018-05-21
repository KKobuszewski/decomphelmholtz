#ifndef __HELMHOLZ_DECOMP__
#define __HELMHOLZ_DECOMP__

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <time.h> /* for ctime() */
#include <unistd.h> /* for sleep() */
#include <fftw3.h> /* for FFTW */


// error codes
#define HELMHOLZ_SUCCESS 0


// global variables
fftw_complex *zin, *zout;
fftw_plan plan_f, plan_b;
int *ixyz2ix, *ixyz2iy, *ixyz2iz;
double *kx, *ky, *kz;
double dkx, dky, dkz;
double *_Acr_x; double *_Acr_y; double *_Acr_z;
double *_Air_x; double *_Air_y; double *_Air_z;
double complex *_zAk_x, *_zAck_x, *_zAik_x;
double complex *_zAk_y, *_zAck_y, *_zAik_y;
double complex *_zAk_z, *_zAck_z, *_zAik_z;

int flag_initialized;

int init_helmholz_decomp(const int nx, const int ny, const int nz)
{
    const int nxyz = nx*ny*nz;
    
    // Prepare for FFTW
    zin  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nxyz);
    zout = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nxyz);
    plan_f = fftw_plan_dft_3d(nx,ny,nz,zin,zout,FFTW_FORWARD ,FFTW_ESTIMATE);
    plan_b = fftw_plan_dft_3d(nx,ny,nz,zin,zout,FFTW_BACKWARD,FFTW_ESTIMATE);
    
    _Acr_x = (double*)         malloc( nxyz * sizeof(double) );
    _Acr_y = (double*)         malloc( nxyz * sizeof(double) );
    _Acr_z = (double*)         malloc( nxyz * sizeof(double) );
    _Air_x = (double*)         malloc( nxyz * sizeof(double) );
    _Air_y = (double*)         malloc( nxyz * sizeof(double) );
    _Air_z = (double*)         malloc( nxyz * sizeof(double) );
    _zAk_x  = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAk_y  = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAk_z  = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAck_x = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAck_y = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAck_z = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAik_x = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAik_y = (double complex*) malloc( nxyz * sizeof(double complex) );
    _zAik_z = (double complex*) malloc( nxyz * sizeof(double complex) );
    
    ixyz2ix = (int*) malloc( nxyz * sizeof(int) );
    ixyz2iy = (int*) malloc( nxyz * sizeof(int) );
    ixyz2iz = (int*) malloc( nxyz * sizeof(int) );
    int ixyz=0;
    for (int ix=0;ix<nx;ix++) for (int iy=0;iy<ny;iy++) for (int iz=0;iz<nz;iz++)
    {
      ixyz2ix[ixyz]=ix; ixyz2iy[ixyz]=iy; ixyz2iz[ixyz]=iz;
      ixyz++;
    }
    
    dkx = 2.0*M_PI/( double )nx;
    dky = 2.0*M_PI/( double )ny;
    dkz = 2.0*M_PI/( double )nz;
    kx = (double*) malloc( nx * sizeof(double) );
    ky = (double*) malloc( ny * sizeof(double) );
    kz = (double*) malloc( nz * sizeof(double) );
    
    for(int ix=0;ix<nx;ix++)
    {
      if(ix<nx/2) kx[ix] = dkx*( double )(ix   );
      else        kx[ix] = dkx*( double )(ix-nx);
    }
    for(int iy=0;iy<ny;iy++)
    {
      if(iy<ny/2) ky[iy] = dky*( double )(iy   );
      else        ky[iy] = dky*( double )(iy-ny);
    }
    for(int iz=0;iz<nz;iz++)
    {
      if(iz<nz/2) kz[iz] = dkz*( double )(iz   );
      else        kz[iz] = dkz*( double )(iz-nz);
    }
    
    return HELMHOLZ_SUCCESS;
}


int finalize_helmholz_decomp()
{
    
    free(_Acr_x); free(_Air_x); free(_zAk_x); free(_zAck_x); free(_zAik_x);
    free(_Acr_y); free(_Air_y); free(_zAk_y); free(_zAck_y); free(_zAik_y);
    free(_Acr_z); free(_Air_z); free(_zAk_z); free(_zAck_z); free(_zAik_z);
    
    free(kx); free(ky); free(kz);
    free(ixyz2ix); free(ixyz2iy); free(ixyz2iz);
    
    fftw_destroy_plan(plan_f); fftw_free(zin);
    fftw_destroy_plan(plan_b); fftw_free(zout);
    
    flag_initialized = 0;
    
    return HELMHOLZ_SUCCESS;
}


int helmholz_decomp(double*  Ar_x, double*  Ar_y, double*  Ar_z,
                    double* Acr_x, double* Acr_y, double* Acr_z,
                    double* Air_x, double* Air_y, double* Air_z,
                    double complex**  Ak_x, double complex**  Ak_y, double complex**  Ak_z,
                    double*  A0_x, double*  A0_y, double*  A0_z,
                    const int nx, const int ny, const int nz)
{
    const int nxyz = nx*ny*nz;
    if (flag_initialized == 0) init_helmholz_decomp(nx,ny,nz);
    
    // Fourier transform: A(r) -> A(k)
    for(int i=0;i<nxyz;i++) zin[i]=Ar_x[i]; fftw_execute(plan_f); for(int i=0;i<nxyz;i++) _zAk_x[i]=zout[i];
    for(int i=0;i<nxyz;i++) zin[i]=Ar_y[i]; fftw_execute(plan_f); for(int i=0;i<nxyz;i++) _zAk_y[i]=zout[i];
    for(int i=0;i<nxyz;i++) zin[i]=Ar_z[i]; fftw_execute(plan_f); for(int i=0;i<nxyz;i++) _zAk_z[i]=zout[i];
    
//     printf("FFTW FORWARD MADE\n");
//     printf("%lf\n",Acr_x[0]);
    
    // Decomposition: A(r) = A0 + Ac(r) + Ai(r)
    // ========================================= Compressible part ===================================================
    for(int k=0;k<nxyz;k++)// x
    {
        int ix=ixyz2ix[k];
        int iy=ixyz2iy[k];
        int iz=ixyz2iz[k];
        double kk=kx[ix]*kx[ix]+ky[iy]*ky[iy]+kz[iz]*kz[iz];
        if   (kk!=0.0) zin[k] = (kx[ix]*_zAk_x[k] + ky[iy]*_zAk_y[k] + kz[iz]*_zAk_z[k])/kk*kx[ix];
        else           zin[k] = 0.0;
    }
    fftw_execute(plan_b); for(int i=0;i<nxyz;i++) Acr_x[i]=zout[i]/nxyz;
    for(int k=0;k<nxyz;k++)// y
    {
	int ix=ixyz2ix[k];
        int iy=ixyz2iy[k];
        int iz=ixyz2iz[k];
        double kk=kx[ix]*kx[ix]+ky[iy]*ky[iy]+kz[iz]*kz[iz];
	if(kk!=0.0) zin[k] = (kx[ix]*_zAk_x[k] + ky[iy]*_zAk_y[k] + kz[iz]*_zAk_z[k])/kk*ky[iy];
	else        zin[k] = 0.0;
    }
    fftw_execute(plan_b); for(int i=0;i<nxyz;i++) Acr_y[i]=zout[i]/nxyz;
    for(int k=0;k<nxyz;k++)// z
    {
	int ix=ixyz2ix[k];
        int iy=ixyz2iy[k];
        int iz=ixyz2iz[k];
        double kk=kx[ix]*kx[ix]+ky[iy]*ky[iy]+kz[iz]*kz[iz];
	if(kk!=0.0) zin[k] = (kx[ix]*_zAk_x[k] + ky[iy]*_zAk_y[k] + kz[iz]*_zAk_z[k])/kk*kz[iz];
	else      zin[k] = 0.0;
    }
    fftw_execute(plan_b); for(int i=0;i<nxyz;i++) Acr_z[i]=zout[i]/nxyz;
    
//     printf("COMPRESSIBLE PART COMPUTED\n");
    
    
    // ========================================= Incompressible part =================================================
    for(int k=0;k<nxyz;k++)// x
    {
	int ix=ixyz2ix[k];
        int iy=ixyz2iy[k];
        int iz=ixyz2iz[k];
        double kk=kx[ix]*kx[ix]+ky[iy]*ky[iy]+kz[iz]*kz[iz];
        if(kk!=0.0) zin[k] = _zAk_x[k] - (kx[ix]*_zAk_x[k] + ky[iy]*_zAk_y[k] + kz[iz]*_zAk_z[k])/kk*kx[ix];
	else      zin[k] = 0.0;
    }
    fftw_execute(plan_b); for(int i=0;i<nxyz;i++) Air_x[i]=zout[i]/nxyz;
    for(int k=0;k<nxyz;k++)// y
    {
	int ix=ixyz2ix[k];
        int iy=ixyz2iy[k];
        int iz=ixyz2iz[k];
        double kk=kx[ix]*kx[ix]+ky[iy]*ky[iy]+kz[iz]*kz[iz];
	if(kk!=0.0) zin[k] = _zAk_y[k] - (kx[ix]*_zAk_x[k] + ky[iy]*_zAk_y[k] + kz[iz]*_zAk_z[k])/kk*ky[iy];
	else      zin[k] = 0.0;
    }
    fftw_execute(plan_b); for(int i=0;i<nxyz;i++) Air_y[i]=zout[i]/nxyz;
    for(int k=0;k<nxyz;k++)// z
    {
	int ix=ixyz2ix[k];
        int iy=ixyz2iy[k];
        int iz=ixyz2iz[k];
        double kk=kx[ix]*kx[ix]+ky[iy]*ky[iy]+kz[iz]*kz[iz];
        if(kk!=0.0) zin[k] = _zAk_z[k] - (kx[ix]*_zAk_x[k] + ky[iy]*_zAk_y[k] + kz[iz]*_zAk_z[k])/kk*kz[iz];
	else      zin[k] = 0.0;
    }
    fftw_execute(plan_b); for(int i=0;i<nxyz;i++) Air_z[i]=zout[i]/nxyz;
    
//     printf("INCOMPRESSIBLE PART COMPUTED\n");
    
    
    // k=0 part (constant in space)
    //double A0_x, A0_y, A0_z;
    *A0_x = Ar_x[0] - Acr_x[0] - Air_x[0];
    *A0_y = Ar_y[0] - Acr_y[0] - Air_y[0];
    *A0_z = Ar_z[0] - Acr_z[0] - Air_z[0];
    
//     printf("CONSTANT PART COMPUTED\n");
    
    
    
    // if pointer is not null get the address of an array storying A_tot(k)
    if (Ak_x) *Ak_x = _zAk_x;
    if (Ak_y) *Ak_y = _zAk_y;
    if (Ak_z) *Ak_z = _zAk_z;
    
    return HELMHOLZ_SUCCESS;
}


int fourier_transform_vector(double *A_x, double *A_y, double *A_z,
                             double complex *Ak_x, double complex *Ak_y, double complex *Ak_z,
                             const int nx, const int ny, const int nz)
{
    const int nxyz = nx*ny*nz;
    
    // Fourier transforms
    // A(r) -> A(k)
    for(int i=0;i<nxyz;i++) zin[i]=A_x[i]; fftw_execute(plan_f); for(int i=0;i<nxyz;i++) Ak_x[i] = zout[i];
    for(int i=0;i<nxyz;i++) zin[i]=A_y[i]; fftw_execute(plan_f); for(int i=0;i<nxyz;i++) Ak_y[i] = zout[i];
    for(int i=0;i<nxyz;i++) zin[i]=A_z[i]; fftw_execute(plan_f); for(int i=0;i<nxyz;i++) Ak_z[i] = zout[i];
    
    return HELMHOLZ_SUCCESS;
}


int get_kvector(double* kx_extern, double* ky_extern, double* kz_extern,
    const int nx, const int ny, const int nz)
{
    for (int ix=0; ix < nx; ix++) kx_extern[ix] = kx[ix];
    for (int iy=0; iy < ny; iy++) ky_extern[iy] = ky[iy];
    for (int iz=0; iz < nz; iz++) kz_extern[iz] = kz[iz];
    
    return HELMHOLZ_SUCCESS;
}


#endif /* __HELMHOLZ_DECOMP__ */