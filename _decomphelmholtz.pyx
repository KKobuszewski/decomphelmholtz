# distutils: language = c++
#cython: boundscheck=False, nonecheck=False, cdivision=True
from __future__ import print_function, division

cimport numpy as np
import numpy as np


cimport cython

from libc.stdio cimport printf


cdef extern from "decomphelmholtz.h" nogil:
    int init_helmholz_decomp(const int nx, const int ny, const int nz);
    int finalize_helmholz_decomp();
    int helmholz_decomp(double*  Ar_x, double*  Ar_y, double*  Ar_z,
                        double* Acr_x, double* Acr_y, double* Acr_z,
                        double* Air_x, double* Air_y, double* Air_z,
                        double complex**  Ak_x, double complex**  Ak_y, double complex**  Ak_z,
                        double*  A0_x, double*  A0_y, double*  A0_z,
                        const int nx, const int ny, const int nz);
    int fourier_transform_vector(double *A_x, double *A_y, double *A_z,
                             double complex *Ak_x, double complex *Ak_y, double complex *Ak_z,
                             const int nx, const int ny, const int nz);
    int get_kvector(double* kx_extern, double* ky_extern, double* ky_extern,
                    const int nx, const int ny, const int nz)


cpdef helmholtz(np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Ax,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Ay,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Az,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Acx,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Acy,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Acz,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Aix,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Aiy,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Aiz):
    
    cdef int nx = <int> Ax.shape[0]
    cdef int ny = <int> Ax.shape[1]
    cdef int nz = <int> Ax.shape[2]
    """
    if Acx is None:
        Acx = np.empty((nx,ny,nz),dtype=np.float64_t,order='c')
    if Acy is None:
        Acy = np.empty((nx,ny,nz),dtype=np.float64_t,order='c')
    if Acz is None:
        Acz = np.empty((nx,ny,nz),dtype=np.float64_t,order='c')
    if Aix is None:
        Aix = np.empty((nx,ny,nz),dtype=np.float64_t,order='c')
    if Aiy is None:
        Aiy = np.empty((nx,ny,nz),dtype=np.float64_t,order='c')
    if Aiz is None:
        Acz = np.empty((nx,ny,nz),dtype=np.float64_t,order='c')
    """
    
    cdef double A0x = 0.0
    cdef double A0y = 0.0
    cdef double A0z = 0.0
    
    err = helmholz_decomp( &Ax[0,0,0], &Ay[0,0,0], &Az[0,0,0],
                          &Acx[0,0,0],&Acy[0,0,0],&Acz[0,0,0],
                          &Aix[0,0,0],&Aiy[0,0,0],&Aiz[0,0,0],
                          NULL, NULL, NULL, &A0x, &A0y, &A0z,
                          nx, ny, nz)
    
    return Acx, Acy, Acz, Aix, Aiy, Aiz, A0x, A0y, A0z


cpdef energy_sectrum_3d(np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Ax,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Ay,
                np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Az):
    
    cdef int nx = <int> Ax.shape[0]
    cdef int ny = <int> Ax.shape[1]
    cdef int nz = <int> Ax.shape[2]
    
    cdef np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Acx = np.empty((nx,ny,nz),dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Acy = np.empty((nx,ny,nz),dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Acz = np.empty((nx,ny,nz),dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Aix = np.empty((nx,ny,nz),dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Aiy = np.empty((nx,ny,nz),dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=3,negative_indices=False,mode='c'] Aiz = np.empty((nx,ny,nz),dtype=np.float64,order='c')
    
    Acx, Acy, Acz, Aix, Aiy, Aiz, A0x, A0y, A0z = helmholtz(Ax, Ay, Az, Acx, Acy, Acz, Aix, Aiy, Aiz)
    
    
    cdef np.ndarray[np.complex128_t,ndim=3,negative_indices=False,mode='c'] Akcx = np.empty((nx,ny,nz),dtype=np.complex128,order='c')
    cdef np.ndarray[np.complex128_t,ndim=3,negative_indices=False,mode='c'] Akcy = np.empty((nx,ny,nz),dtype=np.complex128,order='c')
    cdef np.ndarray[np.complex128_t,ndim=3,negative_indices=False,mode='c'] Akcz = np.empty((nx,ny,nz),dtype=np.complex128,order='c')
    cdef np.ndarray[np.complex128_t,ndim=3,negative_indices=False,mode='c'] Akix = np.empty((nx,ny,nz),dtype=np.complex128,order='c')
    cdef np.ndarray[np.complex128_t,ndim=3,negative_indices=False,mode='c'] Akiy = np.empty((nx,ny,nz),dtype=np.complex128,order='c')
    cdef np.ndarray[np.complex128_t,ndim=3,negative_indices=False,mode='c'] Akiz = np.empty((nx,ny,nz),dtype=np.complex128,order='c')
    
    err = fourier_transform_vector(&Acx[0,0,0],&Acy[0,0,0],&Acz[0,0,0],&Akcx[0,0,0],&Akcy[0,0,0],&Akcz[0,0,0], nx, ny, nz)
    err = fourier_transform_vector(&Aix[0,0,0],&Aiy[0,0,0],&Aiz[0,0,0],&Akix[0,0,0],&Akiy[0,0,0],&Akiz[0,0,0], nx, ny, nz)
    
    ecden_k = 0.5 * ( np.conj(Akcx)*Akcx + np.conj(Akcy)*Akcy + np.conj(Akcz)*Akcz )
    eiden_k = 0.5 * ( np.conj(Akix)*Akix + np.conj(Akcy)*Akiy + np.conj(Akiz)*Akiz )
    
    #Etot = 0.5*(Ax*Ax + Ay*Ay + Az*Az).sum() # * dx
    
    cdef np.ndarray[np.float64_t,ndim=1,negative_indices=False,mode='c'] kx = np.empty(nx,dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=1,negative_indices=False,mode='c'] ky = np.empty(ny,dtype=np.float64,order='c')
    cdef np.ndarray[np.float64_t,ndim=1,negative_indices=False,mode='c'] kz = np.empty(nz,dtype=np.float64,order='c')
    err = get_kvector(&kx[0],&ky[0],&kz[0], nx, ny, nz)
    
    return ecden_k, eiden_k, np.meshgrid(kx,ky,kz,indexing='ij')
