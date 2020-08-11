/// \file       GlobalMacrosTypes.h
/// \brief      Defines macros and typedefs
/// \details    This Header defines macros for indexing (\c IX(i,j)), swapping
///             values (\c SWAP(a,b)) and nested for-loops (\c FOR_EACH_CELL,
///             \c FOR_EACH_INNER_CELL)
/// \date       Jun 15, 2015
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_UTILITY_GLOBALMACROSTYPES_H_
#define ARTSS_UTILITY_GLOBALMACROSTYPES_H_

#include <cstdlib>

//================================= Typedefs =================================
// ***************************************************************************
/// \def real
/// \brief Defines the data type the solver calculates with
///
/// \def read_ptr
/// \brief Defines a constant pointer to a constant real, __restrict__
///        for GPU version
///
/// \def write_ptr
/// \brief Defines a constant pointer to a constant real, __restrict__
///        for GPU version
///
/// \def return_ptr
/// \brief For returning data by class functions
// ***************************************************************************


typedef double real;  // data type for solver (float, double, ...)

// looks like PGI uses __restrict
// https://www.auburn.edu/cosam/departments/physics/department/comp-resources/files/pgi/pgicdkrn.pdf
// p. 14
// nontheless this only a c99 keyword and should not be used
#ifdef __PGI
    typedef const real* __restrict const read_ptr;  // readable for GPU version
    typedef real* __restrict const write_ptr;       // writable ptr GPU version
    typedef const real* __restrict return_ptr;      // returning by class
// const pointer is superfluous, because non-class type
// return values are not modifiable anyway.

#else
    typedef const real* const read_ptr;  // readable
    typedef real* const write_ptr;       // writable
    typedef const real* return_ptr;      // for returning by class functions.
#endif

typedef const real* const aliased_read_ptr;


//================================== Macros ==================================
// ****************************************************************************
/// \def IX(i,j,k,Nx,Ny)
/// \brief Computes the one-dimensional row-major index
///        of a three-dimensional array \a (i,j,k)
///
/// \def xi(i,x,dx)
/// \brief Physical x-coordinates at midpoints calculated by index \a i=0..Nx-1
///
/// \def yj(j,y,dy)
/// \brief Physical y-coordinates at midpoints calculated by index \a j=0..Ny-1
///
/// \def zk(k,z,dz)
/// \brief Physical z-coordinates at midpoints calculated by index \a k=0..Nz-1
///
/// \def tn(n,dt)
/// \brief Simulation time calculated by time steps \a n = 0...Nt
// ***************************************************************************************

#define IX(i,j,k,Nx,Ny) ((i) + (Nx)*(j) + (Nx)*(Ny)*(k)) // row-major index for one dimensional arrays (i=0..Nx-1 columns, j=0..Ny-1 rows, k = 0...Nz-1)
#define getCoordinateI(idx,Nx,Ny,j,k) (idx - k * Nx * Ny - j * Nx)
#define getCoordinateJ(idx,Nx,Ny,k) ((idx - k * Nx * Ny) / Nx)
#define getCoordinateK(idx,Nx,Ny)(idx / (Nx * Ny))


#define xi(i,x,dx) ((x) + ((i)-0.5)*(dx)) // physical xcoords at midpoints calculated by index i=0..Nx-1
#define yj(j,y,dy) ((y) + ((j)-0.5)*(dy)) // physical ycoords at midpoints calculated by index j=0..Ny-1
#define zk(k,z,dz) ((z) + ((k)-0.5)*(dz)) // physical zcoords at midpoints calculated by index k=0..Nz-1
#define tn(n, dt) ((n)*(dt))			  // simulation time calculated by time steps n = 0...Nt

#endif /* ARTSS_UTILITY_GLOBALMACROSTYPES_H_ */
