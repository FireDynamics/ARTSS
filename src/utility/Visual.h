/// \file       Visual.h
/// \brief      Writes all the data needed for post-processing to csv and to a
///             single binary file in VTK legacy format in binary
/// \details    Prepares the arrays in a correct format and writes the
///             structured grid and its zonal variables
///             (\a u,\a v,\a p,\a d,\a T)
/// \date       Jul 12, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef VISUAL_H_
#define VISUAL_H_

#include <string>

#include "../analysis/Solution.h"
#include "../interfaces/SolverI.h"

class Visual {
 public:
    Visual();
    void Visualize(SolverI* solver, real t, const char *fname);

 private:
    void vtkWriteStep(
            const char *fname, int n, read_ptr u, read_ptr v, read_ptr w,
            read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s,
            read_ptr nu_t, read_ptr S_T, read_ptr ua, read_ptr va, read_ptr wa,
            read_ptr pa, read_ptr Ta);

    std::string RemoveExtension(const std::string& filename);
    void vtkPrepareAndWrite(
            const char *fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p,
            read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t,
            read_ptr S_T);
    void vtkPrepareAndWrite(
            const char *fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p,
            read_ptr T, read_ptr s);
    void csvPrepareAndWrite(
            std::string fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p,
            read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t,
            read_ptr S_T);
    void csvPrepareAndWrite(
            std::string fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p,
            read_ptr T, read_ptr s);

    int m_nx = 0, m_ny = 0, m_nz = 0;  // integer due to visit_writer boundary!
    int m_Nx = 0, m_Ny = 0, m_Nz = 0, m_size = 0;
    real m_x1, m_y1, m_z1;
    real m_X1, m_Y1, m_Z1;
    real m_dx, m_dy, m_dz;
    int m_Nt;

    Solution m_solution;
};

#endif /* VISUAL_H_ */
