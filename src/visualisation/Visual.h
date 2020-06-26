/// \file       Visual.h
/// \brief      coordinator of different visualisation methods
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_VISUALISATION_VISUAL_H_
#define ARTSS_VISUALISATION_VISUAL_H_

#include "../analysis/Solution.h"
#include "../interfaces/ISolver.h"

class Visual {
public:
    Visual();

    void visualise(ISolver *solver, real t);

    static void initialiseGrid(float *x_coords, float *y_coords, float *z_coords, int Nx, int Ny, int Nz, real dx, real dy, real dz);

    static void prepareFields(read_ptr *fields, float **vars, int size);

private:
    static std::string RemoveExtension(const std::string &filename);

    std::string m_filename;
    Solution m_solution;
    bool m_save_csv = false;
    bool m_save_vtk = false;

    static std::string createFilename(std::string filename, real t, bool analytical);

    bool m_has_analytical_solution = false;
};

#endif /* ARTSS_VISUALISATION_VISUAL_H_ */
