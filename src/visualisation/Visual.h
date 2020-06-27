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
    explicit Visual(Solution* solution);

    void visualise(ISolver *solver, real t);

    static void initialise_grid(float *x_coords, float *y_coords, float *z_coords, int Nx, int Ny, int Nz, real dx, real dy, real dz);

    static void prepare_fields(read_ptr *fields, float **vars, int size);

private:
    static std::string remove_extension(const std::string &filename);

    std::string m_filename;
    Solution *m_solution;
    bool m_save_csv = false;
    int m_csv_plots = 0;
    bool m_save_vtk = false;
    int m_vtk_plots = 0;
    real m_dt;
    real m_t_end;

    static std::string create_filename(std::string filename, int counter, bool analytical);

    bool m_has_analytical_solution = false;
};

#endif /* ARTSS_VISUALISATION_VISUAL_H_ */
