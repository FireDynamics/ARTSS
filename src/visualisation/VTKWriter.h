/// \file       VTKWriter.h
/// \brief      class to write out vtk files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_VISUALISATION_VTKWRITER_H
#define ARTSS_VISUALISATION_VTKWRITER_H

#include "../interfaces/ISolver.h"
#include "../analysis/Solution.h"

class VTKWriter {
public:
    static void write_numerical(ISolver *solver, const std::string& filename);
    static void write_analytical(Solution *solution, const std::string& filename);

private:
    static void vtkPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T);
    static void vtkPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T);
};

#endif /* ARTSS_VISUALISATION_VTKWRITER_H */
