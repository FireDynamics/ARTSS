/// \file       VTKWriter.h
/// \brief      class to write out vtk files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_VISUALISATION_VTKWRITER_H
#define ARTSS_VISUALISATION_VTKWRITER_H

#include "../analysis/Solution.h"
#include "../field/FieldController.h"

class VTKWriter {
public:
    static void write_numerical(const FieldController &field_controller, const std::string &filename);
    static void write_analytical(const Solution &solution, const std::string &filename);
    static void write_field(const Field &field, const std::string &filename, const std::string &var_name);
    static void write_numerical_debug(const FieldController &field_controller, const std::string &filename);

private:
    static void vtk_prepare_and_write_debug(const char *filename, read_ptr *data,
                                            int size_vars, const char * const *var_names,
                                            int *centering, int *var_dims);
    static void vtk_prepare_and_write(const char *filename,
                                      read_ptr u, read_ptr v, read_ptr w,
                                      read_ptr p,
                                      read_ptr div,
                                      read_ptr T,
                                      read_ptr C,
                                      read_ptr sight,
                                      read_ptr nu_t,
                                      read_ptr source_T);
    static void vtk_prepare_and_write(const char *filename,
                                      read_ptr u, read_ptr v, read_ptr w,
                                      read_ptr p,
                                      read_ptr T);
    static int vtk_counter;
};

#endif /* ARTSS_VISUALISATION_VTKWRITER_H */
