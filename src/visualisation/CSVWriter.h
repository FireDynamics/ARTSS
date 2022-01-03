/// \file       CSVWriter.h
/// \brief      class to write out csv files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_VISUALISATION_CSVWRITER_H
#define ARTSS_VISUALISATION_CSVWRITER_H

#include <string>
#include <vector>

#include "Visual.h"
#include "../domain/DomainData.h"
#include "../analysis/Solution.h"
#include "../field/FieldController.h"

class CSVWriter {
 public:
    static void write_numerical(const FieldController &field_controller, const std::string &filename);
    static void write_analytical(const Solution &solution, const std::string &filename);

 private:
    static void csv_prepare_and_write(const std::string &filename,
                                      return_ptr u, return_ptr v, return_ptr w,
                                      return_ptr p,
                                      return_ptr div,
                                      return_ptr T,
                                      return_ptr C,
                                      return_ptr sight,
                                      return_ptr nu_t,
                                      return_ptr source_T);
    static void csv_prepare_and_write(const std::string &filename,
                                      read_ptr u, read_ptr v, read_ptr w,
                                      read_ptr p,
                                      read_ptr T);
    static void csv_write(const std::string &filename, return_ptr *vars, int size_vars, const std::vector<std::string> &var_names);
};

#endif /* ARTSS_VISUALISATION_CSVWRITER_H */
