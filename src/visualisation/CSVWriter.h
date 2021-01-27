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
#include "../Domain.h"
#include "../analysis/Solution.h"
#include "../field/FieldController.h"

class CSVWriter {
public:
    static void write_numerical(const FieldController& field_controller, const std::string& filename);
    static void write_analytical(const Solution& solution, const std::string& filename);

private:
    static void csvPrepareAndWrite(const std::string& filename, real *u, real* v, real* w, real* p, real* div, real* T, real* C, real* s, real* nu_t, real* S_T);
    static void csvPrepareAndWrite(const std::string& filename, real *u, real* v, real* w, real* p, real* T);

    static void csv_write(const std::string& filename, real **vars, int size_vars, const std::vector<std::string>& var_names);
};

#endif /* ARTSS_VISUALISATION_CSVWRITER_H */
