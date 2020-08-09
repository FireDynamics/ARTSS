/// \file       CSVWriter.h
/// \brief      class to write out csv files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_VISUALISATION_CSVWRITER_H
#define ARTSS_VISUALISATION_CSVWRITER_H

#include <string>
#include "../interfaces/ISolver.h"
#include "../analysis/Solution.h"

class CSVWriter {
public:
    static void write_numerical(ISolver *solver, const std::string& filename);
    static void write_analytical(Solution *solution, const std::string& filename);

private:
    static void csvPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T);
    static void csvPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T);

    static void csv_write(const char *filename, float **vars, int size_vars, const char **var_names);
};

#endif /* ARTSS_VISUALISATION_CSVWRITER_H */
