/// \file       FieldWriter.h
/// \brief      
/// \date       Apr 18, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_VISUAL_FIELDWRITER_H
#define ARTSS_VISUAL_FIELDWRITER_H

#include <string>
#include "../field/FieldController.h"

class FieldWriter {
 public:
    FieldWriter(FieldController &field_controller, real dt);
    void write_out(real t_cur);

 private:
    std::string create_header();

    std::string m_filename = "visualisation.dat";
    FieldController &m_field_controller;
    real m_dt;
};


#endif /* ARTSS_VISUAL_FIELDWRITER_H */
