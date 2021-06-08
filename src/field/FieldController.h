/// \file       FieldController.h
/// \brief      manages everything regarding any field object
/// \date       Aug 12, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_FIELD_FIELDCONTROLLER_H_
#define ARTSS_FIELD_FIELDCONTROLLER_H_


#include "../interfaces/ISource.h"
#include "../utility/GlobalMacrosTypes.h"
#include "Field.h"

class FieldController {
 public:
    FieldController();
    ~FieldController();

    void set_up_boundary();
    void set_up_temporary_fields();
    void update_data(bool sync);
    void replace_data(Field *u, Field *v, Field *w, Field *p, Field *T, Field *C);

    // Getter
    real* get_field_u_data() const { return field_u->data; }
    real* get_field_v_data() const { return field_v->data; }
    real* get_field_w_data() const { return field_w->data; }

    real* get_field_u0_data() const { return field_u0->data; }
    real* get_field_v0_data() const { return field_v0->data; }
    real* get_field_w0_data() const { return field_w0->data; }

    real* get_field_u_tmp_data() const { return field_u_tmp->data; }
    real* get_field_v_tmp_data() const { return field_v_tmp->data; }
    real* get_field_w_tmp_data() const { return field_w_tmp->data; }

    real* get_field_p_data() const { return field_p->data; }
    real* get_field_p0_data() const { return field_p0->data; }
    real* get_field_rhs_data() const { return field_rhs->data; }

    real* get_field_T_data() const { return field_T->data; }
    real* get_field_T0_data() const { return field_T0->data; }
    real* get_field_T_tmp_data() const { return field_T_tmp->data; }
    real* get_field_T_ambient_data() const { return field_T_ambient->data; }

    real* get_field_concentration_data() const { return field_concentration->data; }
    real* get_field_concentration0_data() const { return field_concentration0->data; }
    real* get_field_concentration_tmp_data() const { return field_concentration_tmp->data; }

    real* get_field_sight_data() const { return sight->data; }

    real* get_field_nu_t_data() const { return field_nu_t->data; }

    real* get_field_source_T_data() const { return field_source_T->data; }

    Field *field_u, *field_v, *field_w;          // velocities
    Field *field_u0, *field_v0, *field_w0;
    Field *field_u_tmp, *field_v_tmp, *field_w_tmp;
    Field *field_nu_t;
    Field *field_kappa_t;
    Field *field_gamma_t;

    Field *field_p;                  // pressure
    Field *field_p0;
    Field *field_rhs;

    Field *field_T;                  // temperature
    Field *field_T0;
    Field *field_T_tmp;
    Field *field_T_ambient;

    Field *field_concentration;      // smoke concentration
    Field *field_concentration0;
    Field *field_concentration_tmp;

    Field *field_force_x, *field_force_y, *field_force_z;    // sources
    Field *field_source_T;                // temperature
    Field *field_source_concentration;    // smoke concentration

    Field *sight;

    static void couple_vector(const Field *a, Field *a0, Field *a_tmp, const Field *b, Field *b0, Field *b_tmp, const Field *c, Field *c0, Field *c_tmp, bool sync);
    static void couple_scalar(const Field *a, Field *a0, Field *a_tmp, bool sync);

    void update_device();
    void update_host();
};

#endif /* ARTSS_FIELD_FIELDCONTROLLER_H_ */

