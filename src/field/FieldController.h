/// \file       FieldController.h
/// \brief      manages everything regarding any field object
/// \date       Aug 12, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_FIELD_FIELDCONTROLLER_H_
#define ARTSS_FIELD_FIELDCONTROLLER_H_


#include "../Domain.h"
#include "../interfaces/ISource.h"
#include "../utility/GlobalMacrosTypes.h"
#include "Field.h"

class FieldController {
 public:
    explicit FieldController(Domain const &domain);

    void set_up_boundary();
    void update_data(bool sync);

    // Getter
    Field &get_field_u() { return field_u; }
    Field &get_field_v() { return field_v; }
    Field &get_field_w() { return field_w; }
    real* get_field_u_data() const { return field_u.data; }
    real* get_field_v_data() const { return field_v.data; }
    real* get_field_w_data() const { return field_w.data; }

    Field &get_field_u0() { return field_u0; }
    Field &get_field_v0() { return field_v0; }
    Field &get_field_w0() { return field_w0; }
    real *get_field_u0_data() const { return field_u0.data; }
    real *get_field_v0_data() const { return field_v0.data; }
    real *get_field_w0_data() const { return field_w0.data; }

    Field &get_field_u_tmp() { return field_u_tmp; }
    Field &get_field_v_tmp() { return field_v_tmp; }
    Field &get_field_w_tmp() { return field_w_tmp; }
    real *get_field_u_tmp_data() const { return field_u_tmp.data; }
    real *get_field_v_tmp_data() const { return field_v_tmp.data; }
    real *get_field_w_tmp_data() const { return field_w_tmp.data; }

    Field &get_field_p() { return field_p; }
    Field &get_field_p0() { return field_p0; }
    Field &get_field_rhs() { return field_rhs; }
    real *get_field_p_data() const { return field_p.data; }
    real *get_field_p0_data() const { return field_p0.data; }
    real *get_field_rhs_data() const { return field_rhs.data; }

    Field &get_field_T() { return field_T; }
    Field &get_field_T0() { return field_T0; }
    Field &get_field_T_tmp() { return field_T_tmp; }
    Field &get_field_T_ambient() { return field_T_ambient; }
    real *get_field_T_data() const { return field_T.data; }
    real *get_field_T0_data() const { return field_T0.data; }
    real *get_field_T_tmp_data() const { return field_T_tmp.data; }

    Field &get_field_concentration() { return field_concentration; }
    Field &get_field_concentration0() { return field_concentration0; }
    Field &get_field_concentration_tmp() { return field_concentration_tmp; }
    real * get_field_concentration_data() const { return field_concentration.data; }
    real * get_field_concentration0_data() const { return field_concentration0.data; }
    real * get_field_concentration_tmp_data() const { return field_concentration_tmp.data; }

    Field &get_field_gamma_t() { return field_gamma_t; }
    Field &get_field_kappa_t() { return field_kappa_t; }

    Field &get_field_sight() { return sight; }
    real *get_field_sight_data() const { return sight.data; }

    Field &get_field_nu_t() { return field_nu_t; }
    real *get_field_nu_t_data() const { return field_nu_t.data; }

    Field &get_field_source_T() { return field_source_T; }
    real *get_field_source_T_data() const { return field_source_T.data; }

    real *get_field_force_x_data() const { return field_force_x.data; }
    real *get_field_force_y_data() const { return field_force_y.data; }
    real *get_field_force_z_data() const { return field_force_z.data; }

    real *get_field_kappa_data() const { return field_kappa_t.data; }
    real *get_field_gamma_data() const { return field_gamma_t.data; }

    Field &get_field_source_concentration() { return field_source_concentration; }
    real *get_field_source_concentration_data() const { return field_source_concentration.data; }

    Field &get_field_force_x() { return field_force_x; }
    Field &get_field_force_y() { return field_force_y; }
    Field &get_field_force_z() { return field_force_z; }

    Field field_u, field_v, field_w;          // velocities
    Field field_u0, field_v0, field_w0;
    Field field_u_tmp, field_v_tmp, field_w_tmp;
    Field field_nu_t;
    Field field_kappa_t;
    Field field_gamma_t;

    Field field_p;                  // pressure
    Field field_p0;
    Field field_rhs;

    Field field_T;                  // temperature
    Field field_T0;
    Field field_T_tmp;
    Field field_T_ambient;

    Field field_concentration;      // smoke concentration
    Field field_concentration0;
    Field field_concentration_tmp;

    Field field_force_x, field_force_y, field_force_z;    // sources
    Field field_source_T;                // temperature
    Field field_source_concentration;    // smoke concentration

    Field sight;

    static void couple_vector(
            Field const &a, Field &a0, Field &a_tmp,
            Field const &b, Field &b0, Field &b_tmp,
            Field const &c, Field &c0, Field &c_tmp, bool sync);
    static void couple_scalar(Field const &a, Field &a0, Field &a_tmp, bool sync);

    void update_device();
    void update_host();
};

#endif /* ARTSS_FIELD_FIELDCONTROLLER_H_ */

