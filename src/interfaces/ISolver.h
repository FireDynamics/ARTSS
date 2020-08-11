/// \file       ISolver.h
/// \brief      Interface holds solvers for solving governing equations
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_ISOLVER_H_
#define ARTSS_INTERFACES_ISOLVER_H_

#include <string>

#include "../Field.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "ISource.h"


struct SolverTypes {
    inline static const std::string AdvectionSolver = "AdvectionSolver";
    inline static const std::string AdvectionDiffusionSolver = "AdvectionDiffusionSolver";
    inline static const std::string DiffusionSolver = "DiffusionSolver";
    inline static const std::string DiffusionTurbSolver = "DiffusionTurbSolver";
    inline static const std::string NSSolver = "NSSolver";
    inline static const std::string NSTempSolver = "NSTempSolver";
    inline static const std::string NSTempConSolver = "NSTempConSolver";
    inline static const std::string NSTempTurbSolver = "NSTempTurbSolver";
    inline static const std::string NSTempTurbConSolver = "NSTempTurbConSolver";
    inline static const std::string NSTurbSolver = "NSTurbSolver";
    inline static const std::string PressureSolver = "PressureSolver";
};

class ISolver {
public:
    ISolver();
    virtual ~ISolver();

    virtual void do_step(real t, bool sync) = 0;

    void set_up_boundary(bool sync = true);

    void update_data(bool sync = true);

    void update_sources(real t, bool sync);

    // Getter
    return_ptr get_u() const { return u->data; }
    return_ptr get_v() const { return v->data; }
    return_ptr get_w() const { return w->data; }

    return_ptr get_u0() const { return u0->data; }
    return_ptr get_v0() const { return v0->data; }
    return_ptr get_w0() const { return w0->data; }

    return_ptr get_u_tmp() const { return u_tmp->data; }
    return_ptr get_v_tmp() const { return v_tmp->data; }
    return_ptr get_w_tmp() const { return w_tmp->data; }

    return_ptr get_p() const { return p->data; }
    return_ptr get_p0() const { return p0->data; }
    return_ptr get_rhs() const { return rhs->data; }

    return_ptr get_T() const { return T->data; }
    return_ptr get_T0() const { return T0->data; }
    return_ptr get_T_tmp() const { return T_tmp->data; }

    return_ptr get_concentration() const { return concentration->data; }
    return_ptr get_concentration0() const { return concentration0->data; }
    return_ptr get_concentration_tmp() const { return concentration_tmp->data; }

    return_ptr get_sight() const { return sight->data; }

    return_ptr get_nu_t() const { return nu_t->data; }

    return_ptr get_S_T() const { return S_T->data; }

    Field *u, *v, *w;          // velocities
    Field *u0, *v0, *w0;
    Field *u_tmp, *v_tmp, *w_tmp;
    Field *nu_t;
    Field *kappa_t;
    Field *gamma_t;

    Field *p;                  // pressure
    Field *p0;
    Field *rhs;

    Field *T;                  // temperature
    Field *T0;
    Field *T_tmp;
    Field *T_ambient;

    Field *concentration;      // smoke concentration
    Field *concentration0;
    Field *concentration_tmp;

    Field *f_x, *f_y, *f_z;    // sources
    Field *S_T;                // temperature
    Field *S_concentration;    // smoke concentration

    Field *sight;

    ISource *source_temperature;
    ISource *source_velocity;
    ISource *source_concentration;
protected:
    void set_up();

    void init();

    static void init_c(Field *out, real in);

    static void init_f(Field *out, const real *in);

    static void couple_vector(const Field *a, Field *a0, Field *a_tmp, const Field *b, Field *b0, Field *b_tmp, const Field *c, Field *c0, Field *c_tmp, bool sync);

    static void couple_scalar(const Field *a, Field *a0, Field *a_tmp, bool sync);

private:
    std::string m_string_solver;
    void force_source();
    void temperature_source();
    void momentum_source();
    static void CallRandom(Field* out);
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_INTERFACES_ISOLVER_H_ */
