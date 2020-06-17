/// \file       SolverI.h
/// \brief      Interface holds solvers for solving governing equations
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ATSS_INTERFACES_SOLVERI_H_
#define ATSS_INTERFACES_SOLVERI_H_

#include <string>

#include "../Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "ISource.h"

#ifndef PROFILING
#include <spdlog/logger.h>
#endif

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

    virtual void DoStep(real t, bool sync)=0;

    void SetUpBoundary(bool sync = true);
    void UpdateData(bool sync = true);
    void UpdateSources(real t, bool sync);

    // Getter
    return_ptr GetU() const {
        return u->data;
    }
    return_ptr GetV() const {
        return v->data;
    }
    return_ptr GetW() const {
        return w->data;
    }
    return_ptr GetU0() const {
        return u0->data;
    }
    return_ptr GetV0() const {
        return v0->data;
    }
    return_ptr GetW0() const {
        return w0->data;
    }
    return_ptr GetU_tmp() const {
        return u_tmp->data;
    }
    return_ptr GetV_tmp() const {
        return v_tmp->data;
    }
    return_ptr GetW_tmp() const {
        return w_tmp->data;
    }
    return_ptr GetP() const {
        return p->data;
    }
    return_ptr GetP0() const {
        return p0->data;
    }
    return_ptr GetRhs() const {
        return rhs->data;
    }
    return_ptr GetT() const {
        return T->data;
    }
    return_ptr GetT0() const {
        return T0->data;
    }
    return_ptr GetT_tmp() const {
        return T_tmp->data;
    }
    return_ptr GetC() const {
        return C->data;
    }
    return_ptr GetC0() const {
        return C0->data;
    }
    return_ptr GetC_tmp() const {
        return C_tmp->data;
    }

    return_ptr GetSight() const {
        return sight->data;
    }

    return_ptr GetNu_t() const {
            return nu_t->data;
    }

    return_ptr GetS_T() const {
        return S_T->data;
    }

    Field* u, *v, *w;               // velocities
    Field* u0, *v0, *w0;
    Field* u_tmp, *v_tmp, *w_tmp;
    Field* nu_t;
    Field* kappa_t;
    Field* gamma_t;

    Field* p;                             // pressure
    Field* p0;
    Field* rhs;

    Field* T;                             // temperature
    Field* T0;
    Field* T_tmp;
    Field* T_a;

    Field* C;                             // smoke concentration
    Field* C0;
    Field* C_tmp;

    Field* f_x, *f_y, *f_z;     // sources
    Field* S_T;                             // temperature
    Field* S_C;                             // smoke concentration

    Field* sight;

    ISource* sou_temp;
    ISource* sou_vel;
    ISource* sou_con;

protected:
    void SetUp();
    void Init();
    static void Init_c(Field* out, real in);
    static void Init_f(Field* out, const real* in);
    static void CoupleVector(const Field* a, Field* a0, Field* a_tmp, const Field* b, Field* b0, Field* b_tmp, const Field* c, Field* c0, Field* c_tmp, bool sync);
    static void CoupleScalar(const Field* a, Field* a0, Field* a_tmp, bool sync);

private:
    std::string m_string_solver;
    void ForceSource();
    void TemperatureSource();
    void MomentumSource();
#ifndef PROFILING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ATSS_INTERFACES_SOLVERI_H_ */
