/// \file 		Adaption.h
/// \brief 		Controll class for adaption
/// \date 		Nov 29, 2018
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADAPTION_ADAPTION_H_
#define ARTSS_ADAPTION_ADAPTION_H_


#include "../Field.h"
#include "../interfaces/AdaptionFunctionI.h"
#include "../utility/GlobalMacrosTypes.h"
#ifndef PROFILING
#include <spdlog/logger.h>
#endif

/* enum for different types of dynamic adaption:
 * NO = adaption impossible/no changes
 * UNKNOWN = adaption is possible but not determined
 * YES = adaption shall be carried out!
 */
enum class ADTypes : size_t {
    NO = 0, UNKNOWN = 1, YES = 2
};
enum VectorFieldsTypes : size_t {
    VEL_U = 0, // velocity in x direction
    VEL_V = 1, // velocity in y direction
    VEL_W = 2, // velocity in z direction
    VEL_U0 = 3,
    VEL_V0 = 4,
    VEL_W0 = 5,
    VEL_U_TMP = 6,
    VEL_V_TMP = 7,
    VEL_W_TMP = 8,
    PRESSURE = 9, // pressure
    PRESSURE0 = 10,
    RHS = 11,
    TEMPERATURE = 12, // temperature
    TEMPERATURE0 = 13,
    TEMPERATURE_TMP = 14,
    TEMPERATURE_A = 15,
    CONCENTRATION = 16, // smoke concentration
    CONCENTRATION0 = 17,
    CONCENTRATION_TMP = 18,
    FORCE_X = 19,
    FORCE_Y = 20,
    FORCE_Z = 21,
    SOURCE_T = 22,
    SOURCE_C = 23,
    NU_T = 24, // turb visc (physical diffusion)
    KAPPA_T = 25, // thermal diffusion
    GAMMA_T = 26
};

class Field;

class AdaptionFunctionI;

class Adaption {
public:
    Adaption(Field **fields);

    bool inline isDataExtractionEnabled() { return m_hasDataExtraction; };

    bool inline isDataExtractionBeforeEnabled() { return m_hasDataExtractionBefore; }

    bool inline isDataExtractionAfterEnabled() { return m_hasDataExtractionAfter; }

    bool inline isDataExtractionEndresultEnabled() { return m_hasDataExtractionEndresult; }

    bool inline isTimeMeasuringEnabled() { return m_hasTimeMeasuring; }

    bool inline isWriteFieldEnabled() { return m_hasWriteField; }

    bool inline isWriteRuntimeEnabled() { return m_hasWriteRuntime; }

    real inline getBeforeHeight() { return m_height_before; }

    std::string inline getBeforeName() { return m_filename + "_before.csv"; }

    real inline getAfterHeight() { return m_height_after; }

    std::string inline getAfterName() { return m_filename + "_after.csv"; }

    std::string inline getEndresultName() { return m_filename + "_endresult.csv"; }

    std::string inline getTimeMeasuringName() { return m_filename + "_time.csv"; }

    std::string inline getWriteFieldName(real t_cur) { return m_filename + "_" + std::to_string(t_cur) + ".csv"; }

    std::string inline getWriteRuntimeName() { return m_filename + "_runtime.csv"; }

    void run(real t_cur);

    void expandXDirection(long shift, bool start, size_t *arr_idxExpansion, size_t len_e);

    void expandYDirection(long shift, bool start, size_t *arr_idxExpansion, size_t len_e);

    void reduceXDirection(long shift, bool start, size_t *arr_idxReduction, size_t len_r);

    void reduceYDirection(long shift, bool start, size_t *arr_idxReduction, size_t len_r);

    bool adaptXDirection(const real *f, real checkValue, size_t noBufferCell, real threshold);

    bool adaptXDirection_serial(const real *f, real checkValue, size_t noBufferCell, real threshold);

    bool adaptYDirection(const real *f, real checkValue, size_t noBufferCell, real threshold);

    bool adaptYDirection_serial(const real *f, real checkValue, size_t noBufferCell, real threshold);

    Field **fields;

    void extractData(const std::string filename, real height, real time);

    void extractData(const std::string filename);

    long m_shift_x1, m_shift_x2, m_shift_y1, m_shift_y2, m_shift_z1, m_shift_z2;
private:
    std::shared_ptr<spdlog::logger> m_logger;
    bool isUpdateNecessary();

    void applyChanges();

    AdaptionFunctionI *func;
    bool m_dynamic, m_dynamic_end;
    bool m_reduction;
    bool m_hasDataExtraction;
    bool m_hasDataExtractionBefore=false;
    bool m_hasDataExtractionAfter=false;
    bool m_hasDataExtractionEndresult=false;
    bool m_hasTimeMeasuring=false;
    bool m_hasWriteField=false;
    bool m_hasWriteRuntime=false;

    real m_height_before = 0;
    real m_height_after = 0;

    size_t m_minimal;

    std::string m_filename;
};


#endif /* ARTSS_ADAPTION_ADAPTION_H_ */
