/// \file 		Adaption.h
/// \brief 		Controll class for adaption
/// \date 		Nov 29, 2018
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADAPTION_ADAPTION_H_
#define ARTSS_ADAPTION_ADAPTION_H_


#include "../Field.h"
#include "../interfaces/IAdaptionFunction.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../interfaces/ISolver.h"

/* enum for different types of dynamic adaption:
 * NO = adaption impossible/no changes
 * UNKNOWN = adaption is possible but not determined
 * YES = adaption shall be carried out!
 */
enum class ADTypes : size_t {
    NO = 0, UNKNOWN = 1, YES = 2
};

class Field;

class IAdaptionFunction;

class Adaption {
public:
    Adaption(ISolver *solver);

    bool inline isDataExtractionEnabled() {
        return m_hasDataExtraction;
    };

    bool inline isDataExtractionBeforeEnabled() {
        return m_hasDataExtractionBefore;
    }

    bool inline isDataExtractionAfterEnabled() {
        return m_hasDataExtractionAfter;
    }

    bool inline isDataExtractionEndresultEnabled() {
        return m_hasDataExtractionEndresult;
    }

    bool inline isTimeMeasuringEnabled() {
        return m_hasTimeMeasuring;
    }

    bool inline isWriteFieldEnabled() {
        return m_hasWriteField;
    }

    bool inline isWriteRuntimeEnabled() {
        return m_hasWriteRuntime;
    }

    real inline getBeforeHeight() {
        return m_height_before;
    }

    std::string inline getBeforeName() {
        return m_filename + "_before.csv";
    }

    real inline getAfterHeight() {
        return m_height_after;
    }

    std::string inline getAfterName() {
        return m_filename + "_after.csv";
    }

    std::string inline getEndresultName() {
        return m_filename + "_endresult.csv";
    }

    std::string inline getTimeMeasuringName() {
        return m_filename + "_time.csv";
    }

    std::string inline getWriteFieldName(real t_cur) {
        return m_filename + "_" + std::to_string(t_cur) + ".csv";
    }

    std::string inline getWriteRuntimeName() {
        return m_filename + "_runtime.csv";
    }

    void run(real t_cur);

    static void expandXDirection(long shift, bool start, size_t *arr_idxExpansion, size_t len_e);
    static void expandYDirection(long shift, bool start, size_t *arr_idxExpansion, size_t len_e);
    static void reduceXDirection(long shift, bool start, size_t *arr_idxReduction, size_t len_r);
    static void reduceYDirection(long shift, bool start, size_t *arr_idxReduction, size_t len_r);
    static bool adaptXDirection(const real *f, real checkValue, size_t noBufferCell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);
    static bool adaptXDirection_serial(const real *f, real checkValue, size_t noBufferCell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);
    static bool adaptYDirection(const real *f, real checkValue, size_t noBufferCell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);
    static bool adaptYDirection_serial(const real *f, real checkValue, size_t noBufferCell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);

    void extractData(const std::string &filename, real height, real time);

    void extractData(const std::string &filename);

private:
    bool isUpdateNecessary();

    void applyChanges();

    long m_shift_x1, m_shift_x2, m_shift_y1, m_shift_y2, m_shift_z1, m_shift_z2;
    IAdaptionFunction *func;
    bool m_dynamic, m_dynamic_end;
    bool m_reduction;
    bool m_hasDataExtraction;
    bool m_hasDataExtractionBefore = false;
    bool m_hasDataExtractionAfter = false;
    bool m_hasDataExtractionEndresult = false;
    bool m_hasTimeMeasuring = false;
    bool m_hasWriteField = false;
    bool m_hasWriteRuntime = false;

    real m_height_before = 0;
    real m_height_after = 0;

    size_t m_minimal;

    ISolver *m_solver;

    std::string m_filename;

};


#endif /* ARTSS_ADAPTION_ADAPTION_H_ */
