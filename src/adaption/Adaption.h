/// \file       Adaption.h
/// \brief      Controller class for adaption
/// \date       Nov 29, 2018
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

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
    explicit Adaption(ISolver *solver);

    bool inline is_data_extraction_enabled() { return m_has_data_extraction; };
    bool inline is_data_extraction_before_enabled() { return m_has_data_extraction_before; }
    bool inline is_data_extraction_after_enabled() { return m_has_data_extraction_after; }
    bool inline is_data_extraction_endresult_enabled() { return m_has_data_extraction_endresult; }
    bool inline is_time_measuring_enabled() { return m_has_time_measuring; }
    bool inline is_write_runtime_enabled() { return m_has_write_runtime; }

    real inline get_before_height() { return m_height_before; }
    std::string inline get_before_name() { return m_filename + "_before.csv"; }
    real inline get_after_height() { return m_height_after; }
    std::string inline get_after_name() { return m_filename + "_after.csv"; }
    std::string inline get_endresult_name() { return m_filename + "_endresult.csv"; }
    std::string inline get_time_measuring_name() { return m_filename + "_time.csv"; }
    std::string inline get_write_runtime_name() { return m_filename + "_runtime.csv"; }

    void run(real t_cur);

    static void expand_x_direction(long shift, bool start, size_t *arr_idx_expansion, size_t len_e);
    static void expand_y_direction(long shift, bool start, size_t *arr_idx_expansion, size_t len_e);
    static void reduce_x_direction(long shift, bool start, size_t *arr_idx_reduction, size_t len_r);
    static void reduce_y_Direction(long shift, bool start, size_t *arr_idx_reduction, size_t len_r);
    static bool adapt_x_direction(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);
    static bool adapt_x_direction_serial(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);
    static bool adapt_y_direction(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);
    static bool adapt_y_direction_serial(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce);

    void extractData(const std::string &filename, real height, real time);
    void extractData(const std::string &filename);

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    bool isUpdateNecessary();
    void applyChanges();

    long m_shift_x1, m_shift_x2, m_shift_y1, m_shift_y2, m_shift_z1, m_shift_z2;
    bool m_dynamic, m_dynamic_end;
    bool m_has_data_extraction;
    bool m_has_data_extraction_before = false;
    bool m_has_data_extraction_after = false;
    bool m_has_data_extraction_endresult = false;
    bool m_has_time_measuring = false;
    bool m_has_write_runtime = false;

    real m_height_before = 0;
    real m_height_after = 0;

    size_t m_minimal;

    std::string m_filename;

    ISolver *m_solver;
    IAdaptionFunction *func;
};

#endif /* ARTSS_ADAPTION_ADAPTION_H_ */
