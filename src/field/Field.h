/// \file         Field.h
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FIELD_FIELD_H_
#define ARTSS_FIELD_FIELD_H_

#include <algorithm>
#include <utility>
#include <iostream>
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"

#ifndef ENUM_TYPES
#define ENUM_TYPES
const size_t number_of_field_types = 7;
enum FieldType : int {
    UNKNOWN_FIELD = -1, RHO = 0, U = 1, V = 2, W = 3, P = 4, T = 5, NU = 6
};
#endif

class Field {
 public:
    explicit Field(size_t size);
    explicit Field(FieldType type);
    Field(FieldType type, real val);
    Field(FieldType type, real val, size_t level);
    Field(FieldType type, real val, size_t level, size_t size);
    Field(Field const &original);
    ~Field();

    // getter
    FieldType get_type() const { return m_type; }
    return_ptr get_data() const { return data; }
    size_t get_level() const { return m_level; }
    size_t get_size() const { return m_size; }

    // setter
    inline real &operator[](size_t i) const { return data[i]; }

    // acc functions
    /// \brief update data array on CPU
    void update_host() {
#pragma acc update host(data[:m_size])
    }

    /// \brief update data array on GPU
    void update_dev() {
#pragma acc update device(data[:m_size])
    }

    /// \brief copy data array to GPU
    void copyin() {
#pragma acc enter data copyin(data[:m_size])
#ifdef GPU_DEBUG
    m_gpu_logger->debug("{}{} copyin with data pointer: {}", get_field_type_name(m_type), m_level, static_cast<void *>(data));
#endif
    }

    /// \brief fill data array with the specified number
    /// \param val value to be set
    void set_value(real val) const {
#ifdef GPU_DEBUG
        m_gpu_logger->debug("set value to {} for {}", val, static_cast<void *>(data));
#endif
#pragma acc parallel loop independent present(this->data[:m_size]) async
        for (size_t i = 0; i < m_size; i++) {
            this->data[i] = val;
        }
#pragma acc wait
    }

    /// \brief copy data of Field other. No changes in data pointer
    /// \param Field &other
    void copy_data(const Field &other) const {
        auto other_data = other.data;
#ifdef GPU_DEBUG
        m_gpu_logger->debug("copy data of {} ({}) to {} ({}) for level {}",
                            get_field_type_name(other.get_type()), static_cast<void *>(other_data),
                            get_field_type_name(m_type), static_cast<void *>(data),
                            m_level);
#endif
#pragma acc parallel loop independent present(this->data[:m_size], other_data[:m_size]) async
        for (size_t i = 0; i < m_size; ++i) {
            this->data[i] = other_data[i];
        }
#pragma acc wait
    }

    /// \brief swap pointer of data array of a and b
    /// \warning Use with care! data pointer have to be swapped back in the same
    /// data region otherwise it leads to a `cuCtxSynchronize returned error 700`
    /// \param Field a
    /// \param Field b
    static void swap(Field &a, Field &b) { std::swap(a.data, b.data); }

    Field &operator+=(const real x) {
#pragma acc parallel loop independent present(this->data[:m_size]) async
        for (size_t i = 0; i < m_size; ++i) {
            this->data[i] += x;
        }
#pragma acc wait
        return *this;
    }

    Field &operator+=(const Field &rhs) {
        auto rhs_data = rhs.data;
#pragma acc parallel loop independent present(this->data[:m_size], rhs_data[:m_size]) async
        for (size_t i = 0; i < m_size; ++i) {
            this->data[i] += rhs_data[i];
        }
#pragma acc wait
        return *this;
    }

    Field &operator*=(const real x) {
#pragma acc parallel loop independent present(this->data[:m_size]) async
        for (size_t i = 0; i < m_size; ++i) {
            this->data[i] *= x;
        }
#pragma acc wait
        return *this;
    }

    Field &operator*=(const Field &rhs) {
        auto rhs_data = rhs.data;
#pragma acc parallel loop independent present(this->data[:m_size], rhs_data[:m_size]) async
        for (size_t i = 0; i < m_size; ++i) {
            this->data[i] *= rhs_data[i];
        }
#pragma acc wait
        return *this;
    }

    real *data;
    static std::string get_field_type_name(FieldType f);
    static FieldType match_field(const std::string& string);

 private:
    static int counter;
    size_t const m_level;
    size_t const m_size;
    FieldType const m_type;

#ifdef GPU_DEBUG
    std::shared_ptr<spdlog::logger> m_gpu_logger;
#endif
};

#endif /* ARTSS_FIELD_FIELD_H_ */
