/// \file         Field.h
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FIELD_FIELD_H_
#define ARTSS_FIELD_FIELD_H_

#include <algorithm>
#include <utility>
#include "../utility/GlobalMacrosTypes.h"

#ifndef ENUM_TYPES
#define ENUM_TYPES
const size_t numberOfFieldTypes = 7;
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

    // read_ptr get_data_ro() const { return data; }
    size_t get_level() const { return m_level; }
    size_t get_size() const { return m_size; }

    // setter
    inline real& operator[](size_t i) const { return data[i]; }

    // acc functions
    void update_host() {
        #pragma acc update host(data[:m_size])
    }
    void update_dev() {
        #pragma acc update device(data[:m_size])
    }
    void copyin() {
        #pragma acc enter data copyin(data[:m_size])
    }

    void set_value(real val) const { std::fill(data, data + m_size, val); }
    void copy_data(const Field &other) const {
        // TODO parallelise copy data function for GPU
        std::copy(other.data, other.data + other.m_size, data);
    }
    static void swap(Field &a, Field &b) { std::swap(a.data, b.data); }

    Field &operator+=(const real x) {
#pragma acc parallel loop independent present(this->data[:m_size]) async
        for (size_t i=0; i < m_size; ++i)
            this->data[i] += x;

#pragma acc wait
        return *this;
    }

    Field &operator+=(const Field &rhs) {
        auto rhs_data = rhs.data;
#pragma acc parallel loop independent present(this->data[:m_size], rhs_data[:m_size]) async
        for (size_t i=0; i < m_size; ++i)
            this->data[i] += rhs_data[i];

#pragma acc wait
        return *this;
    }

    Field &operator*=(const real x) {
#pragma acc parallel loop independent present(this->data[:m_size]) async
        for (size_t i=0; i < m_size; ++i)
            this->data[i] *= x;

#pragma acc wait
        return *this;
    }

    Field &operator*=(const Field &rhs) {
        auto rhs_data = rhs.data;
#pragma acc parallel loop independent present(this->data[:m_size], rhs_data[:m_size]) async
        for (size_t i=0; i < m_size; ++i)
            this->data[i] *= rhs_data[i];

#pragma acc wait
        return *this;
    }

    real *data;

 private:
    size_t const m_level;
    size_t const m_size;
    FieldType const m_type;
};

#endif /* ARTSS_FIELD_FIELD_H_ */
