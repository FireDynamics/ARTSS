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
const size_t numberOfFieldTypes = 7;
enum FieldType : int {
    UNKNOWN_FIELD = -1, RHO = 0, U = 1, V = 2, W = 3, P = 4, T = 5, NU = 6
};
#endif

template <typename T>
struct counter
{
    counter()
    {
        objects_created++;
        objects_alive++;
    }

    counter(const counter&)
    {
        objects_created++;
        objects_alive++;
    }

protected:
    virtual ~counter()
    {
        --objects_alive;
    }
    static int objects_created;
    static int objects_alive;
};
template <typename T> int counter<T>::objects_created( 0 );
template <typename T> int counter<T>::objects_alive( 0 );

class Field : counter<Field> {
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
    void update_host() {
#pragma acc update host(data[:m_size])
    }
    void update_dev() {
#pragma acc update device(data[:m_size])
    }
    void copyin() {
#pragma acc enter data copyin(data[:m_size])
#ifndef BENCHMARKING
  #ifdef _OPENACC
    m_gpu_logger->debug("{}{} copyin with data pointer: {}", get_field_type_name(m_type), m_level, static_cast<void *>(data));
  #endif
#endif
    }

    void set_value(real val) const { std::fill(data, data + m_size, val); }

    void copy_data(const Field &other) const {
        auto other_data = other.data;
#pragma acc parallel loop independent present(this->data[:m_size], other_data[:m_size]) async
        for (size_t i = 0; i < m_size; ++i) {
            this->data[i] = other_data[i];
        }
#pragma acc wait
    }

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

#ifndef BENCHMARKING
#ifdef _OPENACC
    std::shared_ptr<spdlog::logger> m_gpu_logger;
#endif
#endif
};

#endif /* ARTSS_FIELD_FIELD_H_ */
