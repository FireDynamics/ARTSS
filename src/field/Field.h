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
    Field(FieldType type);
    Field(FieldType type, real val);
    Field(FieldType type, real val, size_t level);
    Field(FieldType type, real val, size_t level, size_t size);

    ~Field();

    // getter
    FieldType get_type() { return this->m_type; }
    size_t get_level() { return this->m_level; }

    void set_value(real val) { std::fill(data, data + m_size, val); }
    void copy_data(const Field &other) {
        // TODO parallelise copy data function for GPU
        std::copy(other.data, other.data+other.m_size, data);
    }
    static void swap(Field *a, Field *b) { std::swap(a->data, b->data); }

    real *data;

 private:
    size_t m_level;
    size_t m_size;
    FieldType m_type;
};

#endif /* ARTSS_FIELD_FIELD_H_ */

