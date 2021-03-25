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
const size_t numberOfFieldTypes = 6;
enum FieldType : int {
    UNKNOWN_FIELD = -1, RHO = 0, U = 1, V = 2, W = 3, P = 4, T = 5
};
#endif

class Field {
 public:
    Field(FieldType type, real val, size_t level, size_t size);
    explicit Field(Field const &orig);

    ~Field();

    // getter
    FieldType getType() const { return m_type; }
    return_ptr getData() const { return data; }

    // read_ptr get_data_ro() const { return data; }
    size_t getLevel() const { return m_level; }
    size_t getSize() const { return m_size; }

    // setter
    inline real& operator[](size_t i) const { return data[i]; }

    // acc functions
    void updateHost() {
        #pragma acc update host(data[:m_size])
    }
    void updateDev() {
        #pragma acc update device(data[:m_size])
    }
    void copyin() {
        #pragma acc enter data copyin(data[:m_size])
    }

    // basic interface to algorithm
    void setValue(real val) { std::fill(data, data + m_size, val); }
    void copyData(Field const &other) {
        std::copy(other.data, other.data+other.m_size, data);
    }
    static void swap(Field &a, Field &b) { std::swap(a.data, b.data); }

    real *data;

 private:
    size_t const m_level;
    size_t const m_size;

    FieldType const m_type;
};

#endif /* ARTSS_FIELD_FIELD_H_ */

