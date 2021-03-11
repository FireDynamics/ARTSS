/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"
#include "../Domain.h"

// Field::Field(FieldType type, real val) :
//     Field::Field(type, val, 0, Domain::getInstance()->get_size()) {
// }
// 
// Field::Field(FieldType type, real val, size_t level) :
//     Field::Field(type, val, level, Domain::getInstance()->get_size(level)) {
// }

Field::Field(FieldType type, real val, size_t level, size_t size):
    m_level(level), m_size(size), m_type(type) {
    data = new real[m_size];
    set_value(val);
    #pragma acc enter data copyin(this[:1]) create(m_data[:m_size])
}

Field::Field(Field const &orig):
    m_level(orig.get_level()), m_size(orig.get_size()), m_type(orig.get_type()),
    data(new real[orig.get_size()]) {
    this->copy_data(orig);
    #pragma acc enter data copyin(this[:1]) create(m_data[:m_size])
}

Field::~Field() {
    #pragma acc exit data delete(m_data[:m_size])
    delete[] data;
}
