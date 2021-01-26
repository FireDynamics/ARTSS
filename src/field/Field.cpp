/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"

Field::Field(FieldType type, real val) :
    Field::Field(type, val, 0, Domain::getInstance()->get_size()) {
}

Field::Field(FieldType type, real val, size_t level) :
    Field::Field(type, val, level, Domain::getInstance()->get_size(level)) {
}

Field::Field(FieldType type, real val, size_t level, size_t size):
    m_level(level), m_type(type), m_size(size) {
    data = new real[m_size];
    set_value(val);
}

void Field::set_value(real val) {
    std::fill(data, data + m_size, val);
}

void Field::copy_data(const Field &other) {
    std::copy(other.data, other.data+other.m_size, data);
}

Field::~Field() {
    delete[] data;
}
