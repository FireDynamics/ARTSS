/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"
#include "../Domain.h"

Field::Field(FieldType type): m_type(type), m_level(0), m_size(Domain::getInstance()->get_size()) {
    data = new real[m_size];
}

Field::Field(FieldType type, real val) :
        Field::Field(type, val, 0, Domain::getInstance()->get_size()) {
}

Field::Field(FieldType type, real val, size_t level) :
        Field::Field(type, val, level, Domain::getInstance()->get_size(level)) {
}

Field::Field(FieldType type, real val, size_t level, size_t size):
        m_level(level), m_size(size), m_type(type) {
    data = new real[m_size];
    set_value(val);
}

Field::~Field() {
    delete[] data;
}

Field::Field(Field const &original):
        data(new real[original.m_size]),
        m_level(original.m_level), m_size(original.m_size), m_type(original.m_type) {
    this->copy_data(original);
}
