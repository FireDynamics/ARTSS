/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"
#include "../Domain.h"

Field::Field(FieldType type) :
        m_level(0), m_size(Domain::getInstance()->get_size()), m_type(type){
    data = new real[m_size];
#pragma acc enter data create(this)
#pragma acc update device(this)
    // note that the pointer is to host memory, so we overwrite with a
    // pointer to memory allocated on the device.
#pragma acc enter data create(data[:m_size])
}

Field::Field(size_t size) :
        Field::Field(UNKNOWN_FIELD, 0, 0, size) {
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
#pragma acc enter data create(this)
#pragma acc update device(this)
    // note that the pointer is to host memory, so we overwrite with a
    // pointer to memory allocated on the device.
#pragma acc enter data create(data[:m_size])
}

Field::Field(Field const &original):
        data(new real[original.m_size]), m_level(original.m_level),
        m_size(original.m_size), m_type(original.m_type) {
    this->copy_data(original);
#pragma acc enter data copyin(this[:1]) create(data[:m_size])
}

Field::~Field() {
#pragma acc exit data delete(data[:m_size])
#pragma acc exit data delete(this)
    delete[] data;
}
