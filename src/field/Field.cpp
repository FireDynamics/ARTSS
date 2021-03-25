/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"
#include "../Domain.h"

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

Field::Field(Field const &orig):
    data(new real[orig.get_size()]),
    m_level(orig.get_level()), m_size(orig.get_size()), m_type(orig.get_type()) {
    this->copy_data(orig);
#pragma acc enter data copyin(this[:1]) create(m_data[:m_size])
}

Field::~Field() {
#pragma acc exit data delete(data[:m_size])
#pragma acc exit data delete(this)
    delete[] data;
}
