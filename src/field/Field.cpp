/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"
#include "../domain/DomainData.h"

int Field::counter = 0;

Field::Field(FieldType type) :
        m_level(0), m_size(DomainData::getInstance()->get_size()), m_type(type) {
    data = new real[m_size];
#ifndef BENCHMARKING
    m_logger = Utility::create_logger("Field_GPU_" + std::to_string(counter++));
    m_logger->debug("{} level {} create with field pointer: {} data pointer: {}",
                        Mapping::get_field_type_name(m_type), m_level,
                        static_cast<void *>(this), static_cast<void *>(data));
#endif
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
        Field::Field(type, val, 0, DomainData::getInstance()->get_size()) {
}

Field::Field(FieldType type, real val, size_t level) :
        Field::Field(type, val, level, DomainData::getInstance()->get_size(level)) {
}

Field::Field(FieldType type, real val, size_t level, size_t size):
        m_level(level), m_size(size), m_type(type) {
    data = new real[m_size];
#ifndef BENCHMARKING
    m_logger = Utility::create_logger("Field_GPU_" + std::to_string(counter++));
    m_logger->debug("{} level {} create with field pointer: {} data pointer: {}",
                        Mapping::get_field_type_name(m_type), m_level,
                        static_cast<void *>(this), static_cast<void *>(data));
#endif
#pragma acc enter data create(this)
#pragma acc update device(this)
    // note that the pointer is to host memory, so we overwrite with a
    // pointer to memory allocated on the device.
#pragma acc enter data create(data[:m_size])
    set_value(val);
}

Field::Field(const Field &original):
        data(new real[original.m_size]), m_level(original.m_level),
        m_size(original.m_size), m_type(original.m_type) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger("Field_GPU_" + std::to_string(counter++));
    m_logger->debug("{} level {} create with field pointer: {} data pointer: {}",
                        Mapping::get_field_type_name(m_type), m_level,
                        static_cast<void *>(this), static_cast<void *>(data));
#endif
#pragma acc enter data copyin(this)
#pragma acc enter data create(data[:m_size])
    this->copy_data(original);
}

Field::~Field() {
#pragma acc exit data delete(data[:m_size])
#pragma acc exit data delete(this)
    delete[] data;
}

