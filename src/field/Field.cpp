/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Field.h"
#include "../DomainData.h"

inline static const std::vector<std::string> field_type_names = {"rho", "u", "v", "w", "p", "T", "nu"};

Field::Field(FieldType type) :
        m_level(0), m_size(DomainData::getInstance()->get_size()), m_type(type) {
    data = new real[m_size];
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger("Field_GPU_" + std::to_string(counter++));
    m_gpu_logger->info("{} level {} create with field pointer: {} data pointer: {}",
                        get_field_type_name(m_type), m_level,
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
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger("Field_GPU_" + std::to_string(counter++));
    m_gpu_logger->info("{} level {} create with field pointer: {} data pointer: {}",
                        get_field_type_name(m_type), m_level,
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
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger("Field_GPU_" + std::to_string(counter++));
    m_gpu_logger->info("{} level {} create with field pointer: {} data pointer: {}",
                        get_field_type_name(m_type), m_level,
                        static_cast<void *>(this), static_cast<void *>(data));
#endif
#pragma acc enter data copyin(this[:1]) create(data[:m_size])
    this->copy_data(original);
}

Field::~Field() {
#pragma acc exit data delete(data[:m_size])
#pragma acc exit data delete(this)
    delete[] data;
}

//====================================== Matches ===================================================
// *************************************************************************************************
/// \brief  matches string to field_type_names
/// \param  string           string to be matched
// *************************************************************************************************
FieldType Field::match_field(const std::string &string) {
    for (size_t fn = 0; fn < field_type_names.size(); fn++) {
        if (field_type_names[fn] == string) return (FieldType) fn;
    }
    return UNKNOWN_FIELD;
}

std::string Field::get_field_type_name(FieldType f) {
    if (f != FieldType::UNKNOWN_FIELD) {
        return field_type_names[f];
    } else {
        return "UNKNOWN FIELD";
    }
}

