/// \file         Field.cpp
/// \brief        structure including information on cells, lengths, resolution and data for each direction
/// \date         May 25, 2016
/// \author       Severt
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include "Field.h"
#include "Domain.h"

Field::Field(FieldType type, real val): m_type(type){
    size_t m_size = Domain::getInstance()->get_size();
	m_level = 0;
    data = new real[m_size];
    std::fill( data, data + m_size, val);
}

Field::Field(FieldType type, real val, size_t level): m_level(level), m_type(type) {
	size_t m_size = Domain::getInstance()->get_size(level);

	data = new real[m_size];
	std::fill( data, data + m_size, val);
}


Field::~Field(){
	delete[] data;
}


//=============================== Copy Constructor ======================================
Field::Field(const Field & other){
    size_t size = Domain::getInstance()->get_size(other.m_level);
	data = new real[size];
	for (size_t i = 0; i < size; i++){
		data[i] = other.data[i];
	}
	m_type = other.m_type;
	m_level = other.m_level;
}
