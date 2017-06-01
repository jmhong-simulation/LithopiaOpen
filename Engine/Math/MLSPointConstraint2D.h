// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"

class MLSPointConstraint2D
{
public:
	int type_;//type of given constraint (0:value,1:value and normal,2:normal component)
	T value, weight_, divergence_;
	TV2 position_, normal_;
	Array1D<TV2> vector_gradient_;
	T epsilon;

public:
	MLSPointConstraint2D(const TV2& position_input, const T& value_input)
	{
		position_ = position_input;
		value = value_input;
		type_ = 0;
		weight_ = (T)1;
		epsilon = -1;
	}

	MLSPointConstraint2D(const TV2& position_input, const T& value_input, const int& type_input)
	{
		position_ = position_input;
		value = value_input;
		type_ = type_input;
		weight_ = (T)1;
		epsilon = -1;
	}

	MLSPointConstraint2D(const TV2& position_input, const T& value_input, const TV2& normal_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = 1;
		weight_ = (T)1;
		epsilon = -1;
	}

	MLSPointConstraint2D(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = (T)1;
		epsilon = -1;
	}

	MLSPointConstraint2D(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input, const T& weight_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = weight_input;
		epsilon = -1;
	}

	MLSPointConstraint2D(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input, const T& weight_input, const T& divergence_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = weight_input;
		divergence_ = divergence_input;
		epsilon = -1;
	}

	MLSPointConstraint2D(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input, const T& weight_input, const T& divergence_input, const Array1D<TV2>& vector_gradient_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = weight_input;
		divergence_ = divergence_input;
		vector_gradient_ = vector_gradient_input;
		epsilon = -1;
	}

	void Set(const TV2& position_input, const T& value_input)
	{
		position_ = position_input;
		value = value_input;
		type_ = 0;
		weight_ = (T)1;
	}

	void Set(const TV2& position_input, const T& value_input, const int& type_input)
	{
		position_ = position_input;
		value = value_input;
		type_ = type_input;
		weight_ = (T)1;
	}

	void Set(const TV2& position_input, const T& value_input, const TV2& normal_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = 1;
		weight_ = (T)1;
	}

	void Set(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = (T)1;
	}

	void Set(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input, const T& weight_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = weight_input;
	}

	void Set(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input, const T& weight_input, const T& divergence_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = weight_input;
		divergence_ = divergence_input;
	}

	void Set(const TV2& position_input, const T& value_input, const TV2& normal_input, const int& type_input, const T& weight_input, const T& divergence_input, const Array1D<TV2>& vector_gradient_input)
	{
		position_ = position_input;
		value = value_input;
		normal_ = normal_input;
		type_ = type_input;
		weight_ = weight_input;
		divergence_ = divergence_input;
		vector_gradient_ = vector_gradient_input;
	}

};//End of class MLS_CONSTRAINT