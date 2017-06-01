// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../DataStructure/GridUniform2D.h"
#include "../DataStructure/Array2D.h"

class SMOOTHING_UNIFORM_2D
{
public:
	void SmoothIClampJRepeat(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat);
	void SmoothIRepeatJClamp(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat);
	void SmoothIRepeatJClampCell(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat);
	void SmoothDirichletIJ(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat);
	void SmoothSignedDistance(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat, const T th);
};