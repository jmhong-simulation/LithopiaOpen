// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "Geometry/LineSegment2D.h"

class MLSLineConstraint2D
{
public:
	LineSegment2D line_;
	T v0_, v1_;
	T weight_;

	MLSLineConstraint2D()
		: weight_((T)1)
	{}

	MLSLineConstraint2D(const TV2& _p0, const TV2& _p1, const T& _v0, const T& _v1, const T& _w)
		: line_(_p0, _p1), v0_(_v0), v1_(_v1), weight_(_w)
	{}
};