// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "BOX_2D.h"

class LineSegment2D
{
public:
	TV2 p0_, p1_;

	LineSegment2D();

	LineSegment2D(const TV2& _p0, const TV2& _p1);

	void Initialize(const TV2& _p0, const TV2& _p1);

	T getDistance(const TV2& v) const;

	void getClosestPoint(const TV2& v, TV2& c, T& t) const;

	BOX_2D<T> getAABB() const;
};
