// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"

class Sphere2D
{
public:
	TV2 center_;
	T	radius_;

public:
	Sphere2D(const TV2& center_, const T& radius_)
		: center_(center_), radius_(radius_)
	{};

	~Sphere2D(void)
	{};

public:
	T getSignedDistance(const TV2& position) const
	{
		return (position - center_).getMagnitude() - radius_;
	}

	TV2 getNormal(const TV2& position) const
	{
		return position - center_;
	}

	TV2 getUnitNormal(const TV2& position) const
	{
		return (position - center_).getNormalized();
	}
};

typedef Sphere2D Circle;