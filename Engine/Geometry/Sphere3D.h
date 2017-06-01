// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"

class Sphere3D
{
public:
	TV center_;
	T  radius_;

	Sphere3D(const TV& _center, const T& _radius)
		: center_(_center), radius_(_radius)
	{}
};