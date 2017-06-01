// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"

class LINE_SEGMENT
{
public:
	TV p0_, p1_;

	LINE_SEGMENT();

	LINE_SEGMENT(const TV& _p0, const TV& _p1);

	void Initialize(const TV& _p0, const TV& _p1);

	T getDistance(const TV& v);
};
