// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../GENERIC_DEFINITIONS.h"

class RAY
{
public:
	TV origin_;
	TV direction_;
	T t_end_;

public:
	RAY(const TV& _origin, const TV& _direction, const T& _t_end)
		: origin_(_origin), direction_(_direction), t_end_(_t_end)
	{}

	RAY(const TV& _start, const TV& _end)
//		: origin_(_start), t_end_((_end - _start).GetMagnitude()), direction_((_end-_start)/t_end_)
	{
		origin_ = _start;
		t_end_ = (_end - _start).getMagnitude();
		direction_ = (_end - _start) / t_end_;
	}

	const T		GetSphereIntersection(const TV& center, const T& radius) const;
	const TV	GetPosition(const T& t) const;
	const bool	CheckTriangleIntersection(const TV& tri_v0, const TV& tri_v1, const TV& tri_v2, const RAY& segment) const;
};
