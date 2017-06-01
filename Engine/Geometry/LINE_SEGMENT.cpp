// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "LINE_SEGMENT.h"

LINE_SEGMENT::LINE_SEGMENT()
{}

LINE_SEGMENT::LINE_SEGMENT(const TV& _p0, const TV& _p1)
	: p0_(_p0), p1_(_p1)
{}

void LINE_SEGMENT::Initialize(const TV& _p0, const TV& _p1)
{
	p0_ = _p0;
	p1_ = _p1;
}

T LINE_SEGMENT::getDistance(const TV& v)
{
	const T mag = (T)(p1_ - p0_).getMagnitudeDouble();
	const T t = dotProduct(v - p0_, (p1_ - p0_).normalizedDouble());

	if (t <= 0.0f) return (v - p0_).getMagnitude();

	if (t >= mag) return (v - p1_).getMagnitude();

	return (T)(v - ((p1_ - p0_).normalizedDouble()*t + p0_)).getMagnitudeDouble();
}