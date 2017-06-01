// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "LineSegment2D.h"

LineSegment2D::LineSegment2D()
{}

LineSegment2D::LineSegment2D(const TV2& _p0, const TV2& _p1)
	: p0_(_p0), p1_(_p1)
{}

void LineSegment2D::Initialize(const TV2& _p0, const TV2& _p1)
{
	p0_ = _p0;
	p1_ = _p1;
}

T LineSegment2D::getDistance(const TV2& v) const
{
	const T mag = (T)(p1_ - p0_).getMagnitude();
	const T t = dotProduct(v - p0_, (p1_ - p0_).getSafeNormalized());

	if (t <= 0.0f) return (v - p0_).getMagnitude();

	if (t >= mag) return (v - p1_).getMagnitude();

	return (T)(v - ((p1_ - p0_).getSafeNormalized()*t + p0_)).getMagnitude();
}

void LineSegment2D::getClosestPoint(const TV2& v, TV2& c, T& t) const
{
	const T mag = (T)(p1_ - p0_).getMagnitude();

	t = dotProduct(v - p0_, (p1_ - p0_)/mag) / mag;		//TODO: normal safety (check p0 == p1 when initialized)

	t = CLAMP(t, (T)0, (T)1);

	c = p0_ * ((T)1 - t) + p1_ * t;
}

BOX_2D<T> LineSegment2D::getAABB() const
{
	const T min_x = MIN2(p0_.x_, p1_.x_);
	const T max_x = MAX2(p0_.x_, p1_.x_);

	const T min_y = MIN2(p0_.y_, p1_.y_);
	const T max_y = MAX2(p0_.y_, p1_.y_);

	return BOX_2D<T>(min_x, min_y, max_x, max_y);
}