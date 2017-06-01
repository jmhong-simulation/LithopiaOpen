#include "StaticTriangle.h"

TV StaticTriangle::getClosestPosition(const TV& location) const
{
	TV nor = crossProduct(v1_ - v0_, v2_ - v0_); nor.normalize(); // normal unit of plane
	TV p = location - nor*dotProduct(nor, location - v0_);

	TV weights = getBarycentricCoordinates(p, v0_, v1_, v2_);// can be reused for velocity calculation

	// project closest point to the triangle if it's not already inside it

	if (weights.y_ < (T)0 && weights.z_ < (T)0)
	{
		return v0_; //return x1 
	}
	else if (weights.x_ < (T)0 && weights.z_ < (T)0)
	{
		return v1_; //return x2 
	}
	else if (weights.x_ < (T)0 && weights.y_ < (T)0)
	{
		return v2_; //return x3 
	}

	if (weights.x_ < (T)0) // Closest point is on edge x2--x3
	{
		return getClosestPointFromLine(p, v1_, v2_);
	}
	else if (weights.y_ < (T)0) // Closest point is on edge x1--x3
	{
		return getClosestPointFromLine(p, v0_, v2_);
	}
	else if (weights.z_ < (T)0) // Closest point is on edge x1--x2
	{
		return getClosestPointFromLine(p, v0_, v1_);
	}

	return p;
}

T  StaticTriangle::getDistance(const TV& location) const
{
	return (getClosestPosition(location) - location).getMagnitude();
}

T  StaticTriangle::getAspectRatio() const
{
	const T l0 = (v0_ - v1_).getMagnitude();
	const T l1 = (v0_ - v2_).getMagnitude();
	const T l2 = (v1_ - v2_).getMagnitude();

	const T s = (l0 + l1 + l2) * (T)0.5;
	const T aspect_ratio = l0 * l1 * l2 / ((T)8 * (s - l0) * (s - l1) * (s - l2));

	return aspect_ratio;
}

TV StaticTriangle::getBarycentricCoordinates(const TV& location, const TV& v0_, const TV& v1_, const TV& v2_) const
{
	TV u = v1_ - v0_, v = v2_ - v0_, w = location - v0_;

	T u_dot_u = dotProduct(u, u), v_dot_v = dotProduct(v, v), u_dot_v = dotProduct(u, v),
	  u_dot_w = dotProduct(u, w), v_dot_w = dotProduct(v, w);

	T denominator = u_dot_u*v_dot_v - POW2(u_dot_v), one_over_denominator;

	if (abs(denominator) > (T)1e-16)
	{
		one_over_denominator = 1 / denominator;
	}
	else
	{
		one_over_denominator = (T)1e16;
	}

	T a = (v_dot_v*u_dot_w - u_dot_v*v_dot_w)*one_over_denominator, b = (u_dot_u*v_dot_w - u_dot_v*u_dot_w)*one_over_denominator;

	return TV(1 - a - b, a, b);
}

void StaticTriangle::getXMinMax(T& x_min, T& x_max) const
{
	x_min = x_max = v0_.x_;

	x_min = MIN2(x_min, v1_.x_);
	x_max = MAX2(x_max, v1_.x_);

	x_min = MIN2(x_min, v2_.x_);
	x_max = MAX2(x_max, v2_.x_);
}

void StaticTriangle::getYMinMax(T& y_min, T& y_max) const
{
	y_min = y_max = v0_.y_;

	y_min = MIN2(y_min, v1_.y_);
	y_max = MAX2(y_max, v1_.y_);

	y_min = MIN2(y_min, v2_.y_);
	y_max = MAX2(y_max, v2_.y_);
}

void StaticTriangle::getZMinMax(T& z_min, T& z_max) const
{
	z_min = z_max = v0_.z_;

	z_min = MIN2(z_min, v1_.z_);
	z_max = MAX2(z_max, v1_.z_);

	z_min = MIN2(z_min, v2_.z_);
	z_max = MAX2(z_max, v2_.z_);
}

T StaticTriangle::getArea() const
{
	return ABS(crossProduct(v2_ - v0_, v1_ - v0_).getMagnitude()) * (T)0.5;

//	return (v0_.x_ * (v1_.y_ - v2_.z_) + v1_.x_ * (v2_.y_ - v0_.y_) + v2_.x_ * (v0_.y_ - v1_.y_)) * (T)0.5;		// 2D
}

TV StaticTriangle::getClosestPointFromLine(const TV& location, const TV& x1, const TV& x2) const
{
	T p = dotProduct(location - x1, x2 - x1) / dotProduct(x2 - x1, x2 - x1);

	if (p < (T)0)
	{
		return x1;
	}
	else if (p >(T)1)
	{
		return x2;
	}

	return x1 + (x2 - x1)*p;
}