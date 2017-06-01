#include "Triangle2D.h"

BOX_2D<T> Triangle2D::getAABB() const
{
	T min_x = MIN3(v0_.x_, v1_.x_, v2_.x_);
	T max_x = MAX3(v0_.x_, v1_.x_, v2_.x_);

	T min_y = MIN3(v0_.y_, v1_.y_, v2_.y_);
	T max_y = MAX3(v0_.y_, v1_.y_, v2_.y_);
	
	return BOX_2D<T>(min_x, min_y, max_x, max_y);
}

TV3 Triangle2D::getBarycentricCoordinates(const TV2& location) const
{
	TV2 u = v1_ - v0_, v = v2_ - v0_, w = location - v0_;

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

	return TV3(1 - a - b, a, b);
}

bool Triangle2D::isInside(const TV2& p) const
{
	const TV3 bc = getBarycentricCoordinates(p);

	if (bc.x_ < (T)0) return false;
	else if (bc.x_ > (T)1) return false;
	else if (bc.y_ < (T)0) return false;
	else if (bc.y_ > (T)1) return false;
	else if (bc.z_ < (T)0) return false;
	else if (bc.z_ > (T)1) return false;
	else return true;
}

int Triangle2D::getInsideFlag(const TV2& p) const
{
	if (isInside(p) == false) return 0;
	else if (getSignedArea() > (T)0) return 1;
	else return -1;
}