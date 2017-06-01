#include "Quadrilateral2D.h"

BOX_2D<T> Quadrilateral2D::getAABB() const
{
	T min_x = MIN4(v0_.x_, v1_.x_, v2_.x_, v3_.x_);
	T max_x = MAX4(v0_.x_, v1_.x_, v2_.x_, v3_.x_);

	T min_y = MIN4(v0_.y_, v1_.y_, v2_.y_, v3_.y_);
	T max_y = MAX4(v0_.y_, v1_.y_, v2_.y_, v3_.y_);

	return BOX_2D<T>(min_x, min_y, max_x, max_y);
}