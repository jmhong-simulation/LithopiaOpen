// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "BOX_2D.h"

template<class TT>
bool BOX_2D<TT>::ClampInside(TT& x, TT& y)		// clamp (x,y) to be included by this box
{
	bool clamped = false;

	if(x < x_min_)
	{
		x = x_min_;
		clamped = true;
	}
	else if(x > x_max_)
	{
		x = x_max_;
		clamped = true;
	}

	if(y < y_min_)
	{
		y = y_min_;
		clamped = true;
	}
	else if(y > y_max_)
	{
		y = y_max_;
		clamped = true;
	}

	return clamped;
}

template<class TT>
bool BOX_2D<TT>::Inside(const TT& x, const TT& y)
{
	if(x < x_min_) return false;
	else if(x > x_max_) return false;

	if(y < y_min_) return false;
	else if(y > y_max_) return false;

	return true;
}

template<class TT>
bool BOX_2D<TT>::Inside(const Vector2D<TT>& pos) const
{
	if(pos.x_ < x_min_) return false;
	else if(pos.x_ > x_max_) return false;

	if(pos.y_ < y_min_) return false;
	else if(pos.y_ > y_max_) return false;

	return true;
}

template<class TT>
void BOX_2D<TT>::Expand(const TT& width)
{
	x_min_ -= width;
	y_min_ -= width;
	x_max_ += width;
	y_max_ += width;
}

template<class TT>
BOX_2D<TT> BOX_2D<TT>::getExtended(const TT& width) const
{
	return BOX_2D<TT>(x_min_ - width, y_min_ - width, x_max_ + width, y_max_ + width);
}

template<class TT>
Vector2D<TT> BOX_2D<TT>::getMin() const
{
	return Vector2D<TT>(x_min_, y_min_);
}

template<class TT>
Vector2D<TT> BOX_2D<TT>::getMax() const
{
	return Vector2D<TT>(x_max_, y_max_);
}

template<class TT>
void BOX_2D<TT>::Include(const Vector2D<TT>& position)
{
	x_min_ = MIN2(x_min_, (TT)position.x_);
	y_min_ = MIN2(y_min_, (TT)position.y_);

	x_max_ = MAX2(x_max_, (TT)position.x_);
	y_max_ = MAX2(y_max_, (TT)position.y_);
}

template class BOX_2D<int>;
template class BOX_2D<float>;
template class BOX_2D<double>;
