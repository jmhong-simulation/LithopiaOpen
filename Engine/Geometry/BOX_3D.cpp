// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "BOX_3D.h"
#include "CONVENTIONAL_MACROS.h"

template<class TT>
void BOX_3D<TT>::Initialize(const TT& x_min_input, const TT& y_min_input, const TT& z_min_input, const TT& x_max_input, const TT& y_max_input, const TT& z_max_input)
{
	x_min_ = x_min_input;
	y_min_ = y_min_input;
	z_min_ = z_min_input;
	
	x_max_ = x_max_input;
	y_max_ = y_max_input;
	z_max_ = z_max_input;
}

template<class TT>
void BOX_3D<TT>::Initialize(const Vector3D<TT>& center, const TT& width)
{
	const TT half_width = width / 2;

	Initialize(center.x_ - half_width, center.y_ - half_width, center.z_ - half_width,
			   center.x_ + half_width, center.y_ + half_width, center.z_ + half_width);	
}

template<class TT>
void BOX_3D<TT>::Initialize(const BOX_3D<TT>& box)
{
	x_min_ = box.x_min_;
	y_min_ = box.y_min_;
	z_min_ = box.z_min_;

	x_max_ = box.x_max_;
	y_max_ = box.y_max_;
	z_max_ = box.z_max_;
}

template<class TT>
BOX_3D<TT> BOX_3D<TT>::getZResized(const TT& z_min_input, const TT& z_max_input) const
{
	return BOX_3D<TT>(x_min_, y_min_, z_min_input, x_max_, y_max_, z_max_input);
}

template<class TT>
bool BOX_3D<TT>::ClampInside(TT& x, TT& y, TT& z)		// clamp (x,y,z) to be included by this box
{
	bool clamped = false;

	if (x < x_min_)
	{
		x = x_min_;
		clamped = true;
	}
	else if (x > x_max_)
	{
		x = x_max_;
		clamped = true;
	}

	if (y < y_min_)
	{
		y = y_min_;
		clamped = true;
	}
	else if (y > y_max_)
	{
		y = y_max_;
		clamped = true;
	}

	if (z < z_min_)
	{
		z = z_min_;
		clamped = true;
	}
	else if (z > z_max_)
	{
		z = z_max_;
		clamped = true;
	}

	return clamped;
}

template<class TT>
TT BOX_3D<TT>::getSignedDistance(const Vector3D<TT>& p) const
{
	const Vector3D<TT> box_center = Vector3D<TT>((x_min_ + x_max_)*(TT)0.5, (y_min_ + y_max_)*(TT)0.5, (z_min_ + z_max_)*(TT)0.5);
	const Vector3D<TT> box_half_edge = Vector3D<TT>((x_max_ - x_min_)*(TT)0.5, (y_max_ - y_min_)*(TT)0.5, (z_max_ - z_min_)*(TT)0.5);
	const Vector3D<TT> d = (p - box_center).getCompAbs() - box_half_edge;

	return d.getCompMax((TT)0).getMagnitude() + MIN2(MAX3(d.x_, d.y_, d.z_), 0.0);
}

template<class TT>
Vector3D<TT> BOX_3D<TT>::getNormal(const Vector3D<TT>& p, const TT& epsilon) const
{
	return TTV(getSignedDistance(p + TTV(epsilon, 0, 0)) - getSignedDistance(p + TTV(-epsilon, 0, 0)), getSignedDistance(p + TTV(0, epsilon, 0)) - getSignedDistance(p + TTV(0, -epsilon, 0)),
			   getSignedDistance(p + TTV(0, 0, epsilon)) - getSignedDistance(p + TTV(0, 0, -epsilon))).getSafeNormalized();
}

template<class TT>
bool BOX_3D<TT>::isInside(const TT& x, const TT& y, const TT& z)
{
	if (x < x_min_) return false;
	else if (x > x_max_) return false;

	if (y < y_min_) return false;
	else if (y > y_max_) return false;

	if (z < z_min_) return false;
	else if (z > z_max_) return false;

	return true;
}

template<class TT>
bool BOX_3D<TT>::isInside(const Vector3D<TT>& pos) const
{
	if (pos.x_ < x_min_) return false;
	else if (pos.x_ > x_max_) return false;

	if (pos.y_ < y_min_) return false;
	else if (pos.y_ > y_max_) return false;

	if (pos.z_ < z_min_) return false;
	else if (pos.z_ > z_max_) return false;

	return true;
}

template<class TT>
bool BOX_3D<TT>::HasVolume() const
{
	if (k_end_ < k_start_ || j_end_ < j_start_ || i_end_ < i_start_) return false;

	return true;
}

template<class TT>
void BOX_3D<TT>::Extend(const TT& width)
{
	x_min_ -= width;
	y_min_ -= width;
	z_min_ -= width;

	x_max_ += width;
	y_max_ += width;
	z_max_ += width;
}

template<class TT>
void BOX_3D<TT>::Extend(const Vector3D<TT>& width)
{
	x_min_ -= width.x_;
	y_min_ -= width.y_;
	z_min_ -= width.z_;

	x_max_ += width.x_;
	y_max_ += width.y_;
	z_max_ += width.z_;
}

template<class TT>
BOX_3D<TT> BOX_3D<TT>::GetExtended(const TT& width) const
{
	BOX_3D<TT> box(*this);
	box.Extend(width);
	return box;
}

template<class TT>
BOX_3D<TT> BOX_3D<TT>::GetExtended(const Vector3D<TT>& width) const
{
	BOX_3D<TT> box(*this);
	box.Extend(width);
	return box;
}

template<class TT>
void BOX_3D<TT>::Include(const Vector3D<TT>& position)
{
	x_min_ = MIN2(x_min_, (TT)position.x_);
	y_min_ = MIN2(y_min_, (TT)position.y_);
	z_min_ = MIN2(z_min_, (TT)position.z_);

	x_max_ = MAX2(x_max_, (TT)position.x_);
	y_max_ = MAX2(y_max_, (TT)position.y_);
	z_max_ = MAX2(z_max_, (TT)position.z_);
}

template<class TT>
void BOX_3D<TT>::EnlargeMIN(const TT& x, const TT& y, const TT& z)
{
	x_min_ = MIN2(x_min_, x);
	y_min_ = MIN2(y_min_, y);
	z_min_ = MIN2(z_min_, z);
}

template<class TT>
void BOX_3D<TT>::EnlargeMAX(const TT& x, const TT& y, const TT& z)
{
	x_max_ = MAX2(x_max_, x);
	y_max_ = MAX2(y_max_, y);
	z_max_ = MAX2(z_max_, z);
}


template<class TT>
void BOX_3D<TT>::Enlarge(const BOX_3D<TT>& box_to_include)
{
	x_min_ = MIN2(x_min_, box_to_include.x_min_);
	y_min_ = MIN2(y_min_, box_to_include.y_min_);
	z_min_ = MIN2(z_min_, box_to_include.z_min_);

	x_max_ = MAX2(x_max_, box_to_include.x_max_);
	y_max_ = MAX2(y_max_, box_to_include.y_max_);
	z_max_ = MAX2(z_max_, box_to_include.z_max_);
}

template<class TT>
void BOX_3D<TT>::Translate(const Vector3D<TT>& deviation)
{
	x_min_ += deviation.x_;
	y_min_ += deviation.y_;
	z_min_ += deviation.z_;

	x_max_ += deviation.x_;
	y_max_ += deviation.y_;
	z_max_ += deviation.z_;
}

template<class TT>
void BOX_3D<TT>::Scale(const TT s)
{
	x_min_ = x_min_ * s;
	y_min_ = y_min_ * s;
	z_min_ = z_min_ * s;

	x_max_ = x_max_ * s;
	y_max_ = y_max_ * s;
	z_max_ = z_max_ * s;
}

template<class TT>
void BOX_3D<TT>::Scale(const Vector3D<TT>& s)
{
	x_min_ = x_min_ * s.x_;
	y_min_ = y_min_ * s.y_;
	z_min_ = z_min_ * s.z_;

	x_max_ = x_max_ * s.x_;
	y_max_ = y_max_ * s.y_;
	z_max_ = z_max_ * s.z_;
}

template<class TT>
Vector3D<TT> BOX_3D<TT>::GetEdgeLengths() const
{
	return Vector3D<TT>(x_max_ - x_min_, y_max_ - y_min_, z_max_ - z_min_);
}

template<class TT>
Vector3D<TT> BOX_3D<TT>::GetMin() const
{
	return Vector3D<TT>(x_min_, y_min_, z_min_);
}

template<class TT>
Vector3D<TT> BOX_3D<TT>::GetMax() const
{
	return Vector3D<TT>(x_max_, y_max_, z_max_);
}

template<class TT>
Vector3D<TT> BOX_3D<TT>::GetCenter() const
{
	return Vector3D<TT>((x_min_ + x_max_)*0.5f, (y_min_ + y_max_)*0.5f, (z_min_ + z_max_)*0.5f);
}

template<class TT>
TT BOX_3D<TT>::GetMaxUnityScale() const
{
	float scale = MAX3(x_max_ - x_min_, y_max_ - y_min_, z_max_ - z_min_);
	if (scale != 0.0f) scale = 1.0f / scale;

	return scale;
}

template<class TT> std::ostream& operator<<(std::ostream& output, const BOX_3D<TT>& box)
{
	return output << "(" << box.x_min_ << "," << box.y_min_ << "," << box.z_min_ << ")x(" << box.x_max_ << "," << box.y_max_ << "," << box.z_max_ << ")";
}

template<class TT> BOX_3D<TT> GetIntersection(const BOX_3D<TT>& box_a, const BOX_3D<TT>& box_b)
{
	BOX_3D<TT> box;

	box.i_start_ = MAX(box_a.i_start_, box_b.i_start_);
	box.j_start_ = MAX(box_a.j_start_, box_b.j_start_);
	box.k_start_ = MAX(box_a.k_start_, box_b.k_start_);

	box.i_end_ = MIN(box_a.i_end_, box_b.i_end_);
	box.j_end_ = MIN(box_a.j_end_, box_b.j_end_);
	box.k_end_ = MIN(box_a.k_end_, box_b.k_end_);

	return box;
}

template class BOX_3D<int>;
template class BOX_3D<float>;
template class BOX_3D<double>;
template bool BOX_3D<float>::isInside(const Vector3D<float>& pos) const;



