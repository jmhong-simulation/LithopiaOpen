// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

template<class TT>
class Vector4D
{
public:
	TT x_, y_, z_, w_;

	Vector4D()
		: x_(TT()), y_(TT()), z_(TT()), w_(TT())
	{}

	Vector4D(const TT _x, const TT _y, const TT _z, const TT _w)
		: x_(_x), y_(_y), z_(_z), w_(_w)
	{}

	void operator *= (const TT s)
	{
		x_ *= s;
		y_ *= s;
		z_ *= s;
		w_ *= s;
	}

	friend std::ostream& operator<< (std::ostream& stream, const Vector4D<TT>& vec)
	{
		stream << vec.x_ << " " << vec.y_ << " " << vec.z_ << " "<< vec.w_;

		return stream;
	}
};