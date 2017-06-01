// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../DataStructure/Vector3D.h"
#include <iostream>

template<class TT>
class BOX_3D
{
public:

	typedef Vector3D<TT> TTV;

	union
	{
		struct{ TT x_min_, y_min_, z_min_, x_max_, y_max_, z_max_; };
		struct{ TT i_start_, j_start_, k_start_, i_end_, j_end_, k_end_; };
	};

public:
	BOX_3D(void)
	{}

	BOX_3D(const BOX_3D<TT>& box_input)
		: x_min_(box_input.x_min_), y_min_(box_input.y_min_), z_min_(box_input.z_min_),
		  x_max_(box_input.x_max_), y_max_(box_input.y_max_), z_max_(box_input.z_max_)
	{}

	BOX_3D(const Vector3D<TT> min_input, const Vector3D<TT> max_input)
		: x_min_(min_input.x_), y_min_(min_input.y_), z_min_(min_input.z_), x_max_(max_input.x_), y_max_(max_input.y_), z_max_(max_input.z_)
	{}

	BOX_3D(const TT& x_min_input, const TT& y_min_input, const TT& z_min_input,
		   const TT& x_max_input, const TT& y_max_input, const TT& z_max_input)
		   : x_min_(x_min_input), y_min_(y_min_input), z_min_(z_min_input), x_max_(x_max_input), y_max_(y_max_input), z_max_(z_max_input)
	{}

	BOX_3D(const Vector3D<TT>& center, const TT& width)
	{
		Initialize(center, width);
	}

	~BOX_3D(void)
	{}

	void Initialize(const TT& x_min_input, const TT& y_min_input, const TT& z_min_input, const TT& x_max_input, const TT& y_max_input, const TT& z_max_input);
	void Initialize(const Vector3D<TT>& center, const TT& width);
	void Initialize(const BOX_3D<TT>& box);

	bool ClampInside(TT& x, TT& y, TT& z);		// clamp (x,y,z) to be included by this box
	bool isInside(const TT& x, const TT& y, const TT& z);
	bool isInside(const Vector3D<TT>& pos) const;
	void Extend(const TT& width);
	void Extend(const Vector3D<TT>& width);
	BOX_3D<TT> GetExtended(const TT& width) const;
	BOX_3D<TT> GetExtended(const Vector3D<TT>& width) const;
	BOX_3D<TT> getZResized(const TT& z_min_input, const TT& z_max_input) const;
	void Include(const Vector3D<TT>& position);
	void EnlargeMIN(const TT& x, const TT& y, const TT& z);
	void EnlargeMAX(const TT& x, const TT& y, const TT& z);
	void Enlarge(const BOX_3D<TT>& box_to_include);			// enlarge this box to include input box
	bool HasVolume() const;		// has volume

	void Translate(const Vector3D<TT>& deviation);
	void Scale(const TT s);
	void Scale(const Vector3D<TT>& s);

	Vector3D<TT> GetEdgeLengths() const;
	Vector3D<TT> GetMin() const;
	Vector3D<TT> GetMax() const;

	Vector3D<TT> GetCenter() const;
	TT GetMaxUnityScale() const;

	TT getSignedDistance(const TTV& p) const;
	TTV getNormal(const TTV& p, const TT& epsilon) const;
};

template<class TT> std::ostream& operator<<(std::ostream& output, const BOX_3D<TT>& box);
template<class TT> BOX_3D<TT> GetIntersection(const BOX_3D<TT>& box_a, const BOX_3D<TT>& box_b);