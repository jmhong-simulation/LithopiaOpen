// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../GENERIC_DEFINITIONS.h"
#include "../CONVENTIONAL_MACROS.h"

template<class TT>
class BOX_2D
{
public:
	union
	{
		struct{TT x_min_, y_min_,x_max_,y_max_;};
		struct{TT i_start_,j_start_,i_end_,j_end_;};
	};

public:
	BOX_2D(void)
	{}

	BOX_2D(const BOX_2D<TT>& box_input)
		: x_min_(box_input.x_min_), y_min_(box_input.y_min_), x_max_(box_input.x_max_), y_max_(box_input.y_max_)
	{}

	BOX_2D(const Vector2D<TT> min_input, const Vector2D<TT> max_input)
		: x_min_(min_input.x_), y_min_(min_input.y_), x_max_(max_input.x_), y_max_(max_input.y_)
	{}

	BOX_2D(const TT& x_min_input, const TT& y_min_input, const TT& x_max_input, const TT& y_max_input)
		: x_min_(x_min_input), y_min_(y_min_input), x_max_(x_max_input), y_max_(y_max_input)
	{}

	~BOX_2D(void)
	{}

	bool ClampInside(TT& x, TT& y);		// clamp (x,y) to be included by this box
	bool Inside(const TT& x, const TT& y);
	bool Inside(const Vector2D<TT>& pos) const;
	void Expand(const TT& width);
	void Include(const Vector2D<TT>& position);
	Vector2D<TT> getMin() const;
	Vector2D<TT> getMax() const;

	BOX_2D<TT> getExtended(const TT& width) const;
};

template<class TT> std::ostream& operator<<(std::ostream& output, const BOX_2D<TT>& box)
{
	return output << "(" << box.x_min_ << "," << box.y_min_ << ")x(" << box.x_max_ << "," << box.y_max_ << ")";
}
