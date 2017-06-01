// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "GridUniform1D.h"

void GridUniform1D::initialize(const int& i_start_input, const int& i_res_input, const T& x_min_input, const T& x_max_input)
{
	i_start_ = i_start_input;

	i_res_ = i_res_input;

	x_min_ = x_min_input;

	x_max_ = x_max_input;

	i_end_ = i_start_ + i_res_ - 1;

	dx_ = (x_max_ - x_min_) / (T)i_res_;

	one_over_dx_ = (T)1 / dx_;

	one_over_2dx_ = (T)0.5 / dx_;
}

T GridUniform1D::getCellCenter(const int& i) const
{
	return x_min_ + ((T)0.5 + (T)(i - i_start_))*dx_;
}

T GridUniform1D::getNodePosition(const int& i) const
{
	return x_min_ + (T)(i - i_start_)*dx_;
}

int GridUniform1D::getNumAllNodes() const
{
	return i_res_ + 1;
}

int GridUniform1D::getNumAllCells() const
{
	return i_res_;
}