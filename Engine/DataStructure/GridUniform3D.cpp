// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "GridUniform3D.h"
#include "Array1D.h"

GridUniform3D::GridUniform3D()
{}

GridUniform3D::GridUniform3D(const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input,
	const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input)
{
	initialize(i_start_input, j_start_input, k_start_input, i_res_input, j_res_input, k_res_input, x_min_input, y_min_input, z_min_input, x_max_input, y_max_input, z_max_input);
}

GridUniform3D::GridUniform3D(const TV_INT& ijk_start_input, const TV_INT& ijk_res_input,
	const TV& xyz_min_input, const TV& xyz_max_input)
{
	initialize(ijk_start_input.i_, ijk_start_input.j_, ijk_start_input.k_, ijk_res_input.i_, ijk_res_input.j_, ijk_res_input.k_, xyz_min_input.x_, xyz_min_input.y_, xyz_min_input.z_, xyz_max_input.x_, xyz_max_input.y_, xyz_max_input.z_);
}

GridUniform3D::GridUniform3D(const GridUniform3D& grid_input)
{
	initialize(grid_input.i_start_, grid_input.j_start_, grid_input.k_start_, grid_input.i_res_, grid_input.j_res_, grid_input.k_res_,
		grid_input.x_min_, grid_input.y_min_, grid_input.z_min_, grid_input.x_max_, grid_input.y_max_, grid_input.z_max_);
}

GridUniform3D::GridUniform3D(const BOX_3D<T>& bb, const TV_INT& res)
{
	initialize(0, 0, 0, res.x_, res.y_, res.z_, bb.x_min_, bb.y_min_, bb.z_min_, bb.x_max_, bb.y_max_, bb.z_max_);
}

GridUniform3D::~GridUniform3D(void)
{}

void GridUniform3D::initialize(const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input,
	const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input)
{
	i_start_ = i_start_input;
	j_start_ = j_start_input;
	k_start_ = k_start_input;

	i_res_ = i_res_input;
	j_res_ = j_res_input;
	k_res_ = k_res_input;

	x_min_ = x_min_input;
	y_min_ = y_min_input;
	z_min_ = z_min_input;

	x_max_ = x_max_input;
	y_max_ = y_max_input;
	z_max_ = z_max_input;

	i_end_ = i_start_ + i_res_ - 1;
	j_end_ = j_start_ + j_res_ - 1;
	k_end_ = k_start_ + k_res_ - 1;

	dx_ = (x_max_ - x_min_) / (T)i_res_;
	dy_ = (y_max_ - y_min_) / (T)j_res_;
	dz_ = (z_max_ - z_min_) / (T)k_res_;

	one_over_dx_ = (T)1 / dx_;
	one_over_dy_ = (T)1 / dy_;
	one_over_dz_ = (T)1 / dz_;

	one_over_2dx_ = (T)0.5 / dx_;
	one_over_2dy_ = (T)0.5 / dy_;
	one_over_2dz_ = (T)0.5 / dz_;

	ij_res_ = i_res_*j_res_;
	ijk_res_ = ij_res_*k_res_;
}

void GridUniform3D::initialize(const TV_INT& start_input, const TV_INT& res_input, const TV& min_input, const TV& max_input)
{
	initialize(start_input.i_, start_input.j_, start_input.k_, res_input.i_, res_input.j_, res_input.k_,
		min_input.x_, min_input.y_, min_input.z_, max_input.x_, max_input.y_, max_input.z_);
}

void GridUniform3D::initialize(const GridUniform3D& grid_input)
{
	initialize(grid_input.i_start_, grid_input.j_start_, grid_input.k_start_, grid_input.i_res_, grid_input.j_res_, grid_input.k_res_,
		grid_input.x_min_, grid_input.y_min_, grid_input.z_min_, grid_input.x_max_, grid_input.y_max_, grid_input.z_max_);
}

void GridUniform3D::InitializeDualGrid(GridUniform3D& dual_grid_input)
{
	dual_grid_input.initialize(i_start_, j_start_, k_start_, i_res_ - 1, j_res_ - 1, k_res_ - 1, x_min_ + (T)0.5*dx_, y_min_ + (T)0.5*dy_, z_min_ + (T)0.5*dz_,
		x_max_ - (T)0.5*dx_, y_max_ - (T)0.5*dy_, z_max_ - (T)0.5*dz_);
}

void GridUniform3D::initialize(const GridUniform3D& grid_input, const T& resolution_scale)
{
	initialize(grid_input.i_start_, grid_input.j_start_, grid_input.k_start_,
		(int)((T)grid_input.i_res_*resolution_scale), (int)((T)grid_input.j_res_*resolution_scale), (int)((T)grid_input.k_res_*resolution_scale),
		grid_input.x_min_, grid_input.y_min_, grid_input.z_min_, grid_input.x_max_, grid_input.y_max_, grid_input.z_max_);
}

void GridUniform3D::initialize(const GridUniform3D& grid_input, const TV_INT& resolution_scale)
{
	initialize(grid_input.i_start_, grid_input.j_start_, grid_input.k_start_,
		grid_input.i_res_*resolution_scale.x_, grid_input.j_res_*resolution_scale.y_, grid_input.k_res_*resolution_scale.z_,
		grid_input.x_min_, grid_input.y_min_, grid_input.z_min_, grid_input.x_max_, grid_input.y_max_, grid_input.z_max_);
}

void GridUniform3D::initialize(const BOX& box_input, const int& max_res)
{
	const TV edge_length = box_input.GetEdgeLengths();

	if (edge_length.z_ >= edge_length.x_ && edge_length.z_ >= edge_length.y_)
	{
		const int z_res = max_res;
		const T dz = edge_length.z_ / (T)z_res;
		const int y_res = (int)ceil(edge_length.y_ / dz);
		const int x_res = (int)ceil(edge_length.x_ / dz);
		const TV corrected_max(box_input.x_min_ + (T)x_res*dz, box_input.y_min_ + (T)y_res*dz, box_input.z_min_ + (T)z_res*dz); // to make dx_ = dy_ = dz_

		initialize(0, 0, 0, x_res, y_res, z_res, box_input.x_min_, box_input.y_min_, box_input.z_min_, corrected_max.x_, corrected_max.y_, corrected_max.z_);
	}
	else if (edge_length.y_ >= edge_length.x_)
	{
		const int y_res = max_res;
		const T dy = edge_length.y_ / (T)y_res;
		const int z_res = (int)ceil(edge_length.z_ / dy);
		const int x_res = (int)ceil(edge_length.x_ / dy);
		const TV corrected_max(box_input.x_min_ + (T)x_res*dy, box_input.y_min_ + (T)y_res*dy, box_input.z_min_ + (T)z_res*dy); // to make dx_ = dy_ = dz_

		initialize(0, 0, 0, x_res, y_res, z_res, box_input.x_min_, box_input.y_min_, box_input.z_min_, corrected_max.x_, corrected_max.y_, corrected_max.z_);
	}
	else
	{
		const int x_res = max_res;
		const T dx = edge_length.x_ / (T)x_res;
		const int z_res = (int)ceil(edge_length.z_ / dx);
		const int y_res = (int)ceil(edge_length.y_ / dx);
		const TV corrected_max(box_input.x_min_ + (T)x_res*dx, box_input.y_min_ + (T)y_res*dx, box_input.z_min_ + (T)z_res*dx); // to make dx_ = dy_ = dz_

		initialize(0, 0, 0, x_res, y_res, z_res, box_input.x_min_, box_input.y_min_, box_input.z_min_, corrected_max.x_, corrected_max.y_, corrected_max.z_);
	}
}

void GridUniform3D::initialize(const BOX& box_input, const T& dx_input)
{
	const TV edge_length = box_input.GetEdgeLengths();
	const int x_res = (int)ceil(edge_length.x_ / dx_input);
	const int y_res = (int)ceil(edge_length.y_ / dx_input);
	const int z_res = (int)ceil(edge_length.z_ / dx_input);
	const TV corrected_max(box_input.x_min_ + (T)x_res*dx_input, box_input.y_min_ + (T)y_res*dx_input, box_input.z_min_ + (T)z_res*dx_input); // to make dx_ = dy_ = dz_

	initialize(0, 0, 0, x_res, y_res, z_res, box_input.x_min_, box_input.y_min_, box_input.z_min_, corrected_max.x_, corrected_max.y_, corrected_max.z_);
}

void GridUniform3D::initialize(const BOX& box_input, const int& i_res, const int& j_res, const int& k_res)
{
	initialize(0, 0, 0, i_res, j_res, k_res, box_input.x_min_, box_input.y_min_, box_input.z_min_, box_input.x_max_, box_input.y_max_, box_input.z_max_);
}

// 	void InitializeFromBlock(const SCRIPT_BLOCK& block)
// 	{
// 		const T resolution_scale = block.GetFloat("resolution_scale", (T)1);
// 		const TV_INT start = block.GetInt3("start_indices");
// 		const TV_INT res = block.GetInt3("base_grid_resolution");
// 		const TV_INT res_scaled((int)((T)res.x_*resolution_scale), (int)((T)res.y_*resolution_scale), (int)((T)res.z_*resolution_scale));
// 		const TV min = block.GetVector3("base_grid_min");
// 		const TV max = block.GetVector3("base_grid_max");
// 
// 		Initialize(start.i_, start.j_, start.k_, res_scaled.i_, res_scaled.j_, res_scaled.k_, min.i_, min.j_, min.k_, max.i_, max.j_, max.k_);
// 	}

void GridUniform3D::scaleRes(const T scale)
{
	initialize(i_start_, j_start_, k_start_, (int)((T)i_res_*scale), (int)((T)j_res_*scale), (int)((T)k_res_*scale), x_min_, y_min_, z_min_, x_max_, y_max_, z_max_);
}

GridUniform3D GridUniform3D::getPartialGrid(const BOX_3D<int>& partial_ix_box) const
{
	const TV partial_min = getCellMin(partial_ix_box.i_start_, partial_ix_box.j_start_, partial_ix_box.k_start_);
	const TV partial_max = getCellMin(partial_ix_box.i_end_, partial_ix_box.j_end_, partial_ix_box.k_end_);

	return GridUniform3D(partial_ix_box.i_start_, partial_ix_box.j_start_, partial_ix_box.k_start_, partial_ix_box.i_end_ - partial_ix_box.i_start_ + 1, partial_ix_box.j_end_ - partial_ix_box.j_start_ + 1,
		partial_ix_box.k_end_ - partial_ix_box.k_start_ + 1, partial_min.x_, partial_min.y_, partial_min.z_, partial_max.x_, partial_max.y_, partial_max.z_);
}

TV_INT GridUniform3D::getClampedIndex(const int& i, const int& j, const int& k) const
{
	return TV_INT(CLAMP(i, i_start_, i_end_), CLAMP(j, j_start_, j_end_), CLAMP(k, k_start_, k_end_));
}

 TV_INT GridUniform3D::getClampedIndex(const TV_INT& ix) const
{
	return TV_INT(CLAMP(ix.i_, i_start_, i_end_), CLAMP(ix.j_, j_start_, j_end_), CLAMP(ix.k_, k_start_, k_end_));
}

 void GridUniform3D::getClampStartEndIndices(int& l_start, int& m_start, int& n_start, int& l_end, int& m_end, int& n_end)
{
	l_start = MAX2(l_start, i_start_);
	m_start = MAX2(m_start, j_start_);
	n_start = MAX2(n_start, k_start_);

	l_end = MIN2(l_end, i_end_);
	m_end = MIN2(m_end, j_end_);
	n_end = MIN2(n_end, k_end_);
}

int GridUniform3D::get1DIndex(const int& i, const int& j, const int& k) const
{
	return (i - i_start_) + (j - j_start_)*i_res_ + (k - k_start_)*i_res_*j_res_;
}

TV_INT GridUniform3D::get3DIndex(const int& index_1d) const
{
	const int n = index_1d / ij_res_;
	const int n_ij_res = n*ij_res_;
	const int m = (index_1d - n_ij_res) / i_res_;
	const int l = index_1d - i_res_*m - n_ij_res;

	return TV_INT(l, m, n);
}

 int GridUniform3D::getClampedIndexI(const int& i) const
{
	return CLAMP(i, i_start_, i_end_);
}

 int GridUniform3D::getClampedIndexJ(const int& j) const
{
	return CLAMP(j, j_start_, j_end_);
}

 int GridUniform3D::getClampedIndexK(const int& k) const
{
	return CLAMP(k, k_start_, k_end_);
}

 TV GridUniform3D::getCellCenter(const int& i, const int& j, const int& k) const
{
	return TV(x_min_ + ((T)0.5 + (T)(i - i_start_))*dx_, y_min_ + ((T)0.5 + (T)(j - j_start_))*dy_, z_min_ + ((T)0.5 + (T)(k - k_start_))*dz_);
}

 TV GridUniform3D::getCellCenter(const TV_INT& ix) const
{
	return TV(x_min_ + ((T)0.5 + (T)(ix.i_ - i_start_))*dx_, y_min_ + ((T)0.5 + (T)(ix.j_ - j_start_))*dy_, z_min_ + ((T)0.5 + (T)(ix.k_ - k_start_))*dz_);
}

void GridUniform3D::getCenter(const int& i, const int& j, const int& k, TV& position)
{
	position.x_ = x_min_ + ((T)0.5 + (T)(i - i_start_))*dx_;
	position.y_ = y_min_ + ((T)0.5 + (T)(j - j_start_))*dy_;
	position.z_ = z_min_ + ((T)0.5 + (T)(k - k_start_))*dz_;
}

void GridUniform3D::getCenter(const TV_INT& ix, TV& position) const
{
	position.x_ = x_min_ + ((T)0.5 + (T)(ix.i_ - i_start_))*dx_;
	position.y_ = y_min_ + ((T)0.5 + (T)(ix.j_ - j_start_))*dy_;
	position.z_ = z_min_ + ((T)0.5 + (T)(ix.k_ - k_start_))*dz_;
}

 TV_INT GridUniform3D::getCell(const TV& position) const // return a cell index contains the position
{
	return TV_INT(i_start_ + (int)(floor((position.x_ - x_min_)*one_over_dx_)), j_start_ + (int)(floor((position.y_ - y_min_)*one_over_dy_)), k_start_ + (int)(floor((position.z_ - z_min_)*one_over_dz_)));
}

void GridUniform3D::getCell(const TV& position, TV_INT& index) const // return a cell index contains the position
{
	index.i_ = i_start_ + (int)(floor((position.x_ - x_min_)*one_over_dx_));
	index.j_ = j_start_ + (int)(floor((position.y_ - y_min_)*one_over_dy_));
	index.k_ = k_start_ + (int)(floor((position.z_ - z_min_)*one_over_dz_));
}

 void GridUniform3D::getCell(const TV& position, int& i, int& j, int& k) const // return a cell index contains the position
{
	i = i_start_ + (int)(floor((position.x_ - x_min_)*one_over_dx_));
	j = j_start_ + (int)(floor((position.y_ - y_min_)*one_over_dy_));
	k = k_start_ + (int)(floor((position.z_ - z_min_)*one_over_dz_));

	assert(i >= i_start_ && i <= i_end_);
	assert(j >= j_start_ && j <= j_end_);
	assert(k >= k_start_ && k <= k_end_);
}

 void GridUniform3D::getCell(const T& x, const T& y, const T& z, int& i, int& j, int& k) const // return a cell index contains the position
{
	i = i_start_ + (int)(floor((x - x_min_)*one_over_dx_));
	j = j_start_ + (int)(floor((y - y_min_)*one_over_dy_));
	k = k_start_ + (int)(floor((z - z_min_)*one_over_dz_));

	assert(i >= i_start_ && i <= i_end_);
	assert(j >= j_start_ && j <= j_end_);
	assert(k >= k_start_ && k <= k_end_);
}

 TV_INT GridUniform3D::getClampedCell(const TV& position) const
{
	TV_INT index(i_start_ + (int)(floor((position.x_ - x_min_)*one_over_dx_)),
		j_start_ + (int)(floor((position.y_ - y_min_)*one_over_dy_)),
		k_start_ + (int)(floor((position.z_ - z_min_)*one_over_dz_)));

	index.i_ = CLAMP(index.i_, i_start_, i_end_);
	index.j_ = CLAMP(index.j_, j_start_, j_end_);
	index.k_ = CLAMP(index.k_, k_start_, k_end_);

	return index;
}

 void GridUniform3D::getClampedCell(const TV& position, TV_INT& index) const
{
	index.assign(i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_),
		j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_),
		k_start_ + (int)floor((position.z_ - z_min_)*one_over_dz_));

	index.i_ = CLAMP(index.i_, i_start_, i_end_);
	index.j_ = CLAMP(index.j_, j_start_, j_end_);
	index.k_ = CLAMP(index.k_, k_start_, k_end_);
}

 void GridUniform3D::getClampedCell(const TV& position, int& i, int& j, int& k) const
{
	i = i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_);
	j = j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_);
	k = k_start_ + (int)floor((position.z_ - z_min_)*one_over_dz_);

	assert(i_start_ < i_end_);
	assert(j_start_ < j_end_);
	assert(k_start_ < k_end_);

	if (i < i_start_) i = i_start_;
	else if (i > i_end_) i = i_end_;

	if (j < j_start_) j = j_start_;
	else if (j > j_end_) j = j_end_;

	if (k < k_start_) k = k_start_;
	else if (k > k_end_) k = k_end_;
}

 TV_INT GridUniform3D::getLeftBottomCell(const TV& position) const //TODO check performance compared to no return value functions.
{
	return TV_INT(i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_ - (T)0.5),
		j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_ - (T)0.5),
		k_start_ + (int)floor((position.z_ - z_min_)*one_over_dz_ - (T)0.5));
}

 TV_INT GridUniform3D::getRightTopCell(const TV& position) const
{
	return TV_INT(i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_ + (T)0.5),
		j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_ + (T)0.5),
		k_start_ + (int)floor((position.z_ - z_min_)*one_over_dz_ + (T)0.5));
}

void GridUniform3D::getLeftBottomCell(const TV& position, int& i, int& j, int& k) const
{
	i = i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_ - (T)0.5);
	j = j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_ - (T)0.5);
	k = k_start_ + (int)floor((position.z_ - z_min_)*one_over_dz_ - (T)0.5);
}

void GridUniform3D::getLeftBottomCell(const T& x, const T& y, const T& z, int& i, int& j, int& k) const
{
	i = i_start_ + (int)floor((x - x_min_)*one_over_dx_ - (T)0.5);
	j = j_start_ + (int)floor((y - y_min_)*one_over_dy_ - (T)0.5);
	k = k_start_ + (int)floor((z - z_min_)*one_over_dz_ - (T)0.5);
}

void GridUniform3D::getLeftBottomCell(const TV& position, TV_INT& ix) const
{
	return getLeftBottomCell(position, ix.i_, ix.j_, ix.k_);
}

bool GridUniform3D::isInside(const TV& position) const
{
	if (position.x_ <= x_min_) return false;
	else if (position.x_ >= x_max_) return false;
	else if (position.y_ <= y_min_) return false;
	else if (position.y_ >= y_max_) return false;
	else if (position.z_ <= z_min_) return false;
	else if (position.z_ >= z_max_) return false;
	return true;
}

bool GridUniform3D::isInside(const TV& position, const T& width) const
{
	if (position.x_ <= x_min_ + width) return false;
	else if (position.x_ >= x_max_ - width) return false;
	else if (position.y_ <= y_min_ + width) return false;
	else if (position.y_ >= y_max_ - width) return false;
	else if (position.z_ <= z_min_ + width) return false;
	else if (position.z_ >= z_max_ - width) return false;
	return true;
}

bool GridUniform3D::isInside(const int& i, const int& j, const int& k) const
{
	if (i < i_start_) return false;
	else if (i > i_end_) return false;
	else if (j < j_start_) return false;
	else if (j > j_end_) return false;
	else if (k < k_start_) return false;
	else if (k > k_end_) return false;
	else return true;
}

bool GridUniform3D::isInside(const TV_INT& ix) const
{
	if (ix.i_ < i_start_) return false;
	else if (ix.i_ > i_end_) return false;
	else if (ix.j_ < j_start_) return false;
	else if (ix.j_ > j_end_) return false;
	else if (ix.k_ < k_start_) return false;
	else if (ix.k_ > k_end_) return false;
	else return true;
}

bool GridUniform3D::isInside(const TV_INT& ix, const int& inner_width) const
{
	if (ix.i_ < i_start_ + inner_width) return false;
	else if (ix.i_ > i_end_ - inner_width) return false;
	else if (ix.j_ < j_start_ + inner_width) return false;
	else if (ix.j_ > j_end_ - inner_width) return false;
	else if (ix.k_ < k_start_ + inner_width) return false;
	else if (ix.k_ > k_end_ - inner_width) return false;
	else return true;
}

bool GridUniform3D::isInsideI(const int & i)	// check if i_start_ <= i <= i_end_
{
	if (i < i_start_) return false;
	else if (i > i_end_) return false;
	return true;
}

bool GridUniform3D::isInsideJ(const int & j)	// check if j_start_ <= j <= j_end_
{
	if (j < j_start_) return false;
	else if (j > j_end_) return false;
	return true;
}

bool GridUniform3D::isInsideK(const int & k)	// check if k_start_ <= k <= k_end_
{
	if (k < k_start_) return false;
	else if (k > k_end_) return false;
	return true;
}

GridUniform3D GridUniform3D::getCellGrid() const	// returns cell grid when this is a node grid
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_ - 1, j_res_ - 1, k_res_ - 1, x_min_ + (T)0.5*dx_, y_min_ + (T)0.5*dy_, z_min_ + (T)0.5*dz_, x_max_ - (T)0.5*dx_, y_max_ - (T)0.5*dy_, z_max_ - (T)0.5*dz_);
}

GridUniform3D GridUniform3D::getUFaceGrid() const
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_ + 1, j_res_, k_res_, x_min_ - (T)0.5*dx_, y_min_, z_min_, x_max_ + (T)0.5*dx_, y_max_, z_max_);
}

GridUniform3D GridUniform3D::getVFaceGrid() const
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_, j_res_ + 1, k_res_, x_min_, y_min_ - (T)0.5*dy_, z_min_, x_max_, y_max_ + (T)0.5*dy_, z_max_);
}

GridUniform3D GridUniform3D::getWFaceGrid() const
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_, j_res_, k_res_ + 1, x_min_, y_min_, z_min_ - (T)0.5*dz_, x_max_, y_max_, z_max_ + (T)0.5*dz_);
}

GridUniform3D GridUniform3D::getXEdgeGrid() const
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_, j_res_ + 1, k_res_ + 1, x_min_, y_min_ - (T)0.5*dy_, z_min_ - (T)0.5*dz_, x_max_, y_max_ + (T)0.5*dy_, z_max_ + (T)0.5*dz_);
}

GridUniform3D GridUniform3D::getYEdgeGrid() const
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_ + 1, j_res_, k_res_ + 1, x_min_ - (T)0.5*dx_, y_min_, z_min_ - (T)0.5*dz_, x_max_ + (T)0.5*dx_, y_max_, z_max_ + (T)0.5*dz_);
}

GridUniform3D GridUniform3D::getZEdgeGrid() const
{
	return GridUniform3D(i_start_, j_start_, k_start_, i_res_ + 1, j_res_ + 1, k_res_, x_min_ - (T)0.5*dx_, y_min_ - (T)0.5*dy_, z_min_, x_max_ + (T)0.5*dx_, y_max_ + (T)0.5*dy_, z_max_);
}

GridUniform3D GridUniform3D::getEnlarged(const int& width) const
{
	return GridUniform3D(i_start_ - width, j_start_ - width, k_start_ - width, i_res_ + 2 * width, j_res_ + 2 * width, k_res_ + 2 * width,
		x_min_ - (T)width*dx_, y_min_ - (T)width*dy_, z_min_ - (T)width*dz_, x_max_ + (T)width*dx_, y_max_ + (T)width*dy_, z_max_ + (T)width*dz_);
}

void GridUniform3D::enlarge(const int& width)
{
	initialize(i_start_ - width, j_start_ - width, k_start_ - width, i_res_ + 2 * width, j_res_ + 2 * width, k_res_ + 2 * width,
		x_min_ - (T)width*dx_, y_min_ - (T)width*dy_, z_min_ - (T)width*dz_, x_max_ + (T)width*dx_, y_max_ + (T)width*dy_, z_max_ + (T)width*dz_);
}

void GridUniform3D::operator = (const GridUniform3D& grid_input)
{
	initialize(grid_input);
}

void GridUniform3D::translate(const TV& deviation)
{
	x_min_ += deviation.x_;
	y_min_ += deviation.y_;
	z_min_ += deviation.z_;

	x_max_ += deviation.x_;
	y_max_ += deviation.y_;
	z_max_ += deviation.z_;
}

void GridUniform3D::scaleToBox(const BOX& box)
{
	TV_INT ix_min = getClampedCell(box.GetMin());
	TV_INT ix_max = getClampedCell(box.GetMax());

	ix_min = getClampedIndex(ix_min - TV_INT(2, 2, 2));
	ix_max = getClampedIndex(ix_max + TV_INT(2, 2, 2));

	TV cell_center_min = getCellCenter(ix_min);
	TV cell_center_max = getCellCenter(ix_max);

	i_start_ = ix_min.i_;
	j_start_ = ix_min.j_;
	k_start_ = ix_min.k_;
	i_end_ = ix_max.i_;
	j_end_ = ix_max.j_;
	k_end_ = ix_max.k_;

	x_min_ = cell_center_min.x_ - dx_*(T)0.5;
	y_min_ = cell_center_min.y_ - dy_*(T)0.5;
	z_min_ = cell_center_min.z_ - dz_*(T)0.5;
	x_max_ = cell_center_max.x_ + dx_*(T)0.5;
	y_max_ = cell_center_max.y_ + dy_*(T)0.5;
	z_max_ = cell_center_max.z_ + dz_*(T)0.5;
}

int GridUniform3D::getNumAllCells() const
{
	return i_res_ * j_res_ * k_res_;
}

bool GridUniform3D::operator == (const GridUniform3D& grid)	// for debugging only (not optimized)
{
	for (int d = 0; d < 3; ++d)
	{
		if (res_[d] != grid.res_[d]) return false;
		if (ix_start_[d] != grid.ix_start_[d]) return false;
		if (ix_end_[d] != grid.ix_end_[d]) return false;
		if (min_[d] != grid.min_[d]) return false;
		if (max_[d] != grid.max_[d]) return false;
	}
	return true;
}

TV GridUniform3D::getMin() const
{
	return TV(x_min_, y_min_, z_min_);
}

TV GridUniform3D::getMax() const
{
	return TV(x_max_, y_max_, z_max_);
}

BOX_3D<int> GridUniform3D::getIXBox() const
{
	return BOX_3D<int>(i_start_, j_start_, k_start_, i_end_, j_end_, k_end_);
}

BOX_3D<T> GridUniform3D::getMimMax() const
{
	return BOX_3D<T>(getMin(), getMax());
}

TV GridUniform3D::getCellMin(const int& i, const int& j, const int& k) const
{
	const TV center = getCellCenter(i, j, k);

	return TV(center.x_ - 0.5f*dx_, center.y_ - 0.5f*dy_, center.z_ - 0.5f*dz_);
}

TV GridUniform3D::getCellMax(const int& i, const int& j, const int& k) const
{
	const TV center = getCellCenter(i, j, k);

	return TV(center.x_ + 0.5f*dx_, center.y_ + 0.5f*dy_, center.z_ + 0.5f*dz_);
}

BOX_3D<T> GridUniform3D::getCellMinMax(const int& i, const int& j, const int& k) const
{
	return BOX_3D<T>(getCellMin(i, j, k), getCellMax(i, j, k));
}

void GridUniform3D::extend(const int width)
{
	initialize(i_start_ - width, j_start_ - width, k_start_ - width, i_res_ + 2 * width, j_res_ + 2 * width, k_res_ + 2 * width,
		x_min_ - (T)width*dx_, y_min_ - (T)width*dy_, z_min_ - (T)width*dz_, x_max_ + (T)width*dx_, y_max_ + (T)width*dy_, z_max_ + (T)width*dz_);
}

void GridUniform3D::splitInMaxDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids)
{
	// give priority to z and y direction in order to utilize the memory coherency
	if(num_threads <= k_res_) splitInZDirection(num_threads, partial_grids);
	else if(k_res_ >= j_res_ && k_res_ >= i_res_) splitInZDirection(num_threads, partial_grids);
	else if(j_res_ >= i_res_) splitInYDirection(num_threads, partial_grids);
	else splitInXDirection(num_threads, partial_grids);
}

void GridUniform3D::splitInXDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids)
{
	partial_grids.initialize(num_threads);

	const int quotient = i_res_/num_threads;
	const int remainder = i_res_%num_threads;

	int i_start_p = i_start_;
	T x_min_p = x_min_;

	for(int t = 0; t < num_threads; t++)
	{
		const int i_res_p = (t < remainder) ? (quotient + 1) : (quotient);
		const T x_max_p = dx_*(T)i_res_p + x_min_p;
		partial_grids[t].initialize(i_start_p, j_start_, k_start_, i_res_p, j_res_, k_res_, x_min_p, y_min_, z_min_, x_max_p, y_max_, z_max_);
		i_start_p += i_res_p;
		x_min_p = x_max_p;
	}
}

void GridUniform3D::splitInYDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids)
{
	partial_grids.initialize(num_threads);

	const int quotient = j_res_/num_threads;
	const int remainder = j_res_%num_threads;

	int j_start_p = j_start_;
	T y_min_p = y_min_;

	for(int t = 0; t < num_threads; t++)
	{
		const int j_res_p = (t < remainder) ? (quotient + 1) : (quotient);
		const T y_max_p = dy_*(T)j_res_p + y_min_p;
		partial_grids[t].initialize(i_start_, j_start_p, k_start_, i_res_, j_res_p, k_res_, x_min_, y_min_p, z_min_, x_max_, y_max_p, z_max_);
		j_start_p += j_res_p;
		y_min_p = y_max_p;
	}
}

void GridUniform3D::splitInZDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids)
{
	partial_grids.initialize(num_threads);

	const int quotient = k_res_/num_threads;
	const int remainder = k_res_%num_threads;
	
	int k_start_p = k_start_;
	T z_min_p = z_min_;

	for(int t = 0; t < num_threads; t++)
	{
		const int k_res_p = (t < remainder) ? (quotient + 1) : (quotient);
		const T z_max_p = dz_*(T)k_res_p + z_min_p;
		partial_grids[t].initialize(i_start_, j_start_, k_start_p, i_res_, j_res_, k_res_p, x_min_, y_min_, z_min_p, x_max_, y_max_, z_max_p);
		k_start_p += k_res_p;
		z_min_p = z_max_p;
	}
}

void GridUniform3D::splitInZDirectionPartially(const int& num_threads, Array1D<GridUniform3D>& partial_grids)
{
	partial_grids.initialize(num_threads);

	int i_idx = i_end_ - i_start_ + 1;
	int j_idx = j_end_ - j_start_ + 1;
	int k_idx = k_end_ - k_start_ + 1;

	const int quotient  = k_idx / num_threads;
	const int remainder = k_idx % num_threads;
	
	int k_start_p = k_start_;
	T z_min_p = z_min_;

	for(int t = 0; t < num_threads; t++)
	{
		const int k_res_p = (t < remainder) ? (quotient + 1) : (quotient);

		const T z_max_p = dz_*(T)k_res_p + z_min_p;

		partial_grids[t].initialize(i_start_, j_start_, k_start_p, i_idx, j_idx, k_res_p, x_min_, y_min_, z_min_p, x_max_, y_max_, z_max_p);

		k_start_p += k_res_p;
		z_min_p = z_max_p;
	}
}

int GridUniform3D::recommendMaxMultigridLevel(const int& lowest_level_res)
{
	int min_res = MIN3(i_res_, j_res_, k_res_);
	int max_level = 0;
	while(true)// find max_evel
	{
		max_level++;
		min_res /= 2;
		if(min_res < lowest_level_res) break;
	}

	return max_level;
}


T	GridUniform3D::getTrilinearWeight(const TV& deviation) const
{
	// assert ( (T)1 - ABS(deviation.x_) * one_over_dx_ >= (T)0)
	// assert ( (T)1 - ABS(deviation.x_) * one_over_dx_ <= (T)1)
	return MAX2((T)0, ((T)1 - ABS(deviation.x_) * one_over_dx_) * ((T)1 - ABS(deviation.y_) * one_over_dy_) * ((T)1 - ABS(deviation.z_) * one_over_dz_));
}