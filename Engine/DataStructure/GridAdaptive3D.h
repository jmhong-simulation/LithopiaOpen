// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GridUniform3D.h"
#include "Array3D.h"
#include "Geometry/BOX_3D.h"
#include "Geometry/StaticTriangle.h"
#include "Geometry/StaticTriangularSurface.h"

const int RES_TABLE[5] = { 1, 8, 64, 512, 4096 };

class GridAdaptive3D
{
public:
	GridUniform3D grid_;
	Array3D<bool> interfacial_;
	Array3D<int>  level_;			// refinement level of the cell
	Array3D<int>  first_ix_;				// refinement level of the cell

	int max_level_;					// max refinement level

	Array1D<bool> inter_;				

public:
	GridAdaptive3D()
		: max_level_(3)
	{}

	~GridAdaptive3D()
	{}

	void Initialize(const BOX_3D<T>& box_input, const int& max_res)
	{
		grid_.initialize(box_input, max_res);
		grid_.extend(3);	//ghost cells

		grid_.initializeCenterArray(interfacial_);
		interfacial_.AssignAllValues(false);

		grid_.initializeCenterArray(level_);
		level_.AssignAllValues(0);

		grid_.initializeCenterArray(first_ix_);
		first_ix_.AssignAllValues(-1);
	}

	void CheckIntersection(const StaticTriangularSurface& surface_)
	{
		const T th = grid_.dx_ * 1.42f;

		for (int tri = 0; tri < surface_.triangles_.num_elements_; ++tri)
		{
			const BOX_3D<T> tri_aabb = surface_.tri_ops_.getAABB(tri);
			const BOX_3D<int> ix_aabb = BOX_3D<int>(grid_.getCell(tri_aabb.GetMin()), grid_.getCell(tri_aabb.GetMax()));

			const StaticTriangle tri_geo(surface_.tri_ops_.getVertexPosition(tri, 0), surface_.tri_ops_.getVertexPosition(tri, 1), surface_.tri_ops_.getVertexPosition(tri, 2));

			for (int k = ix_aabb.k_start_; k <= ix_aabb.k_end_; ++k)
			for (int j = ix_aabb.j_start_; j <= ix_aabb.j_end_; ++j)
			for (int i = ix_aabb.i_start_; i <= ix_aabb.i_end_; ++i)
			{
				if (tri_geo.getDistance(grid_.getCellCenter(i,j,k)) < th)
				{
					interfacial_(i, j, k) = true;
					level_(i, j, k) = max_level_;
				}
			}
		}

		int num_all_cells = 0;
		for (int i = 0; i < level_.ijk_res_; ++i)
		{
			first_ix_.values_[i] = num_all_cells;

			num_all_cells += RES_TABLE[level_.values_[i]];
		}

		inter_.initialize(num_all_cells);
		inter_.assignAllValues(false);

		for (int tri = 0; tri < surface_.triangles_.num_elements_; ++tri)
		{
			const BOX_3D<T> tri_aabb = surface_.tri_ops_.getAABB(tri);
			const BOX_3D<int> ix_aabb = BOX_3D<int>(grid_.getCell(tri_aabb.GetMin()), grid_.getCell(tri_aabb.GetMax()));

			const StaticTriangle tri_geo(surface_.tri_ops_.getVertexPosition(tri, 0), surface_.tri_ops_.getVertexPosition(tri, 1), surface_.tri_ops_.getVertexPosition(tri, 2));

			for (int k = ix_aabb.k_start_; k <= ix_aabb.k_end_; ++k)
			for (int j = ix_aabb.j_start_; j <= ix_aabb.j_end_; ++j)
			for (int i = ix_aabb.i_start_; i <= ix_aabb.i_end_; ++i)
			{
				const int l = level_(i, j, k);
				if (l == 0)
				{
					if (tri_geo.getDistance(grid_.getCellCenter(i, j, k)) < th)
					{
						interfacial_(i, j, k) = true;
						level_(i, j, k) = max_level_;
					}
				}
				else // if l > 0
				{
					GridUniform3D subgrid;
					subgrid.initialize(grid_.getCellMinMax(i, j, k), POW_OF_TWO(l));

					const T th_sub = th / (T)POW_OF_TWO(l);

					for (int n = subgrid.k_start_; n <= subgrid.k_end_; ++n)
					for (int m = subgrid.j_start_; m <= subgrid.j_end_; ++m)
					for (int l = subgrid.i_start_; l <= subgrid.i_end_; ++l)
					{
						if (tri_geo.getDistance(subgrid.getCellCenter(l, m, n)) < th_sub)
						{
							inter_.values_[first_ix_(i, j, k) + subgrid.get1DIndex(l, m, n)] = true;
						}
					}
				}
			}
		}
	}

};