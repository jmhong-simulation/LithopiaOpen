// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "LITHOPHANE.h"
#include "UV_TO_XYZ.h"
#include "../DataStructure/Vector2D.h"
#include "../DataStructure/GridUniform2D.h"
#include "../Geometry/RAY.h"
#include "../Operations/SMOOTHING_UNIFORM_2D.h"

#include <iostream>

void LITHOPHANE::InitializeCylinder(const IMAGE_2D& image, StaticTriangularSurface& surface, const float inner_radius, const float base_thickness, const float front_print_thickness, const int smoothing_repeat)
{
	// use num_side = 0 to generate a cylinder
	InitializePrism(image, surface, 0, inner_radius, base_thickness, front_print_thickness, smoothing_repeat);
}

void LITHOPHANE::InitializeCircularTruncatedCone(const IMAGE_2D& image, StaticTriangularSurface& surface, const int num_sides, const float upper_inner_radius, const float bottom_inner_radius, const float base_thickness, const float front_print_thickness, const int smoothing_repeat, const int thickness_ceiling, const int thickness_floor)
{
	// use num_side = 0 to generate a cylinder
//	InitializePrism(image, surface, 0, upper_radius, bottom_radius, base_thickness, front_print_thickness, smoothing_repeat);

	const T upper_outer_cylinder_radius = upper_inner_radius + base_thickness + front_print_thickness;
	const T bottom_outer_cylinder_radius = bottom_inner_radius + base_thickness + front_print_thickness;

	const float height = 2.0f*PI*(upper_outer_cylinder_radius + bottom_outer_cylinder_radius)*0.5f;		// use average radius to calculate height
	const float width = height / (T)image.res_y_ * (T)image.res_x_;

	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		height_field_node(i, j) = (image.data_.getIClampedJRepeated(i - 1, j - 1).GetReversedGray() + image.data_.getIClampedJRepeated(i, j - 1).GetReversedGray()
			+ image.data_.getIClampedJRepeated(i, j - 1).GetReversedGray() + image.data_.getIClampedJRepeated(i, j).GetReversedGray()) * 0.25f;
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIClampJRepeat(grid, height_field_node, smoothing_repeat);

	// cylinder shares one j-row
	surface.vertex_positions_.initialize(grid.GetNumAllNodes() - (grid.i_res_ + 1) + grid.GetNumAllNodes() - (grid.i_res_ + 1));

	// define vertices
	int ix = 0;

	// vertices of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		if (num_sides <= 0)	// cylinder
		{
			const T upper_outer_radius = upper_inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;
			const T bottom_outer_radius = bottom_inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;

			surface.vertex_positions_[ix++] = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), upper_outer_radius, bottom_outer_radius, width);
		}
		else // prism
		{
			const T upper_outer_radius = upper_inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;
			const T bottom_outer_radius = bottom_inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;

			surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZ(grid.GetNodeUV(i, j), num_sides, upper_outer_radius, bottom_outer_radius, width);
		}
	}

	// all vertices of inner cylinder plane
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
	{
	//		TV v = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), inner_cylinder_radius, par.width_);
	//		TV v = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), inner_cylinder_radius - height(i, j)*par.back_print_thickness_, par.width_);
	//		if (v.x_ < base_thickness) v = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), par.inner_candle_radius_, par.width_);
//			surface.vertex_positions_[ix++] = v;

			if (num_sides <= 0)	// cylinder
			{
				const T upper_outer_radius = upper_inner_radius;
				const T bottom_outer_radius = bottom_inner_radius;

				TV vp = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), upper_outer_radius, bottom_outer_radius, width);

				if (i >= height_field_node.i_res_ - thickness_ceiling)
				{
					if (i == height_field_node.i_res_ - thickness_ceiling)
						vp = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i-1, j), upper_outer_radius, bottom_outer_radius, width);	// one dx down

					vp.y_ = (T)0;
					vp.z_ = (T)0;
				}

				surface.vertex_positions_[ix++] = vp;
			}
			else // prism
			{
				const T upper_outer_radius = upper_inner_radius;
				const T bottom_outer_radius = bottom_inner_radius;

				surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZ(grid.GetNodeUV(i, j), num_sides, upper_outer_radius, bottom_outer_radius, width);
			}	
	}

	// reducing number of inner triangles
// 	const int v_start_bottom = ix;	// start index of bottom line vertices of inner cylinder
// 	{const int i = 0;
// 	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
// 	{
// 		if (num_sides <= 0)
// 		{
// 			TV vp = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), upper_inner_radius, bottom_inner_radius, width);
// 
// 			surface.vertex_positions_[ix++] = vp;		// cylinder
// 		}
// 		else
// 		{
// 			surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZ(grid.GetNodeUV(i, j), num_sides, upper_inner_radius, bottom_inner_radius, width);			// prism
// 		}
// 	}}
// 
// 	const int v_start_upper = ix;	// start index of upper line vertices of inner cylinder
// 	{const int i = height_field_node.i_res_ - 1;
// 	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
// 	{
// 		if (num_sides <= 0){
// 
// 			TV vp = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), upper_inner_radius, bottom_inner_radius, width);
// 
// 			surface.vertex_positions_[ix++] = vp;		// cylinder
// 		}
// 		else surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZ(grid.GetNodeUV(i, j), num_sides, upper_inner_radius, bottom_inner_radius, width);			// prism
// 	}
// 	}

	LinkedArray<TV_INT> tri;

	// triangles of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	for (int i = 0; i < height_field_node.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, j + 1));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i + 1, j + 1));
	}

	int j = height_field_node.j_res_ - 2;
	for (int i = 0; i < height_field_node.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, 0));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, 0), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i + 1, 0));
	}

	const int offset = height_field_node.i_res_*height_field_node.j_res_ - height_field_node.i_res_;		// one j row is shared

	// triangles of inner cylinder
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
		for (int i = 0; i < height_field_node.i_res_; i++)
	{
		if (i >= height_field_node.i_res_ - thickness_ceiling)
		{
			continue;	// do not use inner ceiling side triangles
		}

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i + 1, j + 1) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
	}

	j = height_field_node.j_res_ - 2;
	for (int i = 0; i < height_field_node.i_res_; i++)
	{
		if (i >= height_field_node.i_res_ - thickness_ceiling)
		{
			continue;	// do not use inner ceiling side triangles
		}

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, 0) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, 0) + offset, height_field_node.get1DIndex(i + 1, 0) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
	}

/*
	for (int j = 0; j < height_field_node.j_res_ - 2; ++j)
	{
		tri.PushBack() = TV_INT(v_start_bottom + j, v_start_bottom + 1 + j, v_start_upper + j);
		tri.PushBack() = TV_INT(v_start_bottom + 1 + j, v_start_upper + j + 1, v_start_upper + j);
	}
	{
		int j = height_field_node.j_res_ - 2;
		tri.PushBack() = TV_INT(v_start_bottom + j, v_start_bottom, v_start_upper + j);
		tri.PushBack() = TV_INT(v_start_bottom, v_start_upper, v_start_upper + j);
	}
	*/

	//	left wall
	int i = 0;

	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1));
	}

	{
		int j = height_field_node.j_res_ - 2;

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, 0) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, 0) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, 0));
	}

	/*
	{int i = 0;
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(v_start_bottom + j, height_field_node.get1DIndex(i, j), v_start_bottom + j + 1);
		tri.PushBack() = TV_INT(v_start_bottom + j + 1, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1));
	}
	

	int j = height_field_node.j_res_ - 2;
	tri.PushBack() = TV_INT(v_start_bottom + j, height_field_node.get1DIndex(i, j), v_start_bottom);
	tri.PushBack() = TV_INT(v_start_bottom, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, 0)); }
	*/

	// right wall
	i = height_field_node.i_res_ - 1;

	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i, j));
	}

	{
		int j = height_field_node.j_res_ - 2;

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, 0) + offset, height_field_node.get1DIndex(i, j));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, 0) + offset, height_field_node.get1DIndex(i, 0), height_field_node.get1DIndex(i, j));
	}

	/*

	{int i = height_field_node.i_res_ - 1;
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(v_start_upper + j, v_start_upper + j + 1, height_field_node.get1DIndex(i, j));
		tri.PushBack() = TV_INT(v_start_upper + j + 1, height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i, j));
	}

	{int j = height_field_node.j_res_ - 2;
	tri.PushBack() = TV_INT(v_start_upper + j, v_start_upper, height_field_node.get1DIndex(i, j));
	tri.PushBack() = TV_INT(v_start_upper, height_field_node.get1DIndex(i, 0), height_field_node.get1DIndex(i, j)); }}
	*/

	tri.CopyToArray(surface.triangles_);

	// rotate cylinder
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV old = surface.vertex_positions_[vix];

		surface.vertex_positions_[vix] = TV(old.y_, old.z_, old.x_);
	}

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();
}

void LITHOPHANE::InitializePlane(const Array2D<T>& height_field_node, StaticTriangularSurface& surface, const float base_depth, const float print_depth, const float print_width)
{
	int res_x = height_field_node.i_res_, res_y = height_field_node.j_res_;

	GridUniform2D grid(height_field_node.i_start_, height_field_node.j_start_, res_x, res_y, 0, 0, print_width, print_width / (T)res_x*(T)res_y);

// 	Array2D<T> height_field_node;
// 	grid.InitializeNodeArray(height_field_node);

//	height_field_node.assignAllValues(0.0f);

// 	for (int j = 0; j < grid.j_res_; ++j)
// 		for (int i = 0; i < grid.i_res_; ++i)
// 			height_field_node(i, height_field_node.j_res_ - 1 - j) = image.data_(i, j).GetReversedGray();

// 	SMOOTHING_UNIFORM_2D smoothing;
// 	smoothing.SmoothDirichletIJ(grid, height_field_node, smoothing_repeat);

	surface.vertex_positions_.initialize(grid.GetNumAllNodes() + grid.GetNumAllNodes());

	// vertices of upper plane
	int ix = 0;
	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			const TV2 uv = grid.GetNode(i, j);

			surface.vertex_positions_[ix++] = TV(uv.u_, uv.v_, base_depth + height_field_node(i, j)*print_depth);
		}

	// vertices of bottom plane
	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			const TV2 uv = grid.GetNode(i, j);

			surface.vertex_positions_[ix++] = TV(uv.u_, uv.v_, 0.0f);
		}

	LinkedArray<TV_INT> tri;

	// triangles of upper plane
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
		for (int i = 0; i < height_field_node.i_res_ - 1; i++)
		{
			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, j + 1));
			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i + 1, j + 1));
		}

	const int offset = height_field_node.i_res_*height_field_node.j_res_;

	// triangles of bottom plane
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
		for (int i = 0; i < height_field_node.i_res_ - 1; i++)
		{
			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i + 1, j + 1) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
		}

	// left wall
	{int i = 0;
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1));
	}}

	// right wall
	{int i = height_field_node.i_res_ - 1;
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i, j));
	}}

	{int j = 0;
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i + 1, j) + offset, height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, j));
	}}

	{int j = height_field_node.j_res_ - 1;
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j) + offset, height_field_node.get1DIndex(i, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i + 1, j) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j));
	}}

	tri.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();
}

void LITHOPHANE::InitializePlane(const IMAGE_2D& image, StaticTriangularSurface& surface, const float base_depth, const float print_depth, const float print_width, const int smoothing_repeat)
{
	int res_x = image.res_x_, res_y = image.res_y_;

	GridUniform2D grid(0, 0, res_x, res_y, 0, 0, print_width, print_width/(T)res_x*(T)res_y);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < grid.j_res_; ++j)
	for (int i = 0; i < grid.i_res_; ++i)
		height_field_node(i, height_field_node.j_res_ - 1 - j) = image.data_(i, j).GetReversedGray();

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothDirichletIJ(grid, height_field_node, smoothing_repeat);

	surface.vertex_positions_.initialize(grid.GetNumAllNodes() + grid.GetNumAllNodes());

	// vertices of upper plane
	int ix = 0;
	for (int j = 0; j < height_field_node.j_res_; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		const TV2 uv = grid.GetNode(i, j);

		surface.vertex_positions_[ix++] = TV(uv.u_, uv.v_, base_depth + height_field_node(i, j)*print_depth);
	}

	// vertices of bottom plane
	for (int j = 0; j < height_field_node.j_res_; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		const TV2 uv = grid.GetNode(i, j);

		surface.vertex_positions_[ix++] = TV(uv.u_, uv.v_, 0.0f);
	}

	LinkedArray<TV_INT> tri;

	// triangles of upper plane
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
	for (int i = 0; i < height_field_node.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, j + 1));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i + 1, j + 1));
	}

	const int offset = height_field_node.i_res_*height_field_node.j_res_;

	// triangles of bottom plane
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
	for (int i = 0; i < height_field_node.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i + 1, j + 1) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
	}

	// left wall
	{int i = 0;
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1));
	}}

	// right wall
	{int i = height_field_node.i_res_ - 1;
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i, j));
	}}

	{int j = 0;
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i + 1, j) + offset, height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, j));
	}}

	{int j = height_field_node.j_res_-1;
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i+1, j) + offset, height_field_node.get1DIndex(i, j) + offset);
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i+1, j) + offset, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i+1, j));
	}}

	tri.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();
}

const TV Cast(const TV2 v, const T z)
{
	return TV(v.x_, v.y_, z);
}

const bool IsInterfacial(const T& phi0, const T& phi1, const T& th)
{
	if (phi0 > th && phi1 > th) return false;
	else if (!(phi0 > th) && !(phi1 > th)) return false;

	return true;
}

const bool IsInterfacial(const Array2D<T>& field, const int i, const int j, const T th)
{
	if (field(i, j) > th)
	{
		if (field.getClamped(i + 1, j) <= th || field.getClamped(i - 1, j) <= th || field.getClamped(i, j + 1) <= th || field.getClamped(i, j - 1) <= th)
		{
			return true;
		}
	}
	else if (field(i, j) <= th)
	{
		if (field.getClamped(i + 1, j) > th || field.getClamped(i - 1, j) > th || field.getClamped(i, j + 1) > th || field.getClamped(i, j - 1) > th)
		{
			return true;
		}
	}

	return false;
}

void LITHOPHANE::InitializeExtrusion(const IMAGE_2D& image, StaticTriangularSurface& surface, const float front_print_thickness, const float width, const int smoothing_repeat)
{
	const float height = width / (T)image.res_x_ * (T)image.res_y_;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-6;

	GridUniform2D grid(0, 0, image.res_x_ - 1, image.res_y_ - 1, 0, 0, width, height);

	const T radius = MAX2(grid.x_max_ - grid.x_min_, grid.y_max_ - grid.y_min_)*(T)0.5 / PI;
	const TV sph_center = TV((grid.x_max_ + grid.x_min_)*(T)0.5, (grid.y_max_ + grid.y_min_)*(T)0.5, radius);
	const TV polar = TV((grid.x_max_ + grid.x_min_)*(T)0.5, (grid.y_max_ + grid.y_min_)*(T)0.5, radius*(T)1.99);

	Array2D<T> height_field_node;
	Array2D<int> u_face, v_face, vix_node;
	grid.InitializeNodeArray(height_field_node);
	grid.InitializeUFaceArray(u_face);
	grid.InitializeVFaceArray(v_face);
	grid.InitializeNodeArray(vix_node);

	height_field_node.assignAllValues(0.0f);
	v_face.assignAllValues(-1);
	u_face.assignAllValues(-1);
	vix_node.assignAllValues(-1);

	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		const T clamped_value = CLAMP(image.data_.getClamped(i, j).GetGray(), 0.0f, 1.0f);

		height_field_node(i, j) = clamped_value;
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothDirichletIJ(grid, height_field_node, smoothing_repeat);
//	smoothing.SmoothSignedDistance(grid, height_field_node, smoothing_repeat, iso_threshold);

	LinkedArray<TV> vertices;
	LinkedArray<TV_INT> triangles;

	// make bottom node vertices
	int offset = 0;
	for (int j = vix_node.j_start_; j <= vix_node.j_end_; ++j)
	for (int i = vix_node.i_start_; i <= vix_node.i_end_; ++i)
	{
		if (height_field_node(i, j) > iso_threshold)
		{
			const TV plane_pos = Cast(grid.GetNode(i, j), 0.0f);

			vertices.PushBack() = plane_pos;

			vix_node(i, j) = offset++;
		}
	}

	// make u-edge vertices
	for (int j = u_face.j_start_; j <= u_face.j_end_; ++j)
	for (int i = u_face.i_start_; i <= u_face.i_end_; ++i)
	{
		const T phi0 = height_field_node(i, j), phi1 = height_field_node(i, j + 1);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1-iso_threshold) / (ABS(phi0-iso_threshold) + ABS(phi1-iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i, j + 1)*(1.0f - alpha), 0.0f);

			vertices.PushBack() = plane_pos;

			u_face(i, j) = offset++;
		}
	}

	// make v-edge vertices
	for (int j = v_face.j_start_; j <= v_face.j_end_; ++j)
	for (int i = v_face.i_start_; i <= v_face.i_end_; ++i)
	{
		const T phi0 = height_field_node(i, j), phi1 = height_field_node(i + 1, j);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i + 1, j)*(1.0f - alpha), 0.0f);

			vertices.PushBack() = plane_pos;

			v_face(i, j) = offset++;
		}
	}

	vertices.CopyToArray(surface.vertex_positions_);	// to check triangle edge length.

	// make bottom faces
	for (int j = grid.j_start_; j <= grid.j_end_; ++j)
	for (int i = grid.i_start_; i <= grid.i_end_; ++i)
	{
		int flag(0);

		if (height_field_node(i, j) > iso_threshold) flag += 1;
		if (height_field_node(i + 1, j) > iso_threshold) flag += 2;
		if (height_field_node(i + 1, j + 1) > iso_threshold) flag += 4;
		if (height_field_node(i, j + 1) > iso_threshold) flag += 8;

		if (flag == 15)		// full
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), vix_node(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), vix_node(i, j + 1)).getReversedCW();
		}
		else if (flag == 5) // n0 n2 (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 10)// n1 n3  (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 14)	// except n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), vix_node(i + 1, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), u_face(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j), vix_node(i + 1, j)).getReversedCW();
		}
		else if (flag == 13)	// except n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), u_face(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
		}
		else if (flag == 11)	// except n2
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), u_face(i + 1, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 7)		// except n3
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 1)		// n0 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 2)		// n1 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
		}
		else if (flag == 4)		// n2 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 8)		// n3 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 3) // n0 n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 6) // n1 n2
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), v_face(i, j)).getReversedCW();
		}
		else if (flag == 12) // n2 n3
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), vix_node(i, j + 1), u_face(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 9) // n3 n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
	}

	// find boundary edges to make side later
	vertices.CopyToArray(surface.vertex_positions_);
	triangles.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();

	LinkedArray<TV2_INT> boundary_edges_linked_array;
	surface.findBoundaryEdges(boundary_edges_linked_array);

	Array1D<TV2_INT> boundary_edges;
	boundary_edges_linked_array.CopyToArray(boundary_edges);

	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		vertices.PushBack() = bottom_pos + TV(0,0,1) * front_print_thickness;
	}

	vertices.CopyToArray(surface.vertex_positions_);

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
		triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
		triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
		triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);
}

void LITHOPHANE::InitializePlanarKeyring(const IMAGE_2D& shape_image, const IMAGE_2D& height_image, StaticTriangularSurface& surface, const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float width, const int smoothing_repeat)
{
	const float height = width / (T)shape_image.res_x_ * (T)shape_image.res_y_;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-6;

	GridUniform2D grid(0, 0, shape_image.res_x_ - 1, shape_image.res_y_ - 1, 0, 0, width, height);

	const T radius = MAX2(grid.x_max_ - grid.x_min_, grid.y_max_ - grid.y_min_)*(T)0.5 / PI;
	const TV sph_center = TV((grid.x_max_ + grid.x_min_)*(T)0.5, (grid.y_max_ + grid.y_min_)*(T)0.5, radius);
	const TV polar = TV((grid.x_max_ + grid.x_min_)*(T)0.5, (grid.y_max_ + grid.y_min_)*(T)0.5, radius*(T)1.99);

	Array2D<T> shape_field_node, thickness_field_node;
	Array2D<int> u_face, v_face, vix_node;
	grid.InitializeNodeArray(shape_field_node);
	grid.InitializeNodeArray(thickness_field_node);
	grid.InitializeUFaceArray(u_face);
	grid.InitializeVFaceArray(v_face);
	grid.InitializeNodeArray(vix_node);

	shape_field_node.assignAllValues(0.0f);
	thickness_field_node.assignAllValues(0.0f);
	v_face.assignAllValues(-1);
	u_face.assignAllValues(-1);
	vix_node.assignAllValues(-1);

	for (int j = 0; j < shape_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < shape_field_node.i_res_ - 1; ++i)
	{
		shape_field_node(i, j) = CLAMP(shape_image.data_.getClamped(i, j).GetReversedGray(), 0.0f, 1.0f);
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothDirichletIJ(grid, shape_field_node, smoothing_repeat);

	for (int j = 0; j < grid.j_res_; ++j)
	for (int i = 0; i < grid.i_res_; ++i)
	{
		thickness_field_node(i, j) = CLAMP(height_image.data_.getClamped(i, j).GetReversedGray(), 0.0f, 1.0f);
	}

	smoothing.SmoothDirichletIJ(grid, thickness_field_node, smoothing_repeat);

	LinkedArray<TV> vertices;
	LinkedArray<TV_INT> triangles;

	// make bottom node vertices
	int offset = 0;
	for (int j = vix_node.j_start_; j <= vix_node.j_end_; ++j)
	for (int i = vix_node.i_start_; i <= vix_node.i_end_; ++i)
	{
		if (shape_field_node(i, j) > iso_threshold)
		{
			const TV plane_pos = Cast(grid.GetNode(i, j), -thickness_field_node(i, j));		// back node vertices are initialized as -thickness field and then modified when generating front vertices

			vertices.PushBack() = plane_pos;

			vix_node(i, j) = offset++;
		}
	}

	// make u-edge vertices
	for (int j = u_face.j_start_; j <= u_face.j_end_; ++j)
	for (int i = u_face.i_start_; i <= u_face.i_end_; ++i)
	{
		const T phi0 = shape_field_node(i, j), phi1 = shape_field_node(i, j + 1);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const T  litho_depth = thickness_field_node(i, j)*(1.0f - alpha) + thickness_field_node(i, j + 1)*alpha;
			
			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i, j + 1)*(1.0f - alpha), -litho_depth);

			vertices.PushBack() = plane_pos;

			u_face(i, j) = offset++;
		}
	}

	// make v-edge vertices
	for (int j = v_face.j_start_; j <= v_face.j_end_; ++j)
	for (int i = v_face.i_start_; i <= v_face.i_end_; ++i)
	{
		const T phi0 = shape_field_node(i, j), phi1 = shape_field_node(i + 1, j);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			// thickness_field_node(i,j) = 1.0f for inter-facial vertices
			const T  litho_depth = thickness_field_node(i, j)*(1.0f - alpha) + thickness_field_node(i + 1, j)*alpha;
			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i + 1, j)*(1.0f - alpha), -litho_depth);

			vertices.PushBack() = plane_pos;

			v_face(i, j) = offset++;
		}
	}

	vertices.CopyToArray(surface.vertex_positions_);	// to check triangle edge length.

	// make bottom faces
	for (int j = grid.j_start_; j <= grid.j_end_; ++j)
	for (int i = grid.i_start_; i <= grid.i_end_; ++i)
	{
		int flag(0);

		if (shape_field_node(i, j) > iso_threshold) flag += 1;
		if (shape_field_node(i + 1, j) > iso_threshold) flag += 2;
		if (shape_field_node(i + 1, j + 1) > iso_threshold) flag += 4;
		if (shape_field_node(i, j + 1) > iso_threshold) flag += 8;

		if (flag == 15)		// full
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), vix_node(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), vix_node(i, j + 1)).getReversedCW();
		}
		else if (flag == 5) // n0 n2 (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 10)// n1 n3  (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 14)	// except n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), vix_node(i + 1, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), u_face(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j), vix_node(i + 1, j)).getReversedCW();
		}
		else if (flag == 13)	// except n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), u_face(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
		}
		else if (flag == 11)	// except n2
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), u_face(i + 1, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 7)		// except n3
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 1)		// n0 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 2)		// n1 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
		}
		else if (flag == 4)		// n2 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 8)		// n3 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 3) // n0 n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 6) // n1 n2
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), v_face(i, j)).getReversedCW();
		}
		else if (flag == 12) // n2 n3
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), vix_node(i, j + 1), u_face(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 9) // n3 n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
	}

	// find boundary edges to make side later
	vertices.CopyToArray(surface.vertex_positions_);
	triangles.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();

	LinkedArray<TV2_INT> boundary_edges_linked_array;
	surface.findBoundaryEdges(boundary_edges_linked_array);

	Array1D<TV2_INT> boundary_edges;
	boundary_edges_linked_array.CopyToArray(boundary_edges);

	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		vertices.PushBack() = TV(bottom_pos.x_, bottom_pos.y_, -bottom_pos.z_*front_print_thickness + base_thickness);
	}

	vertices.CopyToArray(surface.vertex_positions_);

	for (int vix = 0; vix < offset; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		surface.vertex_positions_[vix].z_ = bottom_pos.z_*back_print_thickness;
	}

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
		triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
		triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
		triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);
}

void LITHOPHANE::makeFreePlane(const GridUniform2D& grid, const Array2D<T>& shape_field_node, const Array2D<T>& height_field_node, const T iso_threshold, const T min_edge_length, StaticTriangularSurface& surface)
{
//	Array2D<T> shape_field_node, thickness_field_node;
	Array2D<int> u_face, v_face, vix_node;
//	grid.InitializeNodeArray(shape_field_node);
//	grid.InitializeNodeArray(thickness_field_node);
	grid.InitializeUFaceArray(u_face);
	grid.InitializeVFaceArray(v_face);
	grid.InitializeNodeArray(vix_node);

	v_face.assignAllValues(-1);
	u_face.assignAllValues(-1);
	vix_node.assignAllValues(-1);

	LinkedArray<TV> vertices;
	LinkedArray<TV_INT> triangles;

	// make bottom node vertices
	int offset = 0;
	for (int j = vix_node.j_start_; j <= vix_node.j_end_; ++j)
		for (int i = vix_node.i_start_; i <= vix_node.i_end_; ++i)
		{
			if (shape_field_node(i, j) > iso_threshold)
			{
				const TV plane_pos = Cast(grid.GetNode(i, j), (T)0);		// back node vertices are initialized as -thickness field and then modified when generating front vertices

				vertices.PushBack() = plane_pos;

				vix_node(i, j) = offset++;
			}
		}

	// make u-edge vertices
	for (int j = u_face.j_start_; j <= u_face.j_end_; ++j)
		for (int i = u_face.i_start_; i <= u_face.i_end_; ++i)
		{
			const T phi0 = shape_field_node(i, j), phi1 = shape_field_node(i, j + 1);

			if (IsInterfacial(phi0, phi1, iso_threshold) == true)
			{
				const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				if (alpha < 0.0f || alpha > 1.0f)
				{
					const float temp = 1.0f;

					std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
				}

				const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i, j + 1)*(1.0f - alpha), (T)0);

				vertices.PushBack() = plane_pos;

				u_face(i, j) = offset++;
			}
		}

	// make v-edge vertices
	for (int j = v_face.j_start_; j <= v_face.j_end_; ++j)
		for (int i = v_face.i_start_; i <= v_face.i_end_; ++i)
		{
			const T phi0 = shape_field_node(i, j), phi1 = shape_field_node(i + 1, j);

			if (IsInterfacial(phi0, phi1, iso_threshold) == true)
			{
				const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				if (alpha < 0.0f || alpha > 1.0f)
				{
					const float temp = 1.0f;

					std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
				}

				// thickness_field_node(i,j) = 1.0f for inter-facial vertices
				const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i + 1, j)*(1.0f - alpha), (T)0);

				vertices.PushBack() = plane_pos;

				v_face(i, j) = offset++;
			}
		}

	vertices.CopyToArray(surface.vertex_positions_);	// to check triangle edge length.

	// make bottom faces
	for (int j = grid.j_start_; j <= grid.j_end_; ++j)
		for (int i = grid.i_start_; i <= grid.i_end_; ++i)
		{
			int flag(0);

			if (shape_field_node(i, j) > iso_threshold) flag += 1;
			if (shape_field_node(i + 1, j) > iso_threshold) flag += 2;
			if (shape_field_node(i + 1, j + 1) > iso_threshold) flag += 4;
			if (shape_field_node(i, j + 1) > iso_threshold) flag += 8;

			if (flag == 15)		// full
			{
				if (IsInterfacial(height_field_node(i, j), height_field_node(i + 1, j + 1), iso_threshold))
				{
					triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), vix_node(i, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), vix_node(i, j + 1)).getReversedCW();
				}
				else
				{
					triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j + 1), vix_node(i, j + 1)).getReversedCW();
				}
			}
			else if (flag == 5) // n0 n2 (TODO:check ambiguity)
			{
				triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
			}
			else if (flag == 10)// n1 n3  (TODO:check ambiguity)
			{
				triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
			}
			else if (flag == 14)	// except n0
			{
				if (IsInterfacial(height_field_node(i, j), height_field_node(i + 1, j + 1), iso_threshold))
				{
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), vix_node(i + 1, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
				}
				else
				{
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), vix_node(i + 1, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), u_face(i, j), v_face(i, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j), vix_node(i + 1, j)).getReversedCW();
				}
			}
			else if (flag == 13)	// except n1
			{
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), u_face(i + 1, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
			}
			else if (flag == 11)	// except n2
			{
				if (IsInterfacial(height_field_node(i, j), height_field_node(i + 1, j + 1), iso_threshold))	//TODO: check and finish
				{
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), vix_node(i + 1, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i + 1, j), u_face(i + 1, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), v_face(i, j + 1)).getReversedCW();
				}
				else
				{
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j), u_face(i + 1, j), v_face(i, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i + 1, j)).getReversedCW();
				}
			}
			else if (flag == 7)		// except n3
			{
				triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), u_face(i, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
			}
			else if (flag == 1)		// n0 only
			{
				triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
			}
			else if (flag == 2)		// n1 only
			{
				triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
			}	
			else if (flag == 4)		// n2 only
			{
				triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
			}
			else if (flag == 8)		// n3 only
			{
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
			}
			else if (flag == 3) // n0 n1
			{
				if (IsInterfacial(height_field_node(i, j), height_field_node(i + 1, j + 1), iso_threshold))
				{
					triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), u_face(i, j)).getReversedCW();
				}
				else
				{
					triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i + 1, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j), u_face(i + 1, j), u_face(i, j)).getReversedCW();
				}				
			}
			else if (flag == 6) // n1 n2
			{
				if (IsInterfacial(height_field_node(i, j), height_field_node(i + 1, j + 1), iso_threshold))
				{
					triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), v_face(i, j)).getReversedCW();
				}
				else
				{
					triangles.PushBack() = TV_INT(v_face(i, j), vix_node(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
					triangles.PushBack() = TV_INT(v_face(i, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
				}
			}
			else if (flag == 12) // n2 n3
			{
				if (IsInterfacial(height_field_node(i, j), height_field_node(i + 1, j + 1), iso_threshold))
				{
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face(i + 1, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
				}
				else
				{
					triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), vix_node(i, j + 1), u_face(i + 1, j)).getReversedCW();
					triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face(i + 1, j)).getReversedCW();
				}
			}
			else if (flag == 9) // n3 n0
			{
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
				triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), v_face(i, j + 1)).getReversedCW();
			}
		}

	// find boundary edges to make side later
	vertices.CopyToArray(surface.vertex_positions_);
	triangles.CopyToArray(surface.triangles_);
}

void LITHOPHANE::makeFreeformStamp(const Array2D<T>& shape_field_node, const Array2D<T>& thickness_field_node, StaticTriangularSurface& surface, const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float width, const int smoothing_repeat)
{
	const float height = width / (T)shape_field_node.i_res_ * (T)shape_field_node.j_res_;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-6;

	GridUniform2D node_grid(0, 0, shape_field_node.i_res_, shape_field_node.j_res_, 0, 0, width, height);
	GridUniform2D cell_grid = node_grid.getCellGrid();
	GridUniform2D cell_grid_ghost = cell_grid.getEnlarged(3);

	Array2D<int> bc_map_bottom, bc_map_middle;
	cell_grid_ghost.InitializeCellArray(bc_map_bottom);
	cell_grid_ghost.InitializeCellArray(bc_map_middle);
	bc_map_bottom.assignAllValues(-1);
	bc_map_middle.assignAllValues(-1);

	StaticTriangularSurface bottom_surface, middle_surface, upper_surface;

	makeFreePlane(cell_grid, shape_field_node, thickness_field_node, iso_threshold, min_edge_length, bottom_surface);

	makeFreePlane(cell_grid, shape_field_node, thickness_field_node, iso_threshold, min_edge_length, middle_surface);

	// find boundaries
	Array1D<TV2_INT> bottom_surface_boundary_edges, middle_surface_boundary_edges;
	{bottom_surface.findAdjacentTrianglesOfVertices();
	bottom_surface.findEdgeTrianglesOfTriangles();
	LinkedArray<TV2_INT> boundary_edges_linked_array;
	bottom_surface.findBoundaryEdges(boundary_edges_linked_array);
	boundary_edges_linked_array.CopyToArray(bottom_surface_boundary_edges); }

	{middle_surface.findAdjacentTrianglesOfVertices();
	middle_surface.findEdgeTrianglesOfTriangles();
	LinkedArray<TV2_INT> boundary_edges_linked_array;
	middle_surface.findBoundaryEdges(boundary_edges_linked_array);
	boundary_edges_linked_array.CopyToArray(middle_surface_boundary_edges); }

	// make edge map

	for (int i = 0; i < bottom_surface_boundary_edges.num_elements_; i++)
	{
		const TV edge_center = (bottom_surface.vertex_positions_[bottom_surface_boundary_edges[i].x_] + bottom_surface.vertex_positions_[bottom_surface_boundary_edges[i].y_]) * (T)0.5;

		const TV2_INT ix = cell_grid_ghost.Cell(TV2(edge_center.x_, edge_center.y_));

		bc_map_bottom(ix) = i;
	}

	for (int i = 0; i < middle_surface_boundary_edges.num_elements_; i++)
	{
		const TV edge_center = (middle_surface.vertex_positions_[middle_surface_boundary_edges[i].x_] + middle_surface.vertex_positions_[middle_surface_boundary_edges[i].y_]) * (T)0.5;

		const TV2_INT ix = cell_grid_ghost.Cell(TV2(edge_center.x_, edge_center.y_));

		bc_map_middle(ix) = i;
	}

	// fix middle surface

	for (int i = 0; i < middle_surface.vertex_positions_.num_elements_; i++)
	{
		const TV2 pos_c = TV2(middle_surface.vertex_positions_[i].x_, middle_surface.vertex_positions_[i].y_);
		const TV2 pos_l = pos_c + TV2(-cell_grid.dx_, (T)0);
		const TV2 pos_r = pos_c + TV2(cell_grid.dx_, (T)0);
		const TV2 pos_u = pos_c + TV2((T)0, cell_grid.dx_);
		const TV2 pos_d = pos_c + TV2((T)0, -cell_grid.dx_);

		const T phi_c = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_c), (T)0, (T)1);
		const T phi_l = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_l), (T)0, (T)1);
		const T phi_r = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_r), (T)0, (T)1);
		const T phi_u = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_u), (T)0, (T)1);
		const T phi_d = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_d), (T)0, (T)1);
		
		const T rho_c = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_c), (T)0, (T)1);
		const T rho_l = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_l), (T)0, (T)1);
		const T rho_r = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_r), (T)0, (T)1);
		const T rho_u = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_u), (T)0, (T)1);
		const T rho_d = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_d), (T)0, (T)1);

		TV2 pos;
		int num = 0;

		if (rho_c > iso_threshold && rho_l > iso_threshold && rho_r > iso_threshold && rho_u > iso_threshold && rho_d > iso_threshold)
		{

			if (IsInterfacial(phi_c, phi_l, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_l - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_l - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_l * ((T)1 - alpha);

				num++;
			}

			if (IsInterfacial(phi_c, phi_r, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_r - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_r - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_r * ((T)1 - alpha);

				num++;
			}

			if (IsInterfacial(phi_c, phi_u, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_u - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_u - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_u * ((T)1 - alpha);

				num++;
			}

			if (IsInterfacial(phi_c, phi_d, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_d - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_d - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_d * ((T)1 - alpha);

				num++;
			}

			if (num == 0) pos = pos_c;
			else pos /= (T)num;
		}
		else
			pos = pos_c;

		middle_surface.vertex_positions_[i].x_ = pos.x_;
		middle_surface.vertex_positions_[i].y_ = pos.y_;

		middle_surface.vertex_positions_[i].z_ = base_thickness + front_print_thickness * (phi_c > iso_threshold ? 1.0f : 0.0f);

	}

	for (int i = 0; i < middle_surface.triangles_.num_elements_; i++)
	{
		middle_surface.triangles_[i] = middle_surface.triangles_[i].getReversedCW() + TV_INT(bottom_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_);
	}

	// make side face
	Array1D<TV_INT> bottom_side;
	{LinkedArray<TV_INT> bottom_side_buffer;
	for (int j = cell_grid.j_start_; j <= cell_grid.j_end_; j++)
		for (int i = cell_grid.i_start_; i <= cell_grid.i_end_; i++)
		{
			if (bc_map_bottom(i, j) == -1) continue;

			const TV2_INT bc_bottom = bottom_surface_boundary_edges[bc_map_bottom(i, j)];
			const TV2_INT bc_mid = middle_surface_boundary_edges[bc_map_middle(i, j)] + TV2_INT(bottom_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_);	// offset fix

			bottom_side_buffer.PushBack() = TV_INT(bc_bottom.x_, bc_bottom.y_, bc_mid.y_).getReversedCW();
			bottom_side_buffer.PushBack() = TV_INT(bc_mid.x_, bc_mid.y_, bc_bottom.x_);
		}
		bottom_side_buffer.CopyToArray(bottom_side);}
	
	surface.vertex_positions_.append(bottom_surface.vertex_positions_);
	surface.vertex_positions_.append(middle_surface.vertex_positions_);

	surface.triangles_.append(bottom_surface.triangles_);
	surface.triangles_.append(middle_surface.triangles_);
	surface.triangles_.append(bottom_side);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();
	surface.determineFaceAveragedVertexNormals();

/*	
	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		vertices.PushBack() = TV(bottom_pos.x_, bottom_pos.y_, -bottom_pos.z_*front_print_thickness + base_thickness);
	}

	vertices.CopyToArray(surface.vertex_positions_);

	for (int vix = 0; vix < offset; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		surface.vertex_positions_[vix].z_ = bottom_pos.z_*back_print_thickness;
	}

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
		triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
		triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
		triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);
*/
}


void LITHOPHANE::makeFreeformStampDoubleSided(const Array2D<T>& shape_field_node, const Array2D<T>& thickness_field_node, StaticTriangularSurface& surface, const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float width, const int smoothing_repeat)
{
	const float height = width / (T)shape_field_node.i_res_ * (T)shape_field_node.j_res_;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-6;

	GridUniform2D node_grid(0, 0, shape_field_node.i_res_, shape_field_node.j_res_, 0, 0, width, height);
	GridUniform2D cell_grid = node_grid.getCellGrid();
	GridUniform2D cell_grid_ghost = cell_grid.getEnlarged(3);

	Array2D<int> bc_map_bottom, bc_map_middle;
	cell_grid_ghost.InitializeCellArray(bc_map_bottom);
	cell_grid_ghost.InitializeCellArray(bc_map_middle);
	bc_map_bottom.assignAllValues(-1);
	bc_map_middle.assignAllValues(-1);

	StaticTriangularSurface bottom_surface, upper_surface;

	makeFreePlane(cell_grid, shape_field_node, thickness_field_node, iso_threshold, min_edge_length, bottom_surface);

	makeFreePlane(cell_grid, shape_field_node, thickness_field_node, iso_threshold, min_edge_length, upper_surface);

	// find boundaries
	Array1D<TV2_INT> bottom_surface_boundary_edges, middle_surface_boundary_edges;
	{bottom_surface.findAdjacentTrianglesOfVertices();
	bottom_surface.findEdgeTrianglesOfTriangles();
	LinkedArray<TV2_INT> boundary_edges_linked_array;
	bottom_surface.findBoundaryEdges(boundary_edges_linked_array);
	boundary_edges_linked_array.CopyToArray(bottom_surface_boundary_edges); }

	{upper_surface.findAdjacentTrianglesOfVertices();
	upper_surface.findEdgeTrianglesOfTriangles();
	LinkedArray<TV2_INT> boundary_edges_linked_array;
	upper_surface.findBoundaryEdges(boundary_edges_linked_array);
	boundary_edges_linked_array.CopyToArray(middle_surface_boundary_edges); }

	// make edge map
	for (int i = 0; i < bottom_surface_boundary_edges.num_elements_; i++)
	{
		const TV edge_center = (bottom_surface.vertex_positions_[bottom_surface_boundary_edges[i].x_] + bottom_surface.vertex_positions_[bottom_surface_boundary_edges[i].y_]) * (T)0.5;

		const TV2_INT ix = cell_grid_ghost.Cell(TV2(edge_center.x_, edge_center.y_));

		bc_map_bottom(ix) = i;
	}

	for (int i = 0; i < middle_surface_boundary_edges.num_elements_; i++)
	{
		const TV edge_center = (upper_surface.vertex_positions_[middle_surface_boundary_edges[i].x_] + upper_surface.vertex_positions_[middle_surface_boundary_edges[i].y_]) * (T)0.5;

		const TV2_INT ix = cell_grid_ghost.Cell(TV2(edge_center.x_, edge_center.y_));

		bc_map_middle(ix) = i;
	}

	// fix middle surface to make sharp edges
	for (int i = 0; i < upper_surface.vertex_positions_.num_elements_; i++)
	{
		const TV2 pos_c = TV2(upper_surface.vertex_positions_[i].x_, upper_surface.vertex_positions_[i].y_);
		const TV2 pos_l = pos_c + TV2(-cell_grid.dx_, (T)0);
		const TV2 pos_r = pos_c + TV2(cell_grid.dx_, (T)0);
		const TV2 pos_u = pos_c + TV2((T)0, cell_grid.dx_);
		const TV2 pos_d = pos_c + TV2((T)0, -cell_grid.dx_);

		const T phi_c = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_c), (T)0, (T)1);
		const T phi_l = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_l), (T)0, (T)1);
		const T phi_r = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_r), (T)0, (T)1);
		const T phi_u = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_u), (T)0, (T)1);
		const T phi_d = CLAMP(node_grid.GetClampedLinearInterpolationCell(thickness_field_node, pos_d), (T)0, (T)1);

		const T rho_c = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_c), (T)0, (T)1);
		const T rho_l = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_l), (T)0, (T)1);
		const T rho_r = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_r), (T)0, (T)1);
		const T rho_u = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_u), (T)0, (T)1);
		const T rho_d = CLAMP(node_grid.GetClampedLinearInterpolationCell(shape_field_node, pos_d), (T)0, (T)1);

		TV2 pos;
		int num = 0;

		if (rho_c > iso_threshold && rho_l > iso_threshold && rho_r > iso_threshold && rho_u > iso_threshold && rho_d > iso_threshold)
		{

			if (IsInterfacial(phi_c, phi_l, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_l - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_l - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_l * ((T)1 - alpha);

				num++;
			}

			if (IsInterfacial(phi_c, phi_r, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_r - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_r - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_r * ((T)1 - alpha);

				num++;
			}

			if (IsInterfacial(phi_c, phi_u, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_u - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_u - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_u * ((T)1 - alpha);

				num++;
			}

			if (IsInterfacial(phi_c, phi_d, iso_threshold))
			{
				const T alpha = CLAMP(ABS(phi_d - iso_threshold) / (ABS(phi_c - iso_threshold) + ABS(phi_d - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

				pos += pos_c *alpha + pos_d * ((T)1 - alpha);

				num++;
			}

			if (num == 0) pos = pos_c;
			else pos /= (T)num;
		}
		else
			pos = pos_c;

		upper_surface.vertex_positions_[i].x_ = pos.x_;
		upper_surface.vertex_positions_[i].y_ = pos.y_;

		upper_surface.vertex_positions_[i].z_ = base_thickness + front_print_thickness * (phi_c > iso_threshold ? 1.0f : 0.0f);
	}

	// prepare for appending upper surface triangles to bottom surface triangles
	for (int i = 0; i < upper_surface.triangles_.num_elements_; i++)
	{
		upper_surface.triangles_[i] = upper_surface.triangles_[i].getReversedCW() + TV_INT(bottom_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_);
	}

	// copy boundary vertices of upper surface to generate middle boundary vertices
	Array1D<bool> is_boundary_vertex(upper_surface.vertex_positions_.num_elements_);
	is_boundary_vertex.assignAllValues(false);

	for (int i = 0; i < middle_surface_boundary_edges.num_elements_; i++)
	{
		is_boundary_vertex[middle_surface_boundary_edges[i].x_] = true;
		is_boundary_vertex[middle_surface_boundary_edges[i].y_] = true;
	}

	int count_boundary_vertices = 0;
	for (int i = 0; i < is_boundary_vertex.num_elements_; i++)
	{
		if (is_boundary_vertex[i] == true) count_boundary_vertices++;
	}

	Array1D<int> boundary_vertex_index_map(upper_surface.vertex_positions_.num_elements_);
	boundary_vertex_index_map.assignAllValues(-1);
	Array1D<TV> boundary_vertex_list(count_boundary_vertices);
	int flag = 0;
	for (int i = 0; i < is_boundary_vertex.num_elements_; i++)
	{
		if (is_boundary_vertex[i] == true)
		{
			boundary_vertex_index_map[i] = flag;
			boundary_vertex_list[flag] = upper_surface.vertex_positions_[i];
			boundary_vertex_list[flag].z_ = base_thickness * (T)0.5; // TODO: test only
			flag++;
		}
	}

	// make side face
	Array1D<TV_INT> bottom_side;
	{LinkedArray<TV_INT> bottom_side_buffer;
	for (int j = cell_grid.j_start_; j <= cell_grid.j_end_; j++)
		for (int i = cell_grid.i_start_; i <= cell_grid.i_end_; i++)
		{
			if (bc_map_bottom(i, j) == -1) continue;

			const TV2_INT bc_bottom = bottom_surface_boundary_edges[bc_map_bottom(i, j)];
			const TV2_INT bc_upper = middle_surface_boundary_edges[bc_map_middle(i, j)] + TV2_INT(bottom_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_);	// offset fix
			const TV2_INT bc_mid = TV2_INT(boundary_vertex_index_map[middle_surface_boundary_edges[bc_map_middle(i, j)].x_], boundary_vertex_index_map[middle_surface_boundary_edges[bc_map_middle(i, j)].y_])
				+ TV2_INT(bottom_surface.vertex_positions_.num_elements_ + upper_surface.vertex_positions_.num_elements_, bottom_surface.vertex_positions_.num_elements_ + upper_surface.vertex_positions_.num_elements_);

// 			bottom_side_buffer.PushBack() = TV_INT(bc_bottom.x_, bc_bottom.y_, bc_upper.y_).GetReversedCW();
// 			bottom_side_buffer.PushBack() = TV_INT(bc_upper.x_, bc_upper.y_, bc_bottom.x_);

			bottom_side_buffer.PushBack() = TV_INT(bc_bottom.x_, bc_bottom.y_, bc_mid.y_).getReversedCW();
			bottom_side_buffer.PushBack() = TV_INT(bc_mid.x_, bc_mid.y_, bc_bottom.x_);

			bottom_side_buffer.PushBack() = TV_INT(bc_mid.x_, bc_mid.y_, bc_upper.y_).getReversedCW();
			bottom_side_buffer.PushBack() = TV_INT(bc_upper.x_, bc_upper.y_, bc_mid.x_);


		}
	bottom_side_buffer.CopyToArray(bottom_side); }

	surface.vertex_positions_.append(bottom_surface.vertex_positions_);
	surface.vertex_positions_.append(upper_surface.vertex_positions_);
	surface.vertex_positions_.append(boundary_vertex_list);

	surface.triangles_.append(bottom_surface.triangles_);
	surface.triangles_.append(upper_surface.triangles_);
	surface.triangles_.append(bottom_side);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();
	surface.determineFaceAveragedVertexNormals();

	/*
	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
	const TV bottom_pos = surface.vertex_positions_[vix];

	vertices.PushBack() = TV(bottom_pos.x_, bottom_pos.y_, -bottom_pos.z_*front_print_thickness + base_thickness);
	}

	vertices.CopyToArray(surface.vertex_positions_);

	for (int vix = 0; vix < offset; ++vix)
	{
	const TV bottom_pos = surface.vertex_positions_[vix];

	surface.vertex_positions_[vix].z_ = bottom_pos.z_*back_print_thickness;
	}

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
	triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
	triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
	triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);
	*/
}

void LITHOPHANE::InitializeFreeFormPlanarLithophane(const IMAGE_2D& image, StaticTriangularSurface& surface, const float base_thickness, const float front_print_thickness, const float width, const int smoothing_repeat)
{
	const float height = width / (T)image.res_x_ * (T)image.res_y_;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-6;

	GridUniform2D grid(0, 0, image.res_x_ - 1, image.res_y_ - 1, 0, 0, width, height);

	Array2D<T> hole_field_node, litho_field_node;
	Array2D<int> u_face, v_face, vix_node;
	grid.InitializeNodeArray(hole_field_node);
	grid.InitializeNodeArray(litho_field_node);
	grid.InitializeUFaceArray(u_face);
	grid.InitializeVFaceArray(v_face);
	grid.InitializeNodeArray(vix_node);

	hole_field_node.assignAllValues(0.0f);
	litho_field_node.assignAllValues(0.0f);
	v_face.assignAllValues(-1);
	u_face.assignAllValues(-1);
	vix_node.assignAllValues(-1);

	for (int j = 0; j < hole_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < hole_field_node.i_res_ - 1; ++i)
	{
		const T clamped_value = !image.data_(i, j).IsWhite() ? 1.0f : CLAMP(image.data_.getClamped(i, j).GetReversedGray(), 0.0f, 1.0f);		// tag blue means no hole

//		std::cout << (int)image.data_(i, j).r_ << " " << (int)image.data_(i, j).g_ << " " << (int)image.data_(i, j).b_ << std::endl;

		hole_field_node(i, j) = clamped_value;
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothDirichletIJ(grid, hole_field_node, smoothing_repeat);

	for (int j = 0; j < grid.j_res_; ++j)
	for (int i = 0; i < grid.i_res_; ++i)
		litho_field_node(i, j) = CLAMP(image.data_.getClamped(i, j).GetReversedGray(), 0.0f, 1.0f);

	smoothing.SmoothDirichletIJ(grid, litho_field_node, smoothing_repeat);

	LinkedArray<TV> vertices;
	LinkedArray<TV_INT> triangles;

	// make upper node vertices
	int offset = 0;
	for (int j = vix_node.j_start_; j <= vix_node.j_end_; ++j)
	for (int i = vix_node.i_start_; i <= vix_node.i_end_; ++i)
	{
		if (hole_field_node(i, j) > iso_threshold)
		{
			const T  total_thickness = litho_field_node(i, j)*front_print_thickness + base_thickness;
			const TV plane_pos = Cast(grid.GetNode(i, j), total_thickness);

			vertices.PushBack() = plane_pos;

			vix_node(i, j) = offset++;
		}
	}

	// make upper u-edge vertices
	for (int j = u_face.j_start_; j <= u_face.j_end_; ++j)
	for (int i = u_face.i_start_; i <= u_face.i_end_; ++i)
	{
		const T phi0 = hole_field_node(i, j), phi1 = hole_field_node(i, j + 1);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const T  litho_depth = (litho_field_node(i, j)*(1.0f - alpha) + litho_field_node(i, j + 1)*alpha)*front_print_thickness + base_thickness;
			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i, j + 1)*(1.0f - alpha), litho_depth);

			vertices.PushBack() = plane_pos;

			u_face(i, j) = offset++;
		}
	}

	// make v-edge vertices
	for (int j = v_face.j_start_; j <= v_face.j_end_; ++j)
	for (int i = v_face.i_start_; i <= v_face.i_end_; ++i)
	{
		const T phi0 = hole_field_node(i, j), phi1 = hole_field_node(i + 1, j);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const T  litho_depth = (litho_field_node(i, j)*(1.0f - alpha) + litho_field_node(i + 1, j)*alpha)*front_print_thickness
				+ base_thickness;
			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i + 1, j)*(1.0f - alpha), litho_depth);

			vertices.PushBack() = plane_pos;

			v_face(i, j) = offset++;
		}
	}

	vertices.CopyToArray(surface.vertex_positions_);	// to check triangle edge length.

	// make upper faces
	for (int j = grid.j_start_; j <= grid.j_end_; ++j)
	for (int i = grid.i_start_; i <= grid.i_end_; ++i)
	{
		int flag(0);

		if (hole_field_node(i, j) > iso_threshold) flag += 1;
		if (hole_field_node(i + 1, j) > iso_threshold) flag += 2;
		if (hole_field_node(i + 1, j + 1) > iso_threshold) flag += 4;
		if (hole_field_node(i, j + 1) > iso_threshold) flag += 8;

		if (flag == 15)		// full
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), vix_node(i, j + 1));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), vix_node(i, j + 1));
		}
		else if (flag == 5) // n0 n2 (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j));
		}
		else if (flag == 10)// n1 n3  (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1));
		}
		else if (flag == 14)	// except n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), vix_node(i + 1, j + 1));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), u_face(i, j), v_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j), vix_node(i + 1, j));
		}
		else if (flag == 13)	// except n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), u_face(i + 1, j));
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), vix_node(i + 1, j + 1));
		}
		else if (flag == 11)	// except n2
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j + 1));
			triangles.PushBack() = TV_INT(vix_node(i, j), u_face(i + 1, j), v_face(i, j + 1));
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i + 1, j));
		}
		else if (flag == 7)		// except n3
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), u_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1));
		}
		else if (flag == 1)		// n0 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j));
		}
		else if (flag == 2)		// n1 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j));
		}
		else if (flag == 4)		// n2 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j));
		}
		else if (flag == 8)		// n3 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1));
		}
		else if (flag == 3) // n0 n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), u_face(i, j));
		}
		else if (flag == 6) // n1 n2
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1));
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), v_face(i, j));
		}
		else if (flag == 12) // n2 n3
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), vix_node(i, j + 1), u_face(i + 1, j));
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face(i + 1, j));
		}
		else if (flag == 9) // n3 n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j));
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), v_face(i, j + 1));
		}
	}

	// find boundary edges to make side later
	vertices.CopyToArray(surface.vertex_positions_);
	triangles.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();

	LinkedArray<TV2_INT> boundary_edges_linked_array;
	surface.findBoundaryEdges(boundary_edges_linked_array);

	Array1D<TV2_INT> boundary_edges;
	boundary_edges_linked_array.CopyToArray(boundary_edges);

	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		//		vertices.PushBack() = bottom_pos + TV(0, 0, 1) * front_print_thickness;
		vertices.PushBack() = TV(bottom_pos.x_, bottom_pos.y_, 0.0f);
	}

	vertices.CopyToArray(surface.vertex_positions_);

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
		triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
		triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
		triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);
}

void LITHOPHANE::InitializeExtrusionAndSphericalMapping(const IMAGE_2D& image, StaticTriangularSurface& surface, const float front_print_thickness, const float width, const int smoothing_repeat)
{
	const float height = width / (T)image.res_x_ * (T)image.res_y_;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-4;

	GridUniform2D grid(0, 0, image.res_x_ - 1, image.res_y_ - 1, 0, 0, width, height);

	const T radius = MAX2(grid.x_max_ - grid.x_min_, grid.y_max_ - grid.y_min_)*(T)0.5/PI;
	const TV sph_center = TV((grid.x_max_ + grid.x_min_)*(T)0.5, (grid.y_max_ + grid.y_min_)*(T)0.5, radius);
	const TV polar = TV((grid.x_max_ + grid.x_min_)*(T)0.5, (grid.y_max_ + grid.y_min_)*(T)0.5, radius*(T)1.99);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);
	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		const T clamped_value = CLAMP(image.data_.getClamped(i, j).GetReversedGray(), 0.0f, 1.0f);

		height_field_node(i, j) = clamped_value;
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothDirichletIJ(grid, height_field_node, smoothing_repeat);

	Array2D<int> u_face, v_face, vix_node;

	grid.InitializeUFaceArray(u_face);
	grid.InitializeVFaceArray(v_face);
	grid.InitializeNodeArray(vix_node);

	v_face.assignAllValues(-1);
	u_face.assignAllValues(-1);
	vix_node.assignAllValues(-1);

	LinkedArray<TV> vertices;
	LinkedArray<TV_INT> triangles;

	// make bottom node vertices
	int offset = 0;
	for (int j = vix_node.j_start_; j <= vix_node.j_end_; ++j)
	for (int i = vix_node.i_start_; i <= vix_node.i_end_; ++i)
	{
		if (height_field_node(i, j) > iso_threshold)
		{
			const TV plane_pos = Cast(grid.GetNode(i, j), 0.0f);
			const TV proj_dir = (plane_pos - polar).normalizedDouble();

			RAY ray(polar, proj_dir, radius*100.0f);

			const T t = ray.GetSphereIntersection(sph_center, radius);

			if (t < 0){std::cout << "Wrong projection" << std::endl; exit(1);}

			vertices.PushBack() = polar + proj_dir * t;

			vix_node(i, j) = offset++;
		}
	}

	// make u-edge vertices
	for (int j = u_face.j_start_; j <= u_face.j_end_; ++j)
	for (int i = u_face.i_start_; i <= u_face.i_end_; ++i)
	{
		const T phi0 = height_field_node(i, j), phi1 = height_field_node(i, j + 1);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i, j + 1)*(1.0f - alpha), 0.0f);
			const TV proj_dir = (plane_pos - polar).normalizedDouble();

			RAY ray(polar, proj_dir, radius*100.0f);

			const T t = ray.GetSphereIntersection(sph_center, radius);

			if (t < 0){ std::cout << "Wrong projection" << std::endl; exit(1); }

			vertices.PushBack() = polar + proj_dir * t;

			u_face(i, j) = offset++;
		}
	}

	// make v-edge vertices
	for (int j = v_face.j_start_; j <= v_face.j_end_; ++j)
	for (int i = v_face.i_start_; i <= v_face.i_end_; ++i)
	{
		const T phi0 = height_field_node(i, j), phi1 = height_field_node(i + 1, j);

		if (IsInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const TV plane_pos = Cast(grid.GetNode(i, j)*alpha + grid.GetNode(i + 1, j)*(1.0f - alpha), 0.0f);
			const TV proj_dir = (plane_pos - polar).normalizedDouble();

			RAY ray(polar, proj_dir, radius*100.0f);

			const T t = ray.GetSphereIntersection(sph_center, radius);

			if (t < 0){ std::cout << "Wrong projection" << std::endl; exit(1); }

			vertices.PushBack() = polar + proj_dir * t;

			v_face(i, j) = offset++;
		}
	}

	vertices.CopyToArray(surface.vertex_positions_);	// to check triangle edge length.

	// make bottom faces
	for (int j = grid.j_start_; j <= grid.j_end_; ++j)
	for (int i = grid.i_start_; i <= grid.i_end_; ++i)
	{
		int flag(0);

		if (height_field_node(i, j) > iso_threshold) flag += 1;
		if (height_field_node(i + 1, j) > iso_threshold) flag += 2;
		if (height_field_node(i + 1, j + 1) > iso_threshold) flag += 4;
		if (height_field_node(i, j + 1) > iso_threshold) flag += 8;

		if (flag == 15)		// full
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), vix_node(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), vix_node(i, j + 1)).getReversedCW();
		}
		else if (flag == 5) // n0 n2 (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 10)// n1 n3  (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 14)	// except n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), vix_node(i + 1, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), u_face(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j), vix_node(i + 1, j)).getReversedCW();
		}
		else if (flag == 13)	// except n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), u_face(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i + 1, j), vix_node(i + 1, j + 1)).getReversedCW();
		}
		else if (flag == 11)	// except n2
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), u_face(i + 1, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 7)		// except n3
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 1)		// n0 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 2)		// n1 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), v_face(i, j)).getReversedCW();
		}
		else if (flag == 4)		// n2 only
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), v_face(i, j + 1), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 8)		// n3 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 3) // n0 n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), u_face(i + 1, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 6) // n1 n2
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), vix_node(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i + 1, j), v_face(i, j + 1), v_face(i, j)).getReversedCW();
		}
		else if (flag == 12) // n2 n3
		{
			triangles.PushBack() = TV_INT(vix_node(i + 1, j + 1), vix_node(i, j + 1), u_face(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face(i + 1, j)).getReversedCW();
		}
		else if (flag == 9) // n3 n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
	}

	// find boundary edges to make side later
	vertices.CopyToArray(surface.vertex_positions_);
	triangles.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();

	LinkedArray<TV2_INT> boundary_edges_linked_array;
	surface.findBoundaryEdges(boundary_edges_linked_array);

	Array1D<TV2_INT> boundary_edges;
	boundary_edges_linked_array.CopyToArray(boundary_edges);

	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		const TV proj_dir = (bottom_pos - sph_center).normalizedDouble();

		vertices.PushBack() = bottom_pos - proj_dir * front_print_thickness;
	}
	vertices.CopyToArray(surface.vertex_positions_);

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
		triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
		triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
		triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);
}

void LITHOPHANE::InitializePrism(const IMAGE_2D& image, StaticTriangularSurface& surface, int num_sides, const float inner_radius, const float base_thickness, const float front_print_thickness, const float corner_angle, const int smoothing_repeat)
{
	const T outer_cylinder_radius = inner_radius + base_thickness + front_print_thickness;

	const float height = 2.0f*PI*outer_cylinder_radius;
	const float width = height / (T)image.res_y_ * (T)image.res_x_;

	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < height_field_node.i_res_ - 1; ++i)
	{
		height_field_node(i, j) = (image.data_.getIClampedJRepeated(i - 1, j - 1).GetReversedGray() + image.data_.getIClampedJRepeated(i, j - 1).GetReversedGray()
			+ image.data_.getIClampedJRepeated(i, j - 1).GetReversedGray() + image.data_.getIClampedJRepeated(i, j).GetReversedGray()) * 0.25f;
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIClampJRepeat(grid, height_field_node, smoothing_repeat);

	// cylinder shares one j-row
	surface.vertex_positions_.initialize(grid.GetNumAllNodes() - (grid.i_res_ + 1) + grid.GetNumAllNodes() - (grid.i_res_ + 1));

	// define vertices
	int ix = 0;

	// vertices of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		if (num_sides <= 0)	// cylinder
		{
			const T outer_radius = inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;
			surface.vertex_positions_[ix++] = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), outer_radius, outer_radius, width);
		}
		else // prism
		{
			const T outer_radius = inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;
			surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZLinearCorner(grid.GetNodeUV(i, j), num_sides, corner_angle, outer_radius, outer_radius, width);
		}
	}

	// all vertices of inner cylinder plane
	/*	for (int j = 0; j < height.j_res_-1; ++j)
	for (int i = 0; i < height.i_res_; ++i)
	{
//		TV v = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), inner_cylinder_radius, par.width_);
		TV v = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), inner_cylinder_radius - height(i, j)*par.back_print_thickness_, par.width_);
//		if (v.x_ < base_thickness) v = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), par.inner_candle_radius_, par.width_);

		surface.vertex_positions_[ix++] = v;
	}*/

	const int v_start_bottom = ix;	// start index of bottom line vertices of inner cylinder
	{const int i = 0;
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	{
		if(num_sides <= 0) surface.vertex_positions_[ix++] = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), inner_radius, inner_radius, width);		// cylinder
		else surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZLinearCorner(grid.GetNodeUV(i, j), num_sides, corner_angle, inner_radius, inner_radius, width);			// prism
	}}

	const int v_start_upper = ix;	// start index of upper line vertices of inner cylinder
	{const int i = height_field_node.i_res_ - 1;
	for (int j = 0; j < height_field_node.j_res_ - 1; ++j)
	{
		if (num_sides <= 0) surface.vertex_positions_[ix++] = UV_TO_XYZ::GetCylinderXYZ(grid.GetNodeUV(i, j), inner_radius, inner_radius, width);		// cylinder
		else surface.vertex_positions_[ix++] = UV_TO_XYZ::GetPrismXYZLinearCorner(grid.GetNodeUV(i, j), num_sides, corner_angle, inner_radius, inner_radius, width);			// prism
	}
	}

	LinkedArray<TV_INT> tri;

	// triangles of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	for (int i = 0; i < height_field_node.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, j + 1));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i + 1, j + 1));
	}

	int j = height_field_node.j_res_ - 2;
	for (int i = 0; i < height_field_node.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i, 0));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, 0), height_field_node.get1DIndex(i + 1, j), height_field_node.get1DIndex(i + 1, 0));
	}

	const int offset = height_field_node.i_res_*height_field_node.j_res_ - height_field_node.i_res_;		// one j row is shared

	// triangles of inner cylinder
	/*	for (int j = 0; j < height.j_res_ - 2; j++)
	for (int i = 0; i < height.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j) + offset, height.Get1DIndex(i, j + 1) + offset, height.Get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j + 1) + offset, height.Get1DIndex(i + 1, j + 1) + offset, height.Get1DIndex(i + 1, j) + offset);
	}

	j = height.j_res_ - 2;
	for (int i = 0; i < height.i_res_ - 1; i++)
	{
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j) + offset, height.Get1DIndex(i, 0) + offset, height.Get1DIndex(i + 1, j) + offset);
		tri.PushBack() = TV_INT(height.Get1DIndex(i, 0) + offset, height.Get1DIndex(i + 1, 0) + offset, height.Get1DIndex(i + 1, j) + offset);
	}
	*/

	for (int j = 0; j < height_field_node.j_res_ - 2; ++j)
	{
		tri.PushBack() = TV_INT(v_start_bottom + j, v_start_bottom + 1 + j, v_start_upper + j);
		tri.PushBack() = TV_INT(v_start_bottom + 1 + j, v_start_upper + j + 1, v_start_upper + j);
	}
	{
		int j = height_field_node.j_res_ - 2;
		tri.PushBack() = TV_INT(v_start_bottom + j, v_start_bottom, v_start_upper + j);
		tri.PushBack() = TV_INT(v_start_bottom, v_start_upper, v_start_upper + j);
	}

	//	left wall
	/*	int i = 0;

	for (int j = 0; j < height.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j) + offset, height.Get1DIndex(i, j), height.Get1DIndex(i, j + 1) + offset);
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j + 1) + offset, height.Get1DIndex(i, j), height.Get1DIndex(i, j + 1));
	}

	{
		int j = height.j_res_ - 2;

		tri.PushBack() = TV_INT(height.Get1DIndex(i, j) + offset, height.Get1DIndex(i, j), height.Get1DIndex(i, 0) + offset);
		tri.PushBack() = TV_INT(height.Get1DIndex(i, 0) + offset, height.Get1DIndex(i, j), height.Get1DIndex(i, 0));
	}*/

	{int i = 0;
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(v_start_bottom + j, height_field_node.get1DIndex(i, j), v_start_bottom + j + 1);
		tri.PushBack() = TV_INT(v_start_bottom + j + 1, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, j + 1));
	}

	int j = height_field_node.j_res_ - 2;
	tri.PushBack() = TV_INT(v_start_bottom + j, height_field_node.get1DIndex(i, j), v_start_bottom);
	tri.PushBack() = TV_INT(v_start_bottom, height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i, 0)); }

	// right wall
	/*	i = height.i_res_ - 1;

	for (int j = 0; j < height.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j) + offset, height.Get1DIndex(i, j + 1) + offset, height.Get1DIndex(i, j));
		tri.PushBack() = TV_INT(height.Get1DIndex(i, j + 1) + offset, height.Get1DIndex(i, j + 1), height.Get1DIndex(i, j));
	}

	{
		int j = height.j_res_ - 2;

		tri.PushBack() = TV_INT(height.Get1DIndex(i, j) + offset, height.Get1DIndex(i, 0) + offset, height.Get1DIndex(i, j));
		tri.PushBack() = TV_INT(height.Get1DIndex(i, 0) + offset, height.Get1DIndex(i, 0), height.Get1DIndex(i, j));
	}*/

	{int i = height_field_node.i_res_ - 1;
	for (int j = 0; j < height_field_node.j_res_ - 2; j++)
	{
		tri.PushBack() = TV_INT(v_start_upper + j, v_start_upper + j + 1, height_field_node.get1DIndex(i, j));
		tri.PushBack() = TV_INT(v_start_upper + j + 1, height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i, j));
	}

	{int j = height_field_node.j_res_ - 2;
	tri.PushBack() = TV_INT(v_start_upper + j, v_start_upper, height_field_node.get1DIndex(i, j));
	tri.PushBack() = TV_INT(v_start_upper, height_field_node.get1DIndex(i, 0), height_field_node.get1DIndex(i, j)); }}

	tri.CopyToArray(surface.triangles_);

	// rotate cylinder
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV old = surface.vertex_positions_[vix];

		surface.vertex_positions_[vix] = TV(old.y_, old.z_, old.x_);
	}

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();
}