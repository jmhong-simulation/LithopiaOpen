// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <atomic>
#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"
#include "DataStructure/Array2D.h"
#include "DataStructure/GridUniform2D.h"
#include "Geometry/DynamicContour2D.h"
#include "Parallelism/MultiThreading.h"
#include "DataStructure/LinkedArray.h"

class MarchingTriangles
{
public:
public:
	GridUniform2D cell_grid_, x_edge_grid_, y_edge_grid_, z_edge_grid_;
	Array2D<int> x_edge_vertices_, y_edge_vertices_;

	void polygonize(MT* mt, const int thread_id, const GridUniform2D& node_grid, const Array2D<T>& node_phi_array_, DynamicContour2D& contour);

	void generateLines(MT* mt, const int thread_id, const Array2D<T>& node_phi_array_, DynamicContour2D& contour_);
	void generateEdgeVertices(MT* mt, const int thread_id, const GridUniform2D& node_grid, const GridUniform2D& u_grid_, const Array2D<T>& node_phi_array_, const TV_INT ix_dev, int& num_vertices, LinkedArray<TV2>& vertex_buffer, LinkedArray<TV2>& normal_buffer, Array2D<int>& u_vertices_);
	void generateVertices(MT* mt, const int thread_id, const GridUniform2D& node_grid, const Array2D<T>& node_phi_array_, DynamicContour2D& contour);

	bool isIntersecting(const T phi0, const T phi1)
	{
		if (phi0 <= (T)0 && phi1 > (T)0) return true;
		else if (phi1 <= (T)0 && phi0 > (T)0) return true;
		else return false;
	}
};