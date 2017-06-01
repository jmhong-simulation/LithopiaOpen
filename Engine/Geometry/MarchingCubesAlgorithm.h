// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <atomic>
#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"
#include "DataStructure/Array2D.h"
#include "DataStructure/Array3D.h"
#include "DataStructure/GridUniform3D.h"
#include "Geometry/StaticTriangularSurface.h"
#include "Parallelism/MultiThreading.h"

const static float node_lookup[8][3]={{-1,1,-1},{1,1,-1},{1,-1,-1},{-1,-1,-1},{-1,1,1},{1,1,1},{1,-1,1},{-1,-1,1},};
static const int Nodes_of_Edges[12][2]={{0,1},{1,2},{2,3},{0,3},{4,5},{5,6},{6,7},{7,4},{4,0},{5,1},{2,6},{3,7}};//[edge_number][node_number]
static const int power_of_two[15]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384};

//References: http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/ 

class MarchingCubesAlgorithm
{
public:
	GridUniform3D cell_grid_, x_edge_grid_, y_edge_grid_, z_edge_grid_;
	Array3D<int> x_edge_vertices_, y_edge_vertices_, z_edge_vertices_;

	void polygonize(MT* mt, const int thread_id, const GridUniform3D& node_grid, const Array3D<T>& node_phi_array_, StaticTriangularSurface &surface_);

	void generateTriangles(MT* mt, const int thread_id, const Array3D<T>& node_phi_array_, StaticTriangularSurface &surface_);
	void generateEdgeVertices(MT* mt, const int thread_id, const GridUniform3D& node_grid, const GridUniform3D& u_grid_, const Array3D<T>& node_phi_array_, const TV_INT ix_dev, int& num_vertices, LinkedArray<TV>& vertex_buffer, LinkedArray<TV>& normal_buffer, Array3D<int>& u_vertices_);
	void generateVertices(MT* mt, const int thread_id, const GridUniform3D& node_grid, const Array3D<T>& node_phi_array_, StaticTriangularSurface &surface_);	// generate edge vertices and gather them

	bool isIntersecting(const T phi0, const T phi1);
};//End of Class MARCHING_CUBES_ALGORITHM
