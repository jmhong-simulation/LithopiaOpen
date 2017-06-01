// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <glm/glm.hpp>
#include <atomic>

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/LinkedArray.h"
#include "DataStructure/DynamicArray.h"
#include "DataStructure/Vector3D.h"
#include "Geometry/BOX_3D.h"

class StaticTriangleOperations
{
private:
	Array1D<TV_INT> &triangles_;
	Array1D<TV>	 &vertex_positions_;

	Array1D<TV_INT> &edge_tri_ix_of_triangles_;		// adjacent triangles of triangles

public:
	StaticTriangleOperations(Array1D<TV>& _vertex_positions, Array1D<TV_INT>& _triangles, Array1D<TV_INT>& _edge_tri_ix_of_triangles)
		: vertex_positions_(_vertex_positions), triangles_(_triangles), edge_tri_ix_of_triangles_(_edge_tri_ix_of_triangles)
	{}

	TV  getNormal(const int& tri_index) const;
	TV  getNormalDouble(const int& tri_index) const;

	TV& getVertexPosition(const int& tri_ix, const int& v_ix) const;
	
	TV& getOppositeVertexPosOfEdgeTriangle(const int& tri_ix, const int& edge_number) const;
	TV& getOppositeVertexPos(const int& tri_ix, const int& vi, const int& vj) const;
	
	TV  getButterFlyEdgeVertex(const int& tri_ix, const int& edge_number) const;
	TV  getLinearEdgeVertex(const int& tri_ix, const int& edge_number) const;
	TV  getLoopEdgeVertex(const int& tri_ix, const int& edge_number) const;

	int getEdgeIndex(const int& tri_ix, const int& vi, const int& vj) const;
	int getEdgeIndexOfEdgeTriangle(const int& tri_ix, const int& edge_number) const;
	int getVertexIndex(const int& tri_ix, const int& v_ix) const;

	bool containsVertices(const int& tri_ix, const int& vi, const int& vj) const;
	bool hasDuplicatedEdgeTriangles(const int& tri_ix) const;

	T getSqrEdgeLength(const int& tri_ix, const int& edge_number) const;
	int getShortEdgeIndex(const int& tri_ix, const T& sqr_min_edge_length) const;

	T getVoronoiArea(const int& tri_ix, const int& v0) const;
	void calculateMeanCurvatureHelper(const int& tri_ix, const int& v0, T& area_mix_sum, TV& mean_curvature_sum) const;

	BOX_3D<T> getAABB(const int& tri_ix) const;

	bool checkFlip(const int& tri_ix, const int& v_chage, const TV& new_pos) const;	// check if this triangle flips when v_ix = v_change moves to new_pos
};
