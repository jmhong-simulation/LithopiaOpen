#pragma once

#include "../GENERIC_DEFINITIONS.h"
#include <glm/glm.hpp>
#include <atomic>
#include "../DataStructure/LinkedArray.h"
#include "../DataStructure/DynamicArray.h"
#include "../DataStructure/Vector3D.h"

class VERTEX_OPERATIONS
{
private:
	Array1D<TV> &vertex_positions_;			// vertex position Array1D
	Array1D<TV> &vertex_normals_;

	Array1D<int> &start_ix_adj_tri_of_vertices_;		// 1. count number of adjacent neighbors, 2. store index of neighboring triangles index

public:
	VERTEX_OPERATIONS(Array1D<TV>& _vertex_positions, Array1D<TV>& _vertex_normals, Array1D<int>& _start_ix_adj_tri_of_vertices)
		: vertex_positions_(_vertex_positions), vertex_normals_(_vertex_normals), start_ix_adj_tri_of_vertices_(_start_ix_adj_tri_of_vertices)
	{}

public:
	const TV& GetPosition(const int& v_ix) const;
	const int GetNumTriangles(const int& v_ix) const;
};