#include "VERTEX_OPERATIONS.h"

const TV& VERTEX_OPERATIONS::GetPosition(const int& v_ix) const
{
	return vertex_positions_[v_ix];
}

const int VERTEX_OPERATIONS::GetNumTriangles(const int& v_ix) const
{
	return start_ix_adj_tri_of_vertices_[v_ix + 1] - start_ix_adj_tri_of_vertices_[v_ix] + 1;
}
