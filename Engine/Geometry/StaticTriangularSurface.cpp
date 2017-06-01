#include <iostream>
#include <fstream>
#include <vector>

#include "CONVENTIONAL_MACROS.h"
#include "StaticTriangularSurface.h"
#include "Intersection.h"

using namespace std;

/*
void TriangularSurfaceStatic::ReadTRI(const char *filename)
{
	std::cout << "Start reading TRI file " << filename << std::endl;
	
	int num_vertices(0), num_triangles(0);

	// count number of vertices and triangles
	{ifstream file(filename);
	char c[255];

	while (true)
	{
		file >> c;
		if (file.eof() != 0) break;
		if (strcmp(c, "#") == 0) file.getline(c, 255);	// comments
		else if (strcmp(c, "v") == 0) num_vertices++;	// count vertices
		else if (strcmp(c, "f") == 0) num_triangles++;	// count triangles
	}

	file.close();}

	std::cout << num_vertices << " vertices and " << num_triangles << " triangles are found." << std::endl;

	// read data
	vertex_positions_1d_.Initialize(num_vertices);
	triangles_1d_.Initialize(num_triangles);

	{ifstream file(filename);
	char c[255];		//TODO: find a way not to open the same file twice

	while (true)
	{
		file >> c;
		if (file.eof() != 0) break;
		if (strcmp(c, "#") == 0)		// comments
		{
			file.getline(c, 255);
		}
		else if (strcmp(c, "v") == 0)	// vertex
		{
			T x[3];
			file >> x[0] >> x[1] >> x[2];
			vertex_positions_1d_.Push(TV(x));
		}
		else if (strcmp(c, "f") == 0)
		{
			int v[3];
			file >> v[0] >> v[1] >> v[2];
//			v[0]--; v[1]--; v[2]--;			//NOTE: SMF indices are 1 based and TRI indices are 0 based
			triangles_1d_.Push(TV_INT(v));
		}
	}

	file.close(); }

	std::cout << "Reading complete." << std::endl;
}
*/

/*
void TriangularSurfaceStatic::ReadSMF(const char *filename)
{
	std::cout << "Start reading SMF file " << filename << std::endl;

	int num_vertices(0), num_triangles(0);

	// count number of vertices and triangles
	{ifstream file(filename);
	char c[255];

	while (true)
	{
		file >> c;
		if (file.eof() != 0) break;
		if (strcmp(c, "#") == 0) file.getline(c, 255);	// comments
		else if (strcmp(c, "v") == 0) num_vertices++;	// count vertices
		else if (strcmp(c, "f") == 0) num_triangles++;	// count triangles
	}

	file.close(); }

	std::cout << num_vertices << " vertices and " << num_triangles << " triangles are found." << std::endl;

	// read data
	vertex_positions_1d_.Initialize(num_vertices);
	triangles_1d_.Initialize(num_triangles);

	{ifstream file(filename);
	char c[255];		//TODO: find a way not to open the same file twice

	while (true)
	{
		file >> c;
		if (file.eof() != 0) break;
		if (strcmp(c, "#") == 0) file.getline(c, 255);// comments
		else if (strcmp(c, "v") == 0)	// vertex
		{
			T x[3];
			file >> x[0] >> x[1] >> x[2];
			vertex_positions_1d_.Push(TV(x));
		}
		else if (strcmp(c, "f") == 0)
		{
			TV_INT v;
			file >> v.x_ >> v.y_ >> v.z_;
			v.x_--; v.y_--; v.z_--;			//NOTE: SMF indices are 1 based and TRI indices are 0 based
			triangles_1d_.Push(v);
		}
	}

	file.close(); }

	std::cout << "Reading complete." << std::endl;
}
*/

/*
void TriangularSurfaceStatic::ReadOBJ(const char *filename)
{
	std::cout << "Start reading OBJ file " << filename << std::endl;

	bool flag_vt(false);
	bool flag_vn(false);

	int num_vertices = 0, num_triangles = 0;

	// count # of vertices and triangles
	ifstream file(filename);

	// check if file is correctly opened
	if (file.is_open() == false){std::cout << filename << " does not exist. Program terminated." << std::endl;exit(-1);}

	char c[255];

	while (true)
	{
		file >> c;
		if (file.eof() != 0) break;		// end of file
		if (strcmp(c, "#") == 0) file.getline(c, 255); // comments (less than 255 characters)
		else if (strcmp(c, "v") == 0) num_vertices++;  // count number of vertices
		else if (strcmp(c, "vt") == 0) flag_vt = true;
		else if (strcmp(c, "vn") == 0) flag_vn = true;
		else if (strcmp(c, "f") == 0) num_triangles++; // count number of triangles
	}
	file.clear();
	file.close();

	std::cout << num_vertices << " vertices and " << num_triangles << " triangles are found." << std::endl;

	// read data from here
	file.open(filename);	//TODO: find a way not to open the file twice.

	vertex_positions_1d_.Initialize(num_vertices);
	triangles_1d_.Initialize(num_triangles);

	// check if file is correctly opened
	if (file.is_open() == false){ std::cout << filename << " does not exist. Program terminated." << std::endl; exit(-1); }

	while (true)
	{
		file >> c;

		if (file.eof() != 0) break;		// finish reading if file is ended

		if (strcmp(c, "#") == 0) file.getline(c, 255); // comments (less than 255 characters)
		else if (strcmp(c, "v") == 0) // vertices
		{
			TV &vertex_pos(vertex_positions_1d_.Push());
			file >> vertex_pos.x_ >> vertex_pos.y_ >> vertex_pos.z_;

//			std::cout << vertex_pos.x_ << " " << vertex_pos.y_ << " " << vertex_pos.z_ << std::endl;
		}
		else if (strcmp(c, "vt") == 0) {} //TODO: read texture coordinate
		else if (strcmp(c, "vn") == 0) {} //TODO: read vertex normals
		else if (strcmp(c, "f") == 0)
		{
			int v[3], vt[3], vn[3];
			if (flag_vt == true && flag_vn == true)
			{
				for (int i = 0; i < 3; i++)
				{
					file >> v[i]; file.get(c, 2);
					file >> vt[i]; file.get(c, 2);
					file >> vn[i];

					v[i]--;
					vt[i]--;
					vn[i]--;
				}
			}
			else if (flag_vt == false && flag_vn == true)
			{
				for (int i = 0; i < 3; i++)
				{
					file >> v[i]; file.get(c, 2); file.get(c, 2);
					file >> vn[i];
					v[i]--;
					vn[i]--;
				}
			}
			else if (flag_vt == false && flag_vn == false)
			{
				for (int i = 0; i < 3; i++)
				{
					file >> v[i];
					v[i]--;
				}
			}

			triangles_1d_.Push(TV_INT(v));
	
//			std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
		}
	}
	file.clear();
	file.close();

	std::cout << "Reading complete." << std::endl;
}
*/

void StaticTriangularSurface::findAdjacentTrianglesOfVertices()
{
	const int num_vertices = vertex_positions_.num_elements_;
	const int num_triangles = triangles_.num_elements_;

	start_ix_adj_tri_of_vertices_.initialize(num_vertices+1, 0);		// last value will be used to find the range of adjacent triangles of the last vertex

	// count the number of triangles connected to each vertex and store the number in ver_index_ Array1D
	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const TV_INT &vertex_indices(triangles_[tri_ix]);

		start_ix_adj_tri_of_vertices_[vertex_indices.v0_] ++;
		start_ix_adj_tri_of_vertices_[vertex_indices.v1_] ++;
		start_ix_adj_tri_of_vertices_[vertex_indices.v2_] ++;
	}

	// find the index of 1D adjacent triangle indices
	int accumulated_index = 0;
	for (int ver_ix = 0; ver_ix < num_vertices; ver_ix++)
	{
		const int temp = start_ix_adj_tri_of_vertices_[ver_ix];

		start_ix_adj_tri_of_vertices_[ver_ix] = accumulated_index;
		accumulated_index += temp;
	}

	start_ix_adj_tri_of_vertices_[num_vertices] = accumulated_index;	// last value helps the adj triangles iteration of last vertex

	adj_tri_ix_of_vertices_.initialize(accumulated_index);	// this is the number of all connected triangles

	Array1D<int> tri_count(start_ix_adj_tri_of_vertices_);					// this is a counter TODO: atomic

	// store the indices of adjacent triangles
	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const TV_INT &vertex_indices(triangles_[tri_ix]);

		adj_tri_ix_of_vertices_[tri_count[vertex_indices[0]]++] = tri_ix;
		adj_tri_ix_of_vertices_[tri_count[vertex_indices[1]]++] = tri_ix;
		adj_tri_ix_of_vertices_[tri_count[vertex_indices[2]]++] = tri_ix;
	}
}

void StaticTriangularSurface::determineVertexMeanCurvatures()
{
	// Discrete Differential-Geometry Operators	for Triangulated 2 - Manifolds
	// http://http.cs.berkeley.edu/~jrs/meshpapers/MeyerDesbrunSchroderBarr.pdf

	const int num_vertices = vertex_positions_.num_elements_;
	vertex_mean_curvatures_.initialize(num_vertices);

	for (int ver_ix = 0; ver_ix < num_vertices; ver_ix++)
	{
		// iterate adjacent triangles
		T area_mix_sum = 0.0f;
		TV mean_curvature_sum;
		const int tri_start = start_ix_adj_tri_of_vertices_[ver_ix], tri_end = start_ix_adj_tri_of_vertices_[ver_ix + 1];	// Note the end +1 of ver_index_.Initialize(num_vertices+1) in FindNeighboringTriangles()
		for (int tri = tri_start; tri < tri_end; tri++)
		{
			tri_ops_.calculateMeanCurvatureHelper(adj_tri_ix_of_vertices_[tri], tri_ops_.getVertexIndex(adj_tri_ix_of_vertices_[tri], ver_ix), area_mix_sum, mean_curvature_sum);
		}

		assert(area_mix_sum > 0.0f);

		vertex_mean_curvatures_[ver_ix] = mean_curvature_sum / area_mix_sum * 0.5f;
	}
}

void StaticTriangularSurface::determineFaceAveragedVertexNormals()
{
	const int num_vertices = vertex_positions_.num_elements_;
	vertex_normals_.initialize(num_vertices);

	for (int ver_ix = 0; ver_ix < num_vertices; ver_ix++)
	{
		TV &normal(vertex_normals_[ver_ix]);

		// iterate adjacent triangles
		const int tri_start = start_ix_adj_tri_of_vertices_[ver_ix], tri_end = start_ix_adj_tri_of_vertices_[ver_ix + 1];	// Note the end +1 of ver_index_.Initialize(num_vertices+1) in FindNeighboringTriangles()

		if (tri_start == tri_end)
		{
			// no neighbor (TODO: debug)
			normal = TV(0, 0, 1);

//			std::cout << "No neighbor " << std::endl;
		}
		else for (int tri = tri_start; tri < tri_end; tri++)
		{
//			normal += tri_ops_.GetNormal(adj_tri_ix_of_vertices_[tri])*tri_ops_.GetVoronoiArea(adj_tri_ix_of_vertices_[tri], tri_ops_.GetVertexIndex(adj_tri_ix_of_vertices_[tri], ver_ix));
			normal += tri_ops_.getNormalDouble(adj_tri_ix_of_vertices_[tri]);

// 			if (tri_ops_.GetVoronoiArea(adj_tri_ix_of_vertices_[tri], tri_ops_.GetVertexIndex(adj_tri_ix_of_vertices_[tri], ver_ix)) < 0)
// 			{std::cout << "Negative VoronoiArea" << std::endl;exit(1);}
		}

		normal.normalizeDouble();

//		std::cout << normal << std::endl;
	}
}

void StaticTriangularSurface::findBoundaryEdges(LinkedArray<TV2_INT>& boundary_edges)
{
	const int num_triangles = triangles_.num_elements_;

	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		TV_INT &adj_triangles(edge_tri_ix_of_triangles_[tri_ix]);

		for (int d = 0; d < 3; ++d)
			if (adj_triangles[d] == -1)
				boundary_edges.PushBack() = TV2_INT(triangles_[tri_ix][(d + 1) % 3], triangles_[tri_ix][(d + 2) % 3]);
	}
}

void StaticTriangularSurface::findEdgeTrianglesOfTriangles()
{
	const int num_triangles = triangles_.num_elements_;

	edge_tri_ix_of_triangles_.initialize(num_triangles, TV_INT(-1,-1,-1));

	int num_of_no_neighbor_edges = 0;

	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		TV_INT &adj_triangles(edge_tri_ix_of_triangles_[tri_ix]);

// 		if (triangles_[tri_ix].v0_ == -1)
// 		{
// 			adj_triangles = TV_INT(-1, -1, -1);
// 			continue;
// 		}

		const TV_INT &tri_vertices(triangles_[tri_ix]);
		const int v0(tri_vertices.v0_), v1(tri_vertices.v1_), v2(tri_vertices.v2_);

		// find edge triangle 1 and 2 by searching adjacent triangles of vertex 0
		{const int tri_start = start_ix_adj_tri_of_vertices_[v0], tri_end = start_ix_adj_tri_of_vertices_[v0 + 1];
		for (int tri = tri_start; tri < tri_end; tri++)
		{
			const int v_tri_ix = adj_tri_ix_of_vertices_[tri];

			if (v_tri_ix == tri_ix) continue;		// skip current triangle

			if (tri_ops_.containsVertices(v_tri_ix, v0, v1))
			{
				if (adj_triangles.v2_ != -1) num_of_no_neighbor_edges ++;// std::cout << "Warning: more than one edge triangles" << std::endl;

				adj_triangles.v2_ = v_tri_ix;
			}
//			else if (tri_ops_.ContainsVertices(v_tri_ix, v0, v2))
			if (tri_ops_.containsVertices(v_tri_ix, v0, v2))
			{
				if (adj_triangles.v1_ != -1) num_of_no_neighbor_edges ++; // std::cout << "Warning: more than one edge triangles" << std::endl;

				adj_triangles.v1_ = v_tri_ix;
			}
		}}

		// find edge triangle 0 by searching adjacent triangles of vertex 0
		{const int tri_start = start_ix_adj_tri_of_vertices_[v1], tri_end = start_ix_adj_tri_of_vertices_[v1 + 1];
		for (int adj_tri_ix = tri_start; adj_tri_ix < tri_end; adj_tri_ix++)
		{
			const int v_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

			if (v_tri_ix == tri_ix) continue;		// skip current triangle

			if (tri_ops_.containsVertices(v_tri_ix, v1, v2))
			{
				if (adj_triangles.v0_ != -1) num_of_no_neighbor_edges ++;//  std::cout << "Warning: more than one edge triangles" << std::endl;

				adj_triangles.v0_ = v_tri_ix;
				break;
			}
		}}

		//TODO: option
// 		if (adj_triangles.values_[0] == -1)
// 		{ 
// 			std::cout << "edge triangle not found" << std::endl;
// 			exit(-1);
// 		}
// 
// 		if (adj_triangles.values_[1] == -1)
// 		{
// 			std::cout << "edge triangle not found" << std::endl;
// 			exit(-1);
// 		}
// 
// 		if (adj_triangles.values_[2] == -1)
// 		{
// 			std::cout << "edge triangle not found" << std::endl;
// 			exit(-1);
// 		}
	}

	if (num_of_no_neighbor_edges > 0) std::cout << "Warning: number of edges which have no neighboring triangles " << num_of_no_neighbor_edges << std::endl;
}

TV StaticTriangularSurface::getEdgeVertex(const int mode, const int& tri_ix, const int& edge_number) const
{
	switch (mode)
	{
	case 0:
		return tri_ops_.getLinearEdgeVertex(tri_ix, edge_number);
	case 1:
		return tri_ops_.getButterFlyEdgeVertex(tri_ix, edge_number);
	case 2:
		return tri_ops_.getLoopEdgeVertex(tri_ix, edge_number);
	default:
		assert(false);
		return Vector3D<T>();
	}
}

void StaticTriangularSurface::smoothVertexPositionsLaplacian()
{
	vertex_positions_temp_.initialize(vertex_positions_);			// for loop subdivision or Jacobi iterations

	const int num_vertices = vertex_positions_.num_elements_;

	for (int v_ix = 0; v_ix < num_vertices; v_ix++)
	{
		int num_adj_vertices = 0;
		TV  pos_adj_vertices;

		// find edge triangle 1 and 2 by searching adjacent triangles of vertex 0
		const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
		for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
		{
			const int v_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];
			const TV_INT &v_tri_ix_v(triangles_[v_tri_ix]);

			for (int d = 0; d < 3; d++)
			{
				if (v_ix != v_tri_ix_v[d])
				{
					const int edge_tri_ix_of_v_tri = edge_tri_ix_of_triangles_[v_tri_ix][tri_ops_.getEdgeIndex(v_tri_ix, v_ix, v_tri_ix_v[d])];
					if (edge_tri_ix_of_v_tri > v_tri_ix || edge_tri_ix_of_v_tri == -1)	// or -1
					{
						pos_adj_vertices += vertex_positions_temp_[v_tri_ix_v[d]];
						num_adj_vertices++;
					}
				}
			}
		}	

		if (num_adj_vertices < 3) continue;

		vertex_positions_[v_ix] = pos_adj_vertices / (T)num_adj_vertices;
	}

	vertex_positions_temp_.freeMemory();
}

void StaticTriangularSurface::smoothVertexPositionsLoop()
{
	vertex_positions_temp_.initialize(vertex_positions_);			// for loop subdivision or Jacobi iterations

	const int num_vertices = vertex_positions_.num_elements_;

	for (int v_ix = 0; v_ix < num_vertices; v_ix++)
	{
		int num_adj_vertices = 0;
		TV  pos_adj_vertices;

		// find edge triangle 1 and 2 by searching adjacent triangles of vertex 0
		const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
		for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
		{
			const int v_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];
			const TV_INT &v_tri_ix_v(triangles_[v_tri_ix]);

			for (int d = 0; d < 3; d++)
			{
				if (v_ix != v_tri_ix_v[d])
				{
					const int edge_tri_ix_of_v_tri = edge_tri_ix_of_triangles_[v_tri_ix][tri_ops_.getEdgeIndex(v_tri_ix, v_ix, v_tri_ix_v[d])];
					if (edge_tri_ix_of_v_tri > v_tri_ix || edge_tri_ix_of_v_tri == -1)	// or -1
					{
						pos_adj_vertices += vertex_positions_temp_[v_tri_ix_v[d]];
						num_adj_vertices++;
					}
				}
			}
		}

		T beta;

		if (num_adj_vertices > 3)
		{
			beta = 3.0f / (8.0f*(float)num_adj_vertices);
		}
		else if (num_adj_vertices == 3)
		{
			beta = 3.0f / 16.0f;
		}
		else
		{
			beta = 0.0f;
		}

		vertex_positions_[v_ix] = pos_adj_vertices*beta + (1.0f - beta*(float)num_adj_vertices)*vertex_positions_[v_ix];
	}

	vertex_positions_temp_.freeMemory();
}

void StaticTriangularSurface::generateEdgeVertices(const int mode)
{
	new_vertex_postitions_.Reset();

	edge_v_ix_of_triangles_.initialize(triangles_.num_elements_, TV_INT(-1, -1, -1));

	const int ver_ix_offset = vertex_positions_.num_elements_;
	const int num_triangles = triangles_.num_elements_;

	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const int *vix = triangles_.values_[tri_ix].values_;
		const int *etri = edge_tri_ix_of_triangles_.values_[tri_ix].values_;

		// Generate edge vertices
		// don't generate edge vertex if edge triangle already has it.
		for (int d = 0; d < 3; d++)
		{
			if (etri[d] != -1 && etri[d] < tri_ix) continue;		// let earlier triangle generate edge vertex

			const int new_evix = ver_ix_offset + new_vertex_postitions_.num_elements_;		// index of new edge vertex
			new_vertex_postitions_.PushBack() = getEdgeVertex(mode, tri_ix, d);
			edge_v_ix_of_triangles_.values_[tri_ix].values_[d] = new_evix;

			if (etri[d] != -1)
				edge_v_ix_of_triangles_.values_[etri[d]].values_[tri_ops_.getEdgeIndex(etri[d], vix[(d + 1) % 3], vix[(d + 2) % 3])] = new_evix;		// update the edge vertex of edge triangle
		}
	}
}

void StaticTriangularSurface::generateEdgeVerticesPlaneCut(const PLANE& cut_plane)
{
	new_vertex_postitions_.Reset();

	edge_v_ix_of_triangles_.initialize(triangles_.num_elements_, TV_INT(-1, -1, -1));

	const int ver_ix_offset = vertex_positions_.num_elements_;
	const int num_triangles = triangles_.num_elements_;

	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const int *vix = triangles_.values_[tri_ix].values_;
		const int *etri = edge_tri_ix_of_triangles_.values_[tri_ix].values_;

		// Generate edge vertices
		// don't generate edge vertex if edge triangle already has it.
		for (int d = 0; d < 3; d++)
		{
			if (etri[d] != -1 && etri[d] < tri_ix) continue;		// let earlier triangle generate edge vertex

			// check if this edge intersects the plane
			T t;
			TV q;
			if (Intersection::checkLinePlane(vertex_positions_[vix[(d + 1) % 3]], vertex_positions_[vix[(d + 2) % 3]], cut_plane, t, q) == 0) continue;

			const int new_evix = ver_ix_offset + new_vertex_postitions_.num_elements_;		// index of new edge vertex
			new_vertex_postitions_.PushBack() = q;
			edge_v_ix_of_triangles_.values_[tri_ix].values_[d] = new_evix;

			if (etri[d] != -1)
				edge_v_ix_of_triangles_.values_[etri[d]].values_[tri_ops_.getEdgeIndex(etri[d], vix[(d + 1) % 3], vix[(d + 2) % 3])] = new_evix;		// update the edge vertex of edge triangle
		}
	}
}

void StaticTriangularSurface::splitTriangles()
{
	new_triangles_.Reset();

	const int num_triangles = triangles_.num_elements_;

	// split triangles
	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const int e0_ix = edge_v_ix_of_triangles_[tri_ix].v0_;
		const int e1_ix = edge_v_ix_of_triangles_[tri_ix].v1_;
		const int e2_ix = edge_v_ix_of_triangles_[tri_ix].v2_;
		
		const int v0_ix = triangles_.values_[tri_ix].values_[0];
		const int v1_ix = triangles_.values_[tri_ix].values_[1];
		const int v2_ix = triangles_.values_[tri_ix].values_[2];

		if (e0_ix != -1 && e1_ix != -1 && e2_ix != -1)	// 1 to 4 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v0_ix, e2_ix, e1_ix);
			new_triangles_.PushBack() = TV_INT(e2_ix, v1_ix, e0_ix);
			new_triangles_.PushBack() = TV_INT(e1_ix, e0_ix, v2_ix);
			new_triangles_.PushBack() = TV_INT(e2_ix, e0_ix, e1_ix);
		}
		else if (e0_ix != -1 && e1_ix == -1 && e2_ix == -1)	// edge 0 only, 1 to 2 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v0_ix, v1_ix, e0_ix);
			new_triangles_.PushBack() = TV_INT(v0_ix, e0_ix, v2_ix);
		}
		else if (e0_ix == -1 && e1_ix != -1 && e2_ix == -1)	// edge 1 only, 1 to 2 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v0_ix, v1_ix, e1_ix);
			new_triangles_.PushBack() = TV_INT(v1_ix, v2_ix, e1_ix);
		}
		else if (e0_ix == -1 && e1_ix == -1 && e2_ix != -1)	// edge 2 only, 1 to 2 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v0_ix, e2_ix, v2_ix);
			new_triangles_.PushBack() = TV_INT(v1_ix, v2_ix, e2_ix);
		}
		else if (e0_ix == -1 && e1_ix != -1 && e2_ix != -1)	// not edge 0 only, 1 to 3 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v0_ix, e2_ix, e1_ix);

			//TODO: ambiguity
			new_triangles_.PushBack() = TV_INT(e1_ix, e2_ix, v1_ix);
			new_triangles_.PushBack() = TV_INT(e1_ix, v1_ix, v2_ix);
		}
		else if (e0_ix != -1 && e1_ix == -1 && e2_ix != -1)	// not edge 1 only, 1 to 3 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v1_ix, e0_ix, e2_ix);

			//TODO: ambiguity
			new_triangles_.PushBack() = TV_INT(e2_ix, e0_ix, v0_ix);
			new_triangles_.PushBack() = TV_INT(v0_ix, e0_ix, v2_ix);
		}
		else if (e0_ix != -1 && e1_ix != -1 && e2_ix == -1)	// not edge 2 only, 1 to 3 subdivision
		{
			new_triangles_.PushBack() = TV_INT(v2_ix, e1_ix, e0_ix);

			//TODO: ambiguity
			new_triangles_.PushBack() = TV_INT(e1_ix, v0_ix, e0_ix);
			new_triangles_.PushBack() = TV_INT(e0_ix, v0_ix, v1_ix);
		}
		else // not divided
		{
			new_triangles_.PushBack() = TV_INT(v0_ix, v1_ix, v2_ix);
		}
	}
}

void StaticTriangularSurface::updateVertexAndTriangleLists()
{
	//TODO: selective subdivision

	const int ver_ix_offset = vertex_positions_.num_elements_;

	// update vertex_positions and triangles
	TV *merged_vertex_positions = new TV[vertex_positions_.num_elements_ + new_vertex_postitions_.num_elements_];

	for (int ver_ix = 0; ver_ix < vertex_positions_.num_elements_; ver_ix++)
		merged_vertex_positions[ver_ix] = vertex_positions_.values_[ver_ix];	//TODO: use memcpy or realloc

	new_vertex_postitions_.CopyToPartialArray(merged_vertex_positions, ver_ix_offset);

	SWAP(vertex_positions_.values_, merged_vertex_positions, TV*);
	delete [] merged_vertex_positions;

	vertex_positions_.num_elements_ = vertex_positions_.num_elements_ + new_vertex_postitions_.num_elements_;

	new_triangles_.CopyToArray(triangles_);		// non-divided triangles are already inside new_triangles_.
}

void StaticTriangularSurface::advanceOneSubdivisionStep(const int mode)
{
	findAdjacentTrianglesOfVertices();
	findEdgeTrianglesOfTriangles();
	applySubdivision(mode);
	findAdjacentTrianglesOfVertices();
	determineFaceAveragedVertexNormals();
}

void StaticTriangularSurface::applyLaplacianSmoothing()
{
	findAdjacentTrianglesOfVertices();
	findEdgeTrianglesOfTriangles();

	smoothVertexPositionsLaplacian();

	findAdjacentTrianglesOfVertices();
	determineFaceAveragedVertexNormals();
}

void StaticTriangularSurface::applySubdivision(const int mode)
{
	generateEdgeVertices(mode);		// add new edge vertices to new_vertex_postitions_

	if(mode == 2) smoothVertexPositionsLoop();

	splitTriangles();				// add to new_triangles_ (non-divided triangles are also copied to new_triangles)

	updateVertexAndTriangleLists();	// update vertex_positions_ and triangles_
}

void StaticTriangularSurface::applyPlaneCutSuvdivision(const PLANE& cut_plane)
{
	generateEdgeVerticesPlaneCut(cut_plane);

	splitTriangles();				// add to new_triangles_ (non-divided triangles are also copied to new_triangles)

	updateVertexAndTriangleLists();	// update vertex_positions_ and triangles_
}

void StaticTriangularSurface::copyRenderingData(Array1D<TV>& positions, Array1D<TV>& normals) const
{
	const int num_triangles = triangles_.num_elements_;

	positions.initialize(num_triangles * 3);
	normals.initialize(num_triangles * 3);

	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const TV_INT &vertices(triangles_[tri_ix]);

		const int tri_ix_three = 3 * tri_ix;

		positions[tri_ix_three] = vertex_positions_[vertices.v0_];
		positions[tri_ix_three + 1] = vertex_positions_[vertices.v1_];
		positions[tri_ix_three + 2] = vertex_positions_[vertices.v2_];

		normals[tri_ix_three] = vertex_normals_[vertices.v0_];
		normals[tri_ix_three + 1] = vertex_normals_[vertices.v1_];
		normals[tri_ix_three + 2] = vertex_normals_[vertices.v2_];
	}
}

void StaticTriangularSurface::copyRenderingData(LinkedArray<TV>& positions, LinkedArray<TV>& normals) const
{
	const int num_triangles = triangles_.num_elements_;

	positions.Reset();
	normals.Reset();

	for (int tri_ix = 0; tri_ix < num_triangles; tri_ix++)
	{
		const TV_INT &vertices(triangles_[tri_ix]);

		if (vertices.v0_ == -1) continue;

		positions.PushBack() = vertex_positions_[vertices.v0_];
		positions.PushBack() = vertex_positions_[vertices.v1_];
		positions.PushBack() = vertex_positions_[vertices.v2_];

		if (use_face_normal_ == false)	// use vertex normal
		{
			normals.PushBack() = vertex_normals_[vertices.v0_];
			normals.PushBack() = vertex_normals_[vertices.v1_];
			normals.PushBack() = vertex_normals_[vertices.v2_];
		}
		else
		{
			const TV face_normal = tri_ops_.getNormalDouble(tri_ix);

			normals.PushBack() = face_normal;
			normals.PushBack() = face_normal;
			normals.PushBack() = face_normal;
		}
	}
}

void StaticTriangularSurface::copyRenderingDataShortEdges(LinkedArray<TV>& positions, const T min_edge_length) const
{
	const T sqr_length = POW2(min_edge_length);

	for (int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
	{
		const TV_INT &vertices(triangles_[tri_ix]);
		
		if (vertices.v0_ == -1) continue;

		const int edge_ix = tri_ops_.getShortEdgeIndex(tri_ix, sqr_length);

		if (edge_ix == -1) continue;

		positions.PushBack() = vertex_positions_[vertices.values_[(edge_ix + 1) % 3]];
		positions.PushBack() = vertex_positions_[vertices.values_[(edge_ix + 2) % 3]];
	}
}

void StaticTriangularSurface::copyRenderingDataShortEdgeTriangles(LinkedArray<TV>& positions, LinkedArray<TV>& normals, const T min_edge_length, const T normal_offset) const
{
	const T sqr_length = POW2(min_edge_length);

	for (int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
	{
		const TV_INT &vertices(triangles_[tri_ix]);

		if (vertices.v0_ == -1) continue;

		const int edge_ix = tri_ops_.getShortEdgeIndex(tri_ix, sqr_length);

		if (edge_ix == -1) continue;

		positions.PushBack() = vertex_positions_[vertices.v0_] + normal_offset*vertex_normals_[vertices.v0_];
		positions.PushBack() = vertex_positions_[vertices.v1_] + normal_offset*vertex_normals_[vertices.v1_];
		positions.PushBack() = vertex_positions_[vertices.v2_] + normal_offset*vertex_normals_[vertices.v2_];

		normals.PushBack() = vertex_normals_[vertices.v0_];
		normals.PushBack() = vertex_normals_[vertices.v1_];
		normals.PushBack() = vertex_normals_[vertices.v2_];
	}
}

void StaticTriangularSurface::getTriangle(const int& tri_ix, StaticTriangle& triangle) const
{
	const TV_INT &vertices(triangles_[tri_ix]);

	triangle.v0_ = vertex_positions_[vertices.v0_];
	triangle.v1_ = vertex_positions_[vertices.v1_];
	triangle.v2_ = vertex_positions_[vertices.v2_];
}

void StaticTriangularSurface::getTriangleUV(const int& tri_ix, TV2 uv_arr[3]) const
{
	const TV_INT &vertices(uv_triangles_[tri_ix]);

	uv_arr[0] = vertex_uv_[vertices.v0_];
	uv_arr[1] = vertex_uv_[vertices.v1_];
	uv_arr[2] = vertex_uv_[vertices.v2_];
}

void StaticTriangularSurface::copyRenderingDataHighCurvatureVertices(LinkedArray<TV>& positions, LinkedArray<TV>& colors, const T kappa, const T normal_offset) const
{
	for (int v_ix = 0; v_ix < vertex_positions_.num_elements_; v_ix++)
	{
		if (vertex_mean_curvatures_[v_ix].getMagnitude() > kappa)
		{
			positions.PushBack() = vertex_positions_[v_ix] + normal_offset*vertex_normals_[v_ix];
			colors.PushBack() = TV(1, 1, 0);
		}
	}
}

void StaticTriangularSurface::writeCTM(const char *filename)
{
	CTMcontext context = ctmNewContext(CTM_EXPORT);

	// Define our mesh representation to OpenCTM
	ctmDefineMesh(context, (float*)vertex_positions_.values_, vertex_positions_.num_elements_, (unsigned int*)triangles_.values_, triangles_.num_elements_, 0);

	if(vertex_uv_.num_elements_ > 0)
		ctmAddUVMap(context, (float*)vertex_uv_.values_, "litho", "1920x430_texture.bmp");

	// Save the OpenCTM file
	ctmSave(context, filename);
	// Free the context
	ctmFreeContext(context);
}

// http://www.sgh1.net/b4/cpp-read-stl-file
void StaticTriangularSurface::readFile(const char *filename)
{
	const bool use_cout = true;

	using namespace std;

	cout << "# Reading " << filename << std::endl;
	ifstream pFile(filename, ios::in | ios::binary);

	if (pFile.is_open())
	{
		//ASCII header
		char header_info[80] = "";
//		pFile.setf(ios::left);
//		pFile.width(sizeof(unsigned char) * 80);
//		pFile >> name.c_str();
		pFile.read(header_info, 80);
		if (use_cout) cout << header_info << std::endl;

		// vertex positions
		unsigned int size = 0;
		pFile.read((char*)&size, sizeof(unsigned int));
		vertex_positions_.initialize(size);
		for (unsigned int v_ix = 0; v_ix < size; ++v_ix)
		{
			TV &vp = vertex_positions_[v_ix];

			pFile.read((char*)&(vp.x_), sizeof(float));
			pFile.read((char*)&(vp.y_), sizeof(float));
			pFile.read((char*)&(vp.z_), sizeof(float));

			if (use_cout) cout << vp << std::endl;
		}

		// triangle indices
		size = 0;
		pFile.read((char*)&size, sizeof(unsigned int));
		triangles_.initialize(size);
		for (unsigned int tri_ix = 0; tri_ix < size; ++tri_ix)
		{			
			unsigned int x, y, z;

			pFile.read((char*)&(x), sizeof(unsigned int));
			pFile.read((char*)&(y), sizeof(unsigned int));
			pFile.read((char*)&(z), sizeof(unsigned int));

			triangles_[tri_ix].assign(x, y, z);

			if (use_cout) cout << triangles_[tri_ix] << std::endl;
		}

		// vertex normals
		size = 0;
		pFile.read((char*)&size, sizeof(unsigned int));
		vertex_normals_.initialize(size);
		for (unsigned int v_ix = 0; v_ix < size; ++v_ix)
		{
			TV &vn = vertex_normals_[v_ix];

			pFile.read((char*)&(vn.x_), sizeof(float));
			pFile.read((char*)&(vn.y_), sizeof(float));
			pFile.read((char*)&(vn.z_), sizeof(float));

			if (use_cout) cout << vn << std::endl;
		}

		pFile.close();

		cout << "Finished" << endl;
	}
}

// note that size and index are unsigned int
void StaticTriangularSurface::writeFile(const char *filename)
{
	using namespace std;

	cout << "# Writing " << filename << std::endl;
	ofstream pFile(filename, ios::out | ios::binary);

	if (pFile.is_open())
	{
		//ASCII header
		string name = "Custom Triangular Surface File";
		pFile.setf(ios::left);
		pFile.width(sizeof(unsigned char) * 80);
		pFile << name.c_str();

		// vertex positions
		unsigned int size = vertex_positions_.num_elements_;
		pFile.write((char*)&size, sizeof(unsigned int));
		for (unsigned int v_ix = 0; v_ix < size; ++v_ix)
		{
			const TV &vp = vertex_positions_[v_ix];

			pFile.write((char*)&(vp.x_), sizeof(float));
			pFile.write((char*)&(vp.y_), sizeof(float));
			pFile.write((char*)&(vp.z_), sizeof(float));
		}

		// triangle indices
		size = triangles_.num_elements_;
		pFile.write((char*)&size, sizeof(unsigned int));
		for (unsigned int tri_ix = 0; tri_ix < size; ++tri_ix)
		{
			const TV_INT &tri = triangles_[tri_ix];
			const unsigned int x = tri.x_, y = tri.y_, z = tri.z_;

			pFile.write((char*)&(x), sizeof(unsigned int));
			pFile.write((char*)&(y), sizeof(unsigned int));
			pFile.write((char*)&(z), sizeof(unsigned int));
		}

		// vertex normals
// 		size = vertex_normals_.num_elements_;
// 		pFile.write((char*)&size, sizeof(unsigned int));
// 		for (unsigned int v_ix = 0; v_ix < size; ++v_ix)
// 		{
// 			const TV &vn = vertex_normals_[v_ix];
// 
// 			pFile.write((char*)&(vn.x_), sizeof(float));
// 			pFile.write((char*)&(vn.y_), sizeof(float));
// 			pFile.write((char*)&(vn.z_), sizeof(float));
// 		}

		pFile.close();

		cout << "Finished" << endl;
	}
}

void StaticTriangularSurface::writeSTL(const char *filename)
{
	using namespace std;

	cout << "# Writing " << filename << std::endl;
	ofstream pFile(filename, ios::out | ios::binary);

	if (pFile.is_open())
	{
		string name = "Binary STL file";
		pFile.setf(ios::left);
		pFile.width(sizeof(unsigned char)* 80);
		//header
		pFile << name.c_str();
		unsigned int size = triangles_.num_elements_;
		//size
		pFile.write((char*)&size, sizeof(size));

		for (int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
		{
//			Vector3D<T> face_normal = tri_ops_.GetNormal(tri_ix);
			Vector3D<T> face_normal = tri_ops_.getNormalDouble(tri_ix);

			pFile.write((char*)&(face_normal.x_), sizeof(float));
			pFile.write((char*)&(face_normal.y_), sizeof(float));
			pFile.write((char*)&(face_normal.z_), sizeof(float));

			//std::cout << face_normal.x_ << " " << face_normal.y_ << " " << face_normal.z_ << std::endl;

			for (int v_ix = 0; v_ix < 3; v_ix++)
			{
				Vector3D<T> &vertex = vertex_positions_[triangles_[tri_ix][v_ix]];

				pFile.write((char*)&vertex.x_, sizeof(float));
				pFile.write((char*)&vertex.y_, sizeof(float));
				pFile.write((char*)&vertex.z_, sizeof(float));

				//std::cout << vertex.x_ << " "<<vertex.y_ << " "<<vertex.z_ << std::endl;
			}

			const unsigned short attribute_byte_count = 0;

			pFile.write((char*)&attribute_byte_count, sizeof(unsigned short));
		}

		pFile.close();

		cout << "Finished" << endl;
	}
}

void StaticTriangularSurface::writeOBJ(const char *filename)
{
	std::cout << "# Writing " << filename << std::endl;
	std::ofstream file(filename);

//	int frame=0;

	// vertex positions
	for(int v_ix = 0; v_ix < vertex_positions_.num_elements_; v_ix ++)
	{
		file << "v " << vertex_positions_[v_ix].x_ << " " << vertex_positions_[v_ix].y_ << " " << vertex_positions_[v_ix].z_ << std::endl;
	}

	// vertex texture coordinates
	for (int v_ix = 0; v_ix < vertex_positions_.num_elements_; v_ix++)
	{
		file << "vt " << "1" << " " << "0" << std::endl;
	}

	// vertex normals
	for (int v_ix = 0; v_ix < vertex_normals_.num_elements_; v_ix++)
	{
		file << "vn " << vertex_normals_[v_ix].v0_ << " " << vertex_normals_[v_ix].v1_ << " " << vertex_normals_[v_ix].v2_ << std::endl;
	}

	// vertex velocity
// 	for (int v_ix = 0; v_ix < pos_of_vertices_.num_elements_; v_ix++)
// 	{
// 		file << "vv " << vertices_[index - 1]->velocity_.values_[0] * scale << " " << vertices_[index - 1]->velocity_.values_[1] * scale << " " << vertices_[index - 1]->velocity_.values_[2] * scale << std::endl;
// 	}

	for(int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
	{
		const int index[3] = { triangles_[tri_ix].v0_+1, triangles_[tri_ix].v1_+1, triangles_[tri_ix].v2_+1};

		file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << std::endl;
	}

	file.close();
}

void StaticTriangularSurface::backup()
{
	vertex_positions_backup_.initialize(vertex_positions_);
	triangles_backup_.initialize(triangles_);
}

void StaticTriangularSurface::restore()
{
	vertex_positions_.initialize(vertex_positions_backup_);
	triangles_.initialize(triangles_backup_);

	findAdjacentTrianglesOfVertices();
	determineFaceAveragedVertexNormals();
}

void StaticTriangularSurface::replaceVertex(const int& tri_ix, const int& old_v, const int& new_v)
{
	TV_INT &v_ix(triangles_[tri_ix]);

	if (v_ix.v0_ == old_v)
	{
		v_ix.v0_ = new_v;
	}

	if (v_ix.v1_ == old_v)
	{
		v_ix.v1_ = new_v;
	}

	if (v_ix.v2_ == old_v)
	{
		v_ix.v2_ = new_v;
	}
}

void StaticTriangularSurface::emptifyTriangle(const int& tri_ix)
{
	triangles_[tri_ix] = TV_INT(-1, -1, -1);
}

void StaticTriangularSurface::replaceAdjTriangle(const int& v_ix, const int& old_tri_ix, const int& new_tri_ix)
{
	const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
	for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
	{
		int &tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

		if (tri_ix == -1) continue;

		if (tri_ix == old_tri_ix) tri_ix = new_tri_ix;
	}
}

void StaticTriangularSurface::replaceAdjTriangle(const int& v_ix, const int& old_tri_ix1, const int& old_tri_ix2, const int& new_tri_ix)
{
	const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
	for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
	{
		int &tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

		if (tri_ix == -1) continue;

		if (tri_ix == old_tri_ix1 || tri_ix == old_tri_ix2) tri_ix = new_tri_ix;
	}
}

void StaticTriangularSurface::replaceVertexOfAdjTriangles(const int& v_ix, const int& old_v_ix, const int& new_v_ix)
{
	const int t_start = start_ix_adj_tri_of_vertices_[v_ix], t_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
	for (int t = t_start; t < t_end; t++)
	{
		const int v_tri_ix = adj_tri_ix_of_vertices_[t];

		if (v_tri_ix == -1) continue;

//		triangles_flag_[v_tri_ix] = -1;	// tag triangle not to remove in edge collapse

		replaceVertex(v_tri_ix, old_v_ix, new_v_ix);
	}
}

void StaticTriangularSurface::replaceEdgeTriangle(const int& tri_ix, const int& old_tri_ix, const int& new_tri_ix)
{
	TV_INT &edge_tris(edge_tri_ix_of_triangles_[tri_ix]);

	if (edge_tris.v0_ == old_tri_ix) edge_tris.v0_ = new_tri_ix;
	if (edge_tris.v1_ == old_tri_ix) edge_tris.v1_ = new_tri_ix;
	if (edge_tris.v2_ == old_tri_ix) edge_tris.v2_ = new_tri_ix;
}

void StaticTriangularSurface::tagAdjTriangles(const int& v_ix, const int& tri_flag, const int recursive_count)
{
	if (recursive_count > 1) return;

	const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
	for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
	{
		const int v_adj_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

		if (v_adj_tri_ix == -1) continue;

		triangles_flag_[v_adj_tri_ix] = MIN2(triangles_flag_[v_adj_tri_ix], tri_flag);	// tag triangle not to remove in edge collapse
		//TODO: -2 is deleted triangle, overwriting tri_flag=-1 makes problem

		// two neighbor triangles
		tagAdjTriangles(triangles_[v_adj_tri_ix].v0_, tri_flag, recursive_count + 1);
		tagAdjTriangles(triangles_[v_adj_tri_ix].v1_, tri_flag, recursive_count + 1);
		tagAdjTriangles(triangles_[v_adj_tri_ix].v2_, tri_flag, recursive_count + 1);
	}
}

bool StaticTriangularSurface::containInvalidAdjTriangles(const int& v_ix)
{
	const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix + 1];
	for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
	{
		const int v_adj_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

		if (v_adj_tri_ix == -1) continue;

		if (triangles_flag_[v_adj_tri_ix] < 0) return true;	// if this triangle is fixed or deleted elsewhere.
	}

	return false;
}

bool StaticTriangularSurface::checkAllTrianglesContainVIX(const int&v_ix)
{
	for (int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
	{
		for(int d = 0; d < 3; d ++)
		if (triangles_[tri_ix][d] == v_ix)
		{
			std::cout << "triangle " << tri_ix << " " << d << " contains " << v_ix << std::endl;
			return true;
		}
	}

	return false;

}

void StaticTriangularSurface::removeEdge(const int& tri_ix0, const int& tri_ix1, const int& edge_ix0, const int& edge_ix1, const int& v_ix0, const int& v_ix1)
{
	vertices_flag_[v_ix1] = -1;	// deleted

	// merge v0 and v1
	// remove v1 from triangles (replace v1 with v0)
	replaceVertexOfAdjTriangles(v_ix1, v_ix1, v_ix0);

	vertex_positions_[v_ix0] = (vertex_positions_[v_ix0] + vertex_positions_[v_ix1])*0.5f;

	// removing triangles
	removeEdgeSharingTriangles(tri_ix0, tri_ix1, edge_ix0, edge_ix1);
}

void StaticTriangularSurface::removeEdgeSharingTriangles(const int &tri_ix0, const int & tri_ix1, const int& edge_ix0, const int& edge_ix1)
{
	// remove t0 and t1 from adj tri list of all vertices
	for (int d = 0; d < 3; d++) replaceAdjTriangle(triangles_[tri_ix0][d], tri_ix0, -1);
	for (int d = 0; d < 3; d++) replaceAdjTriangle(triangles_[tri_ix1][d], tri_ix1, -1);

	// fix adj tri of tri connectivity
// 	ReplaceEdgeTriangle(edge_tri_ix_of_triangles_[tri_ix0][(edge_ix0 + 1) % 3], tri_ix0, edge_tri_ix_of_triangles_[tri_ix0][(edge_ix0 + 2) % 3]);
// 	ReplaceEdgeTriangle(edge_tri_ix_of_triangles_[tri_ix0][(edge_ix0 + 2) % 3], tri_ix0, edge_tri_ix_of_triangles_[tri_ix1][(edge_ix0 + 1) % 3]);
//   
// 	ReplaceEdgeTriangle(edge_tri_ix_of_triangles_[tri_ix1][(edge_ix1 + 1) % 3], tri_ix1, edge_tri_ix_of_triangles_[tri_ix0][(edge_ix1 + 2) % 3]);
// 	ReplaceEdgeTriangle(edge_tri_ix_of_triangles_[tri_ix1][(edge_ix1 + 2) % 3], tri_ix1, edge_tri_ix_of_triangles_[tri_ix1][(edge_ix1 + 1) % 3]);

	// remove triangles t0 t1
	emptifyTriangle(tri_ix0);
	emptifyTriangle(tri_ix1);

	triangles_flag_[tri_ix0] = -2;	// -2 means removed
	triangles_flag_[tri_ix1] = -2;
}

bool StaticTriangularSurface::checkEdgeConnectivity()
{
	for (int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
	{
		if (triangles_[tri_ix].v0_ == -1) continue;

		if (edge_tri_ix_of_triangles_[tri_ix].v0_ == -1) return false;
		if (edge_tri_ix_of_triangles_[tri_ix].v1_ == -1) return false;
		if (edge_tri_ix_of_triangles_[tri_ix].v2_ == -1) return false;
	}

	return true;
}

int StaticTriangularSurface::removeShortEdges(const T& min_edge_length, const T kappa)
{
	determineVertexMeanCurvatures();

	triangles_flag_.initialize(triangles_.num_elements_, 0);		// -1: fixed (skip to remove this time), -2: removed
	vertices_flag_.initialize(vertex_positions_.num_elements_, 0);	// -1: removed

	const T sqr_length = POW2(min_edge_length);

	int num_removed_tri = 0;

	for (int tri_ix = 0; tri_ix < triangles_.num_elements_; tri_ix++)
	{
		if (triangles_flag_[tri_ix] < 0) continue;		// skip this triangle if this is already deleted or fixed.

		const int short_edge_ix = tri_ops_.getShortEdgeIndex(tri_ix, sqr_length);	// find an edge that is shorter than min_edge_length and shortest among three edges

		if (short_edge_ix == -1) continue;					// this triangle doesn't have shorter edge than sqr_length

		const int tri_ix1 = edge_tri_ix_of_triangles_[tri_ix][short_edge_ix];		// assert(t1 != -1)

		if (tri_ix1 < 0) continue;		// deleted elsewhere
		if (triangles_flag_[tri_ix1] < 0) continue;

		const TV_INT &v_ix(triangles_[tri_ix]);
		const int v_ix0(v_ix[(short_edge_ix + 1) % 3]), v_ix1(v_ix[(short_edge_ix + 2) % 3]);
		const int edge_ix1 = tri_ops_.getEdgeIndex(tri_ix1, v_ix0, v_ix1);

		// check if this edge contains high curvature vertices
		const T kappa_v0 = vertex_mean_curvatures_[v_ix0].getMagnitudeDouble();
		const T kappa_v1 = vertex_mean_curvatures_[v_ix1].getMagnitudeDouble();
 		if (kappa_v0 > (double)kappa) continue;
 		if (kappa_v1 > (double)kappa) continue;

		// check if two-neighbor triangles are fixed or removed
		if (containInvalidAdjTriangles(v_ix[0]) == true) continue;
		if (containInvalidAdjTriangles(v_ix[1]) == true) continue;
		if (containInvalidAdjTriangles(v_ix[2]) == true) continue;
		if (containInvalidAdjTriangles(triangles_[tri_ix1][edge_ix1]) == true) continue;

		const TV new_pos = (vertex_positions_[v_ix0] + vertex_positions_[v_ix1]) * (T)0.5;

		// check triangle flip http://graphics.stanford.edu/courses/cs468-10-fall/LectureSlides/08_Simplification.pdf page 47
		bool flip_expected = false;

		{const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix0], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix0 + 1];
		for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
		{
			const int v_adj_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

			if (v_adj_tri_ix == -1) continue;
// 			if (v_adj_tri_ix == tri_ix) continue;
// 			if (v_adj_tri_ix == tri_ix1) continue;

			if (tri_ops_.checkFlip(v_adj_tri_ix, tri_ops_.getVertexIndex(v_adj_tri_ix, v_ix0), new_pos) == true)
			{
				flip_expected = true;
				break;
			}
		}}

		if (flip_expected == true) continue;

		{const int adj_tri_start = start_ix_adj_tri_of_vertices_[v_ix1], adj_tri_end = start_ix_adj_tri_of_vertices_[v_ix1 + 1];
		for (int adj_tri_ix = adj_tri_start; adj_tri_ix < adj_tri_end; adj_tri_ix++)
		{
			const int v_adj_tri_ix = adj_tri_ix_of_vertices_[adj_tri_ix];

			if (v_adj_tri_ix == -1) continue;
// 			if (v_adj_tri_ix == tri_ix) continue;
// 			if (v_adj_tri_ix == tri_ix1) continue;

			if (tri_ops_.checkFlip(v_adj_tri_ix, tri_ops_.getVertexIndex(v_adj_tri_ix, v_ix1), new_pos) == true)
			{
				flip_expected = true;
				break;
			}
		}}

		if (flip_expected == true) continue;

		// tag these triangles so that they are not deleted this turn.
		tagAdjTriangles(v_ix[0], -1, 0);
		tagAdjTriangles(v_ix[1], -1, 0);
		tagAdjTriangles(v_ix[2], -1, 0);
 		tagAdjTriangles(triangles_[tri_ix1][edge_ix1], -1, 0);

		removeEdge(tri_ix, tri_ix1, short_edge_ix, edge_ix1, v_ix0, v_ix1);

		num_removed_tri++;

// 		// remove only one triangle at a time (test)
// 		break;
	}

	std::cout << num_removed_tri << " triangles are removed " << std::endl;

	if (num_removed_tri == 0) return num_removed_tri;

// 	if (CheckEdgeConnectivity()) std::cout << "Edge connectivity ok" << std::endl;
// 	else std::cout << "Edge connectivity NOT ok " << std::endl;

	triangles_.compactArray(triangles_flag_, -2);
	vertex_positions_.compactArray(vertices_flag_, -1);

	// update old vertex indices to new vertex indices
	for (int tri = 0; tri < triangles_.num_elements_; tri++)
	{
		TV_INT &vertices(triangles_[tri]);

		vertices.v0_ = vertices_flag_[vertices.v0_];
		vertices.v1_ = vertices_flag_[vertices.v1_];
		vertices.v2_ = vertices_flag_[vertices.v2_];
	}

	findAdjacentTrianglesOfVertices();
	determineFaceAveragedVertexNormals();
	findEdgeTrianglesOfTriangles();

	return num_removed_tri;
}

BOX_3D<T> StaticTriangularSurface::getAABB() const
{
	if (vertex_positions_.num_elements_ == 0) return BOX_3D<T>(0, 0, 0, 0, 0, 0);

	TV min(vertex_positions_[0]), max(vertex_positions_[0]);
	
	for (int v_ix = 1; v_ix < vertex_positions_.num_elements_; ++v_ix)
	{
		min.setComponentWiseMin(vertex_positions_[v_ix]);
		max.setComponentWiseMax(vertex_positions_[v_ix]);
	}

	return BOX_3D<T>(min, max);
}
