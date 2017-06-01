// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <list>
#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"
#include "DataStructure/LinkedArray.h"
#include "Geometry/Quadrilateral2D.h"

class DynamicContour2D
{
public:

	T max_dl_, min_dl_;

	Array1D<TV2>		vertex_positions_;
	Array1D<TV2>		vertex_velocities_;
	Array1D<TV2_INT>	lines_;

	// vertex-lines connectivity
	Array1D<int> start_ix_adj_lines_of_vertices_;		// 1. count number of adjacent neighbors, 2. store index of neighboring lines index
	Array1D<int> adj_line_ix_of_vertices_;				// stores the indices of connected lines	//TODO: atomic

	DynamicContour2D()
		: max_dl_(1e8), min_dl_(0.0f)
	{}

	~DynamicContour2D()
	{
		reset();
	}

	Quadrilateral2D getAdvanceQuad(const int line_ix, const T& dt) const
	{
		const TV2 v0 = vertex_positions_[lines_[line_ix].v0_], v1 = vertex_positions_[lines_[line_ix].v1_];
		const TV2 v0_new = v0 + vertex_velocities_[lines_[line_ix].v0_] * dt, v1_new = v1 + vertex_velocities_[lines_[line_ix].v1_] * dt;

		return Quadrilateral2D(v1, v0, v0_new, v1_new);
	}

	void reset()
	{
		vertex_positions_.freeMemory();
		lines_.freeMemory();
	}

	T getLineLength(const int line_ix) const
	{
		return (vertex_positions_[lines_[line_ix].v0_] - vertex_positions_[lines_[line_ix].v1_]).getMagnitude();
	}

	TV2& getP0(const int& line_ix)
	{
		return vertex_positions_[lines_[line_ix].x_];
	}

	TV2& getP1(const int& line_ix)
	{
		return vertex_positions_[lines_[line_ix].y_];
	}

	// split line if it is longer than max_dl
	int splitHelper(const T max_dl)
	{
// 		for (int l = 0; l < lines_.num_elements_; l++)
// 		{
// 			std::cout << lines_[l].x_ << " " << lines_[l].y_ << std::endl;
// 		}
// 		std::cout<<" ---------------------"<<std::endl;

		LinkedArray<TV2> new_pts;
		LinkedArray<TV2_INT> new_lines;

		const int temp = vertex_positions_.num_elements_;
		int count = vertex_positions_.num_elements_;

		for (int l = 0; l < lines_.num_elements_; l++)
		{
			if (getLineLength(l) > max_dl)
			{
				// add new vertex
				new_pts.PushBack() = (getP0(l) + getP1(l)) * (T)0.5;

				// splint this line to 2 sub lines
				const  int p1_temp = lines_[l].y_;
				lines_[l].y_ = count;

				new_lines.PushBack() = TV2_INT(count, p1_temp);

				count++;
			}
		}

		vertex_positions_.resize(count);
		lines_.resize(count);

		new_pts.CopyToPartialArray(vertex_positions_, temp);
		new_lines.CopyToPartialArray(lines_, temp);

		return count - temp;

// 		int num_deleted = 0;
// 		for (std::list<DynamicLine*>::iterator itr = lines_.begin(); itr != lines_.end(); itr++)
// 		{
// 			if ((*itr)->getLength() > max_dl)
// 			{
// 				splitLine(*itr);
// 				SAFE_DELETE(*itr);
// 				num_deleted++;
// 			}
// 		}
// 
// 		lines_.remove(nullptr);
// 
// 		return num_deleted;
	}

	void split()
	{
		while (true)
		{
			if (splitHelper(max_dl_) == 0) break;
		}
	}

	void findAdjacentLinesOfVertices()
	{
		const int num_vertices = vertex_positions_.num_elements_;
		const int num_lines = lines_.num_elements_;

		start_ix_adj_lines_of_vertices_.initialize(num_vertices + 1, 0);		// last value will be used to find the range of adjacent triangles of the last vertex

		// count the number of triangles connected to each vertex and store the number in ver_index_ Array1D
		for (int line_ix = 0; line_ix < num_lines; line_ix++)
		{
			const TV2_INT &vertex_indices(lines_[line_ix]);

			start_ix_adj_lines_of_vertices_[vertex_indices.v0_] ++;
			start_ix_adj_lines_of_vertices_[vertex_indices.v1_] ++;
		}

		// find the index of 1D adjacent triangle indices
		int accumulated_index = 0;
		for (int ver_ix = 0; ver_ix < num_vertices; ver_ix++)
		{
			const int temp = start_ix_adj_lines_of_vertices_[ver_ix];

			start_ix_adj_lines_of_vertices_[ver_ix] = accumulated_index;
			accumulated_index += temp;
		}

		start_ix_adj_lines_of_vertices_[num_vertices] = accumulated_index;	// last value helps the adj triangles iteration of last vertex

		adj_line_ix_of_vertices_.initialize(accumulated_index);	// this is the number of all connected triangles

		Array1D<int> line_count(start_ix_adj_lines_of_vertices_);					// this is a counter TODO: atomic

		// store the indices of adjacent triangles
		for (int line_ix = 0; line_ix < num_lines; line_ix++)
		{
			const TV2_INT &vertex_indices(lines_[line_ix]);

			adj_line_ix_of_vertices_[line_count[vertex_indices.v0_]++] = line_ix;
			adj_line_ix_of_vertices_[line_count[vertex_indices.v1_]++] = line_ix;
		}

		// how to iterate
		for (int v_ix = 0; v_ix < vertex_positions_.num_elements_; v_ix++)
		{
			std::cout << "Lines connected to vertex " << v_ix << " are lines # ";

			const int adj_line_start = start_ix_adj_lines_of_vertices_[v_ix], adj_line_end = start_ix_adj_lines_of_vertices_[v_ix + 1];
			for (int adj_line_ix = adj_line_start; adj_line_ix < adj_line_end; adj_line_ix++)
			{
				const int v_line_ix = adj_line_ix_of_vertices_[adj_line_ix];
				
				std::cout << v_line_ix << " ";
			}

			std::cout << std::endl;
		}
	}

	void removeShortLines()
	{
		while (removeShortLinesHelper() != 0)
		{}
	}

	int removeShortLinesHelper()
	{
		Array1D<int> vertex_flags_;
		vertex_flags_.initialize(vertex_positions_.num_elements_, -1);		// flag -1 means no-change.

		// iterate lines and find short lines to remove.
		int line_remove = 0;
		for (int line_ix = 0; line_ix < lines_.num_elements_; line_ix++)
		{
			if (getLineLength(line_ix) < min_dl_)
			{
				const int v0ix = lines_[line_ix].v0_, v1ix = lines_[line_ix].v1_;

				if (vertex_flags_[v0ix] < 0 || vertex_flags_[v1ix] < 0) continue;	// do not delete lines containing deleted or moved particles

				vertex_positions_[v0ix] = (vertex_positions_[v0ix] + vertex_positions_[v1ix]) * (T)0.5;	// move v0 to middle point

				vertex_flags_[v0ix] = -2;	// tag moved partciel
				vertex_flags_[v1ix] = v0ix;	// let other lines containing v1 replace v1 to v0

				lines_[line_ix].v0_ = -1;	// invalidate this line
				lines_[line_ix].v1_ = -1;

				line_remove++;
			}
		}

//		std::cout<<"# of lines to be removed "<< line_remove <<std::endl;

		{
			// compact line list
			int count = 0;
//			TV2_INT *line_temp = new TV2_INT[lines_.num_elements_ - line_remove];
			Array1D<TV2_INT> line_temp(lines_.num_elements_ - line_remove);

			for (int line_ix = 0; line_ix < lines_.num_elements_; line_ix++)
			{
				if (lines_[line_ix].v0_ == -1) continue;	// skip deleted lines

				if (vertex_flags_[lines_[line_ix].v0_] >= 0) lines_[line_ix].v0_ = vertex_flags_[lines_[line_ix].v0_];		// replace deleted points
				if (vertex_flags_[lines_[line_ix].v1_] >= 0) lines_[line_ix].v1_ = vertex_flags_[lines_[line_ix].v1_];

				line_temp[count ++] = lines_[line_ix];	// copy lines
			}

			SWAP(lines_.values_, line_temp.values_, TV2_INT*);
			
			assert(count == lines_.num_elements_ - line_remove);

			lines_.num_elements_ = lines_.num_elements_ - line_remove;

//			SAFE_DELETE_ARRAY(line_temp);
		}

//		return;

		{
			// compact vertex list
			// swap flag
			int count = 0;	// count the number of non-deleted vertices
			for (int i = 0; i < vertex_flags_.num_elements_; i++)
			{
				if (vertex_flags_[i] != -1) vertex_flags_[i] = -1;		// this non -1 flag is not necessary any more because line v ixes are updated already.
				else vertex_flags_[i] = count++;
			}

//			std::cout << "# of deleted vertices " << vertex_flags_.num_elements_ - count << std::endl;

			TV2 *temp_pos = new TV2[count];

			// copy vertices to new array
			for (int i = 0; i < vertex_flags_.num_elements_; i++)
			{
				if (vertex_flags_[i] == -1) continue;

				temp_pos[vertex_flags_[i]] = vertex_positions_[i];

				assert(vertex_flags_[i] != -1);
			}

			SWAP(vertex_positions_.values_, temp_pos, TV2*);
			vertex_positions_.num_elements_ = count;

			SAFE_DELETE_ARRAY(temp_pos);
		}

		{
			// update line v ixes
			for (int line_ix = 0; line_ix < lines_.num_elements_; line_ix++)
			{
				TV2_INT &line = lines_[line_ix];

				const int new_v0_ix = vertex_flags_[line.v0_];
				const int new_v1_ix = vertex_flags_[line.v1_];

				assert(new_v0_ix != -1);
				assert(new_v1_ix != -1);

				line.v0_ = new_v0_ix;
				line.v1_ = new_v1_ix;
			}
		}

//		findAdjacentLinesOfVertices();

		return line_remove;
	}

	int removeIntersectingLines(const T dt)
	{
		Array1D<int> vertex_flags_;
		vertex_flags_.initialize(vertex_positions_.num_elements_, -1);		// flag -1 means no-change.

		// iterate lines and find short lines to remove.
		int line_remove = 0;
		for (int line_ix = 0; line_ix < lines_.num_elements_; line_ix++)
		{
			if (isIntersecting(line_ix, dt) == true)
			{
				const int v0ix = lines_[line_ix].v0_, v1ix = lines_[line_ix].v1_;

				if (vertex_flags_[v0ix] < 0 || vertex_flags_[v1ix] < 0) continue;	// do not delete lines containing deleted or moved particles

				vertex_positions_[v0ix] = (vertex_positions_[v0ix] + vertex_positions_[v1ix]) * (T)0.5;	// move v0 to middle point

				vertex_flags_[v0ix] = -2;	// tag moved partciel
				vertex_flags_[v1ix] = v0ix;	// let other lines containing v1 replace v1 to v0

				lines_[line_ix].v0_ = -1;	// invalidate this line
				lines_[line_ix].v1_ = -1;

				line_remove++;
			}
		}

		//		std::cout<<"# of lines to be removed "<< line_remove <<std::endl;

		{
			// compact line list
			int count = 0;
			//			TV2_INT *line_temp = new TV2_INT[lines_.num_elements_ - line_remove];
			Array1D<TV2_INT> line_temp(lines_.num_elements_ - line_remove);

			for (int line_ix = 0; line_ix < lines_.num_elements_; line_ix++)
			{
				if (lines_[line_ix].v0_ == -1) continue;	// skip deleted lines

				if (vertex_flags_[lines_[line_ix].v0_] >= 0) lines_[line_ix].v0_ = vertex_flags_[lines_[line_ix].v0_];		// replace deleted points
				if (vertex_flags_[lines_[line_ix].v1_] >= 0) lines_[line_ix].v1_ = vertex_flags_[lines_[line_ix].v1_];

				line_temp[count++] = lines_[line_ix];	// copy lines
			}

			SWAP(lines_.values_, line_temp.values_, TV2_INT*);

			assert(count == lines_.num_elements_ - line_remove);

			lines_.num_elements_ = lines_.num_elements_ - line_remove;

			//			SAFE_DELETE_ARRAY(line_temp);
		}

		//		return;

		{
			// compact vertex list
			// swap flag
			int count = 0;	// count the number of non-deleted vertices
			for (int i = 0; i < vertex_flags_.num_elements_; i++)
			{
				if (vertex_flags_[i] != -1) vertex_flags_[i] = -1;		// this non -1 flag is not necessary any more because line v ixes are updated already.
				else vertex_flags_[i] = count++;
			}

			//			std::cout << "# of deleted vertices " << vertex_flags_.num_elements_ - count << std::endl;

			TV2 *temp_pos = new TV2[count];

			// copy vertices to new array
			for (int i = 0; i < vertex_flags_.num_elements_; i++)
			{
				if (vertex_flags_[i] == -1) continue;

				temp_pos[vertex_flags_[i]] = vertex_positions_[i];

				assert(vertex_flags_[i] != -1);
			}

			SWAP(vertex_positions_.values_, temp_pos, TV2*);
			vertex_positions_.num_elements_ = count;

			SAFE_DELETE_ARRAY(temp_pos);
		}

		{
			// update line v ixes
			for (int line_ix = 0; line_ix < lines_.num_elements_; line_ix++)
			{
				TV2_INT &line = lines_[line_ix];

				const int new_v0_ix = vertex_flags_[line.v0_];
				const int new_v1_ix = vertex_flags_[line.v1_];

				assert(new_v0_ix != -1);
				assert(new_v1_ix != -1);

				line.v0_ = new_v0_ix;
				line.v1_ = new_v1_ix;
			}
		}

		//		findAdjacentLinesOfVertices();

		return line_remove;
	}

	bool isIntersecting(const int line_ix, const T dt)
	{
		const TV2 v0 = vertex_positions_[lines_[line_ix].v0_], v1 = vertex_positions_[lines_[line_ix].v1_];
		const TV2 d0 = vertex_velocities_[lines_[line_ix].v0_] * dt, d1 = vertex_velocities_[lines_[line_ix].v1_] * dt;

		const TV2 l = (v1 - v0).getNormalized();
		const T length = (v1 - v0).getMagnitude();

		if (dotProduct(d0, l) - dotProduct(d1, l) >= length) return true;
		return false;
	}
};