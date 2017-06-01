// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "GridUniform2D.h"
#include "Array1D.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// void GRID_UNIFORM_2D::SplitInHeight(const int& num_threads, Array1D<GRID_UNIFORM_2D>& partial_grids)
// {
// 	partial_grids.Initialize(num_threads);
// 
// 	int quotient = j_res_ / num_threads;
// 	int remainder = j_res_ % num_threads;
// 	int j_start_p = j_start_;
// 	T y_min_p = y_min_;
// 
// 	for(int i=0; i<num_threads; i++)
// 	{
// 		int y_height = i < remainder ? quotient+1 : quotient;
// 
// 		T y_max_p = dy_ * y_height + y_min_p;
// 
// 		partial_grids[i].Initialize(i_start_, j_start_p, i_res_, y_height, x_min_, y_min_p, x_max_, y_max_p);
// 
// 		j_start_p += y_height;
// 		y_min_p = y_max_p;
// 	}
// }
// 


