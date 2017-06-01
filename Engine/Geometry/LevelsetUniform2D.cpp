// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "LevelsetUniform2D.h"

void LevelsetUniform2D::reinitializeFSM(MT* mt, const int& thread_id, const int& itr)
{
	fixInterfacialCells(mt, thread_id);

	for(int t = 0; t < itr; ++t)
	{
		const TV2_INT j_range = mt->getParallelRange(thread_id, grid_.j_start_, grid_.j_end_);
		const int j_start = j_range.t_min_, j_end = j_range.t_max_, i_start = grid_.i_start_, i_end = grid_.i_end_;

		switch(sweep_direction_ % 2)
		{
		case 0:
			sweep(i_start, j_start, i_end, j_end);	mt->sync();
			sweep(i_end, j_end, i_start, j_start);	mt->sync();
			break;

		case 1:
			sweep(i_end, j_start, i_start, j_end);	mt->sync();
			sweep(i_start, j_end, i_end, j_start);	mt->sync();
			break;
		}

		fillGhostCells(mt, thread_id, phi_);

		ONE_THREAD_WORK(++sweep_direction_);
	}
}

/*
void LevelsetUniform2D::reinitializeFSMThreaded(const int& itr)
{
	multithreading_->RunThreads(&LevelsetUniform2D::reinitializeFSM, this, itr);
}
*/
/*
void LevelsetUniform2D::updateInterfacialPhi(const int i, const int j)
{
	static const T h(grid_.dx_), hh(h*h), hh2(hh*(T)2);

	T &phi_center = phi_(i, j);
	T a, b, a1, a2, update;

	if (phi_center == (T)0)		// ABS(phi_center) < 1e-8
	{
		return;	
	}
	if (phi_center > (T)0)
	{
		// See definitions of u^h_{x min}, u^h_{y min} in page 605.
		a = MIN2(phi_ghost_(i - 1, j), phi_ghost_(i + 1, j));
		b = MIN2(phi_ghost_(i, j - 1), phi_ghost_(i, j + 1));

		if (a <= (T)0)
		{
			a = -a / (phi_center - a) * h;
		}

		if (b <= (T)0)
		{
			b = -b / (phi_center - b) * h;
		}

		// Follows the n dimension implementation explained in 606.
		// See "First we order the a_k's in ~".
		INCREASING_SORT2(a, b, a1, a2);

		update = a1 + grid_.dx_;	// dx_ = dy_
		if (update > a2)
		{
			update = (T)0.5*(a1 + a2 + sqrt(hh2 - POW2(a1 - a2)));
		}

		phi_center = MIN2(update, phi_center);
	}
	else // if(phi_center <= (T)0)
	{
		// See Equation (2.8) in page 607 for the implementation of negative parts.
		//					a = MAX2(phi_(i-1,j), phi_(i+1,j));
		//					b = MAX2(phi_(i,j-1), phi_(i,j+1));
		a = MAX2(phi_ghost_(i - 1, j), phi_ghost_(i + 1, j));
		b = MAX2(phi_ghost_(i, j - 1), phi_ghost_(i, j + 1));

		if (a > (T)0)
		{
			a = a / (-phi_center + a) * h;
		}

		if (b > (T)0)
		{
			b = b / (-phi_center + b) * h;
		}

		INCREASING_SORT2(a, b, a1, a2);

		update = a2 - h;
		if (update < a1)
		{
			update = (T)0.5*(a2 + a1 - sqrt(hh2 - POW2(a2 - a1)));
		}

		phi_center = MAX2(update, phi_center);
	}// end of else of if(phi_center > (T)0)
}
*/

void LevelsetUniform2D::reinitializeInterfacialCells(MT* mt, const int& thread_id)
{
	fillGhostCells(mt, thread_id, phi_);

	// copy ghost values
	BEGIN_1D_ITERATION(phi_.getNumAllValues())
	{
		phi_ghost_.values_[p] = phi_.values_[p];	//TODO: memset
	}
	END_1D_ITERATION;

	const int i_res = phi_.i_res_;

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		const int ix = phi_ghost_.get1DIndex(i, j);

		const T phi_ij = phi_ghost_.values_[ix];

		if (phi_ij <= (T)1e-8 && phi_ij >= -(T)1e-8) continue;		// skip so small phi that it approaches to zero

		T dist_x = 1e8, dist_y = 1e8;

		if (phi_ij > (T)0)
		{
			if (phi_ghost_.values_[ix + 1] <= (T)0)	dist_x = MIN2(dist_x, phi_ij / (phi_ij - phi_ghost_.values_[ix + 1]) * grid_.dx_);
			else if (phi_ghost_.values_[ix - 1] <= (T)0) dist_x = MIN2(dist_x, phi_ij / (phi_ij - phi_ghost_.values_[ix - 1]) * grid_.dx_);
			else if (phi_ghost_.values_[ix + i_res] <= (T)0) dist_y = MIN2(dist_y, phi_ij / (phi_ij - phi_ghost_.values_[ix + i_res]) * grid_.dy_);
			else if (phi_ghost_.values_[ix - i_res] <= (T)0) dist_y = MIN2(dist_y, phi_ij / (phi_ij - phi_ghost_.values_[ix - i_res]) * grid_.dy_);

			const T dist = MIN3(dist_x, dist_y, dist_x*dist_y / sqrt(dist_x*dist_x + dist_y*dist_y));

			if (dist < 1e8) phi_.values_[ix] = dist;
		}
		else // if(phi(i,j) <= (T)0)
		{
			if (phi_ghost_.values_[ix + 1] > (T)0) dist_x = MIN2(dist_x, -phi_ij / (-phi_ij + phi_ghost_.values_[ix + 1]) * grid_.dx_);
			else if (phi_ghost_.values_[ix - 1] > (T)0) dist_x = MIN2(dist_x, -phi_ij / (-phi_ij + phi_ghost_.values_[ix - 1]) * grid_.dx_);
			else if (phi_ghost_.values_[ix + i_res] > (T)0) dist_y = MIN2(dist_y, -phi_ij / (-phi_ij + phi_ghost_.values_[ix + i_res]) * grid_.dy_);
			else if (phi_ghost_.values_[ix - i_res] > (T)0) dist_y = MIN2(dist_y, -phi_ij / (-phi_ij + phi_ghost_.values_[ix - i_res]) * grid_.dy_);

			const T dist = MIN3(dist_x, dist_y, dist_x*dist_y / sqrt(dist_x*dist_x + dist_y*dist_y));

			if (dist < 1e8) phi_.values_[ix] = -dist;
		}
	}
	END_GRID_ITERATION_2D;
}

void LevelsetUniform2D::fixInterfacialCells(MT* mt, const int& thread_id)
{
	//NOTE: instead of modify interfacial phi values, move surface front-and-back a little bit

	const T nb_width = (T)10000*grid_.dx_;	//TODO: make this global.

	BEGIN_1D_ITERATION(fixed_.getNumAllValues())
	{
		fixed_.values_[p] = false;		// use memset

		phi_ghost_.values_[p] = phi_.values_[p];
	}
	END_1D_ITERATION;

	int ix;
	int i_res = grid_ghost_.i_res_;

	T* phi_values = phi_.values_;
	bool* fix_values = fixed_.values_;

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		ix = (i - grid_ghost_.i_start_) + (j - grid_ghost_.j_start_) * i_res;

		bool &fix_center(fix_values[ix]);

		if(fix_center == true) continue;

		T& phi_center(phi_values[ix]);

		// fix values of the cells too far from interface
		if(ABS(phi_center) > nb_width)//TODO: check the performance of narrow banding
		{
			fix_center = true;
			continue;
		}

		// fix (do not update) interfacial cells
		if(phi_center > (T)0)
		{
			//TODO: use pointer operations to access neighbor phis.
			if	   (phi_values[ix+1]     <= (T)0) { fix_center = true; fix_values[ix+1]     = true; }
			else if(phi_values[ix-1]     <= (T)0) { fix_center = true; fix_values[ix-1]     = true; }
			else if(phi_values[ix+i_res] <= (T)0) { fix_center = true; fix_values[ix+i_res] = true; }
			else if(phi_values[ix-i_res] <= (T)0) { fix_center = true; fix_values[ix-i_res] = true; }
		}
		else // if(phi(i,j) <= (T)0)
		{
			if     (phi_values[ix+1]     > (T)0) { fix_center = true; fix_values[ix+1]      = true; }
			else if(phi_values[ix-1]     > (T)0) { fix_center = true; fix_values[ix-1]      = true; }
			else if(phi_values[ix+i_res] > (T)0) { fix_center = true; fix_values[ix+i_res]  = true; }
			else if(phi_values[ix-i_res] > (T)0) { fix_center = true; fix_values[ix-i_res]  = true; }
		}

		// assign large magnitude values to non interfacial cells
		if(fix_center == false)
		{
			if(phi_center > (T)0) phi_center = nb_width;
			else phi_center = -nb_width;
		}
	}
	END_GRID_ITERATION_2D;

	fillGhostCells(mt, thread_id, phi_);
}

/*
void LevelsetUniform2D::CuthillFastSweepMethod(const int& thread_id, const int& direction)
{
	const T one_over_three((T)1/(T)3), one_over_two((T)1/(T)2);
	const T h(grid_.dx_), hh(h*h), hh2(hh*(T)2);

	T a, b, a1, a2, update;

	int ip, jp;
	int i_start, j_start, i_end, j_end;

	if(direction == 0)
	{
		i_start = grid_.i_start_;
		i_end = grid_.i_end_;
		j_start = grid_.j_start_;
		j_end = grid_.j_end_;

		jp = 1; ip = 1;
	}

	if(direction == 1)
	{
		i_start = grid_.i_end_;
		i_end = grid_.i_start_;
		j_start = grid_.j_end_;
		j_end = grid_.j_start_;

		jp = -1; ip = -1;
	}

	if(direction == 2)
	{
		i_start = grid_.i_start_;
		i_end = grid_.i_end_;
		j_start = grid_.j_end_;
		j_end = grid_.j_start_;

		jp = -1; ip = 1;
	}

	if(direction == 3)
	{
		i_start = grid_.i_end_;
		i_end = grid_.i_start_;
		j_start = grid_.j_start_;
		j_end = grid_.j_end_;

		jp = 1; ip = -1;
	}

	for(int j_s = j_start; j_s != j_end + jp; j_s += jp)
	{
		int t  = 0;

		while(true)
		{
			int i = i_start + (thread_id + multithreading_->num_threads_ * t)*ip;
			int j = j_s     + (thread_id + multithreading_->num_threads_ * t)*jp;

			++t;
			if(i > grid_.i_end_   || j > grid_.j_end_  ) break;
			if(i < grid_.i_start_ || j < grid_.j_start_) break;

			T &phi_center(phi_(i,j));
			if(fixed_(i,j) == false)
			{
				if(phi_center > (T)0)
				{
					a = MIN(phi_(i-1,j), phi_(i+1,j));
					b = MIN(phi_(i,j-1), phi_(i,j+1));

					INCREASING_SORT2(a, b, a1, a2);

					update = a1 + h;
					if(update > a2)
					{
						update = one_over_two*(a1 + a2 + sqrt(hh2 - POW2(a1-a2)));
					}

					phi_center = MIN(update, phi_center);
				}
				else // if(phi_center <= (T)0)
				{
					a = MAX(phi_(i-1,j), phi_(i+1,j));
					b = MAX(phi_(i,j-1), phi_(i,j+1));

					INCREASING_SORT2(a, b, a1, a2);

					update = a2 - h;
					if(update < a1)
					{
						update = one_over_two*(a2 + a1 - sqrt(hh2 - POW2(a2-a1)));
					}

					phi_center = MAX(update, phi_center);
				}
			}		
		}

		multithreading_->Sync(thread_id);
	}

	for(int i_s = i_start; i_s != i_end + ip; i_s += ip)
	{
		int t  = 0;

		while(true)
		{
			int i = i_s     + (thread_id + multithreading_->num_threads_ * t)*ip;
			int j = j_start + (thread_id + multithreading_->num_threads_ * t)*jp;

			++t;
			if(i > grid_.i_end_   || j > grid_.j_end_  ) break;
			if(i < grid_.i_start_ || j < grid_.j_start_) break;

			T &phi_center(phi_(i,j));
			if(fixed_(i,j) == false)
			{
				if(phi_center > (T)0)
				{
					a = MIN(phi_(i-1,j), phi_(i+1,j));
					b = MIN(phi_(i,j-1), phi_(i,j+1));

					INCREASING_SORT2(a, b, a1, a2);

					update = a1 + h;
					if(update > a2)
					{
						update = one_over_two*(a1 + a2 + sqrt(hh2 - POW2(a1-a2)));
					}

					phi_center = MIN(update, phi_center);
				}
				else // if(phi_center <= (T)0)
				{
					a = MAX(phi_(i-1,j), phi_(i+1,j));
					b = MAX(phi_(i,j-1), phi_(i,j+1));

					INCREASING_SORT2(a, b, a1, a2);

					update = a2 - h;
					if(update < a1)
					{
						update = one_over_two*(a2 + a1 - sqrt(hh2 - POW2(a2-a1)));
					}

					phi_center = MAX(update, phi_center);
				}
			}		
		}

		multithreading_->Sync(thread_id);
	}

	multithreading_->Sync(thread_id);
}
*/

void LevelsetUniform2D::sweep(const int& i_start, const int& j_start, const int& i_end, const int& j_end)
{
	const T one_over_three((T)1/(T)3), one_over_two((T)1/(T)2);
	const T h(grid_.dx_), hh(h*h), hh2(hh*(T)2);

	T a, b;
	T a1, a2;
	T update;

	// traverse orders
	int ip(1), jp(1);        			// increasing order
	if(i_start > i_end) ip = -1;		// decreasing order
	if(j_start > j_end) jp = -1;		// decreasing order

	int i, j, ix;
	int i_res = grid_ghost_.i_res_;

	T* phi_values = phi_.values_;
	bool* fix_values = fixed_.values_;

	j = j_start;
	while(j != j_end + jp)
	{
		i = i_start;
		ix = (i - grid_ghost_.i_start_) + (j - grid_ghost_.j_start_) * i_res;
		while(i != i_end + ip)
		{
			T &phi_center(phi_values[ix]);

			if(fix_values[ix] == false)
			{
				if(phi_center > (T)0)
				{
					// See definitions of u^h_{x min}, u^h_{y min} in page 605.
					a = MIN2(phi_values[ix-1], phi_values[ix+1]);
					b = MIN2(phi_values[ix-i_res], phi_values[ix+i_res]);

					// Follows the n dimension implementation explained in 606.
					// See "First we order the a_k's in ~".
					INCREASING_SORT2(a, b, a1, a2);

					update = a1 + h;
					if(update > a2)
					{
						update = one_over_two*(a1 + a2 + sqrt(hh2 - POW2(a1-a2)));
					}

					phi_center = MIN2(update, phi_center);
				}
				else // if(phi_center <= (T)0)
				{
					// See Equation (2.8) in page 607 for the implementation of negative parts.
//					a = MAX2(phi_(i-1,j), phi_(i+1,j));
//					b = MAX2(phi_(i,j-1), phi_(i,j+1));
					a = MAX2(phi_values[ix-1], phi_values[ix+1]);
					b = MAX2(phi_values[ix-i_res], phi_values[ix+i_res]);

					INCREASING_SORT2(a, b, a1, a2);

					update = a2 - h;
					if(update < a1)
					{
						update = one_over_two*(a2 + a1 - sqrt(hh2 - POW2(a2-a1)));
					}

					phi_center = MAX2(update, phi_center);
				}// end of else of if(phi_center > (T)0)
			}

			i += ip;
			ix += ip; 
		}
		j += jp;
	}
}// end of Sweep

/*
template <class TT> 
void LevelsetUniform2D::ExtrapolatingSweep(FIELD_UNIFORM_2D<TT>& value_ghost, LevelsetUniform2D& object_levelset, const int& i_start, const int& j_start, const int& i_end, const int& j_end)
{
	// traverse orders
	int ip(1), jp(1);					// increasing order
	if(i_start > i_end) ip = -1;		// decreasing order
	if(j_start > j_end) jp = -1;		// decreasing order

	const T half_dx = (T)0.5*grid_.dx_;

	TT weighted_val_sum, v_min;
	T weight_sum, phi_min;
	T phi_diff;							// spatial partial derivative
	int i, j;

	j = j_start - jp;
	while(true)
	{
		if(j == j_end) break;
		else j += jp;

		i = i_start - ip;
		while(true)
		{
			if(i == i_end) break;
			else i += ip;

			const T &phi_center(phi_(i,j,k));

			if(phi_center <=(T)0 && object_levelset(i,j,k) > (T)0) continue;			// extrapolate values at positive levelset cells and object cells
//			if(phi_center < -half_dx && object_levelset(i,j,k) > (T)0) continue;		// extrapolate values at positive levelset cells and object cells

			weighted_val_sum = TT();
			weight_sum = (T)0;

			for(int d = 0; d < 3; d ++)
			{
				TV2_INT ixp(i,j,k), ixm(i,j,k);
				ixp.values_[d] ++; ixm.values_[d] --;		//TODO: use 1D array index number

				if(phi_(ixp) < phi_(ixm))
				{
					phi_min = phi_(ixp);
					v_min = value_ghost(ixp);
				}
				else
				{
					phi_min = phi_(ixm);
					v_min = value_ghost(ixm);
				}

				phi_diff = phi_center - phi_min;

				if(phi_diff >= (T)0)
				{
					weighted_val_sum += v_min*phi_diff;
					weight_sum += phi_diff;
				}
			}

			if(weight_sum != 0) value_ghost(i,j,k) = weighted_val_sum / weight_sum;
			//Note: "else value_ghost(i,j,k) = TT()" is unnecessary. If you would like to erase the velocities outside narraw band, it is more efficient erase them before starting sweepings.
			//TODO: make a 'clearbeforesweeping' option.
		}// end i
	}// end j

}// end of Sweep
*/

/*
template<class TT> 
void LevelsetUniform2D::FastSweepingExtrapolation(const int& thread_id, FIELD_UNIFORM_2D<TT>& value, FIELD_UNIFORM_2D<TT>& value_ghost, LevelsetUniform2D& object_levelset)
{	
	value_ghost.FillGhostCellsFrom(thread_id, value.array_, true);

	INIT_GRIDRANGE_2D(value.partial_grids_[thread_id], i_start, j_start, i_end, j_end);

	static int extrapolating_sweep_direction = 0;// to sweep 2 directions in tern.

	switch(extrapolating_sweep_direction % 4)	
	{
	case 0:
		ExtrapolatingSweep(value_ghost, object_levelset, i_start, j_start, i_end, j_end); multithreading_->Sync(thread_id);
		ExtrapolatingSweep(value_ghost, object_levelset, i_end, j_end, i_start, j_start); multithreading_->Sync(thread_id);
		break;
	case 1:
		ExtrapolatingSweep(value_ghost, object_levelset, i_end, j_start, i_start, j_end); multithreading_->Sync(thread_id);
		ExtrapolatingSweep(value_ghost, object_levelset, i_start, j_end, i_end, j_start); multithreading_->Sync(thread_id);
		break;
	case 2:
		ExtrapolatingSweep(value_ghost, object_levelset, i_start, j_end, i_end, j_start); multithreading_->Sync(thread_id);
		ExtrapolatingSweep(value_ghost, object_levelset, i_end, j_start, i_start, j_end); multithreading_->Sync(thread_id);
		break;
	case 3:
		ExtrapolatingSweep(value_ghost, object_levelset, i_end, j_end, i_start, j_start); multithreading_->Sync(thread_id);
		ExtrapolatingSweep(value_ghost, object_levelset, i_start, j_start, i_end, j_end); multithreading_->Sync(thread_id);
		break;
	}

	if(thread_id == 0) extrapolating_sweep_direction ++;

	value.CopyAllValuesFrom(thread_id, value_ghost);
}
*/

void LevelsetUniform2D::computeNormals(MT* mt, const int& thread_id)
{
	// fill ghost cells of phi
	BEGIN_GRID_ITERATION_2D(grid_ghost_)
	{
		phi_(i, j) = phi_(grid_.ClampedIndex(i, j));
	}
	END_GRID_ITERATION_2D;

	// calculate normal vectors
	BEGIN_GRID_ITERATION_2D(grid_)
	{
		TV2 &nor(normal_(i, j));

		nor.x_ = (phi_(i + 1, j) - phi_(i - 1, j))*grid_.one_over_2dx_;
		nor.y_ = (phi_(i, j + 1) - phi_(i, j - 1))*grid_.one_over_2dy_;

		nor.normalize();
	}
	END_GRID_ITERATION_2D;

	// fill ghost cells of normal
	BEGIN_GRID_ITERATION_2D(grid_ghost_)
	{
		normal_(i, j) = normal_(grid_.ClampedIndex(i, j));
	}
	END_GRID_ITERATION_2D;
}

void LevelsetUniform2D::computeCurvatures(MT* mt, const int& thread_id)
{
//	normal_.FillGhostCellsFrom(thread_id, normal_.array_, false);

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		T curv((normal_(i + 1, j).x_ - normal_(i - 1, j).x_)*grid_.one_over_2dx_);
		curv += (normal_(i, j + 1).y_ - normal_(i, j - 1).y_)*grid_.one_over_2dy_;		

		curvature_(i, j) = curv;
	}
	END_GRID_ITERATION_2D;

	fillGhostCells(mt, thread_id, curvature_);
}

T LevelsetUniform2D::computeCurvature(const int& i, const int& j)
{
	return (normal_(i + 1, j).x_ - normal_(i - 1, j).x_)*grid_ghost_.one_over_2dx_ + (normal_(i, j + 1).y_ - normal_(i, j - 1).y_)*grid_ghost_.one_over_2dy_;
}

void LevelsetUniform2D::computeNormalThreaded(MT* mt)
{
	mt->runWithID(&LevelsetUniform2D::computeNormals, this);
}