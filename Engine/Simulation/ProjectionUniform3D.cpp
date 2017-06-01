//#include "SYSTEM_MANAGER_3D.h"
#include "Geometry/BOX_3D.h"
#include "Parallelism/MultiThreading.h"
#include "ProjectionUniform3D.h"
//#include "VERTEX.h"
//#include "DEFORMABLE_OBJECT_3D.h"

//#define R2B(rigid_body_index) (-rigid_body_index-3)		// rigid body index to boundary index
//#define B2R(boundary_index)   (-boundary_index-3)		// boundary index to rigid body index

ProjectionUniform3D::ProjectionUniform3D(void)
{
}

ProjectionUniform3D::~ProjectionUniform3D(void)
{
}

/*
void ProjectionUniform3D::Update(const int& thread_id, const PROJECTION_MODE& mode, SIM_DATA_UNIFORM_3D& data_)
{
	DetermineDivergence(thread_id, mode, data_);
//	DetermineDivergence2ndOrder(thread_id, data_.grid_, data_.boundary_array_, data_.velocity_array_, data_.rigid_body_list_, data_.divergence_array_);
//	DeterminePressureMatrixFreeGaussSeidel(thread_id, data_.grid_, data_.pressure_array_, data_.boundary_array_, data_.divergence_array_);
	DeterminePressureCPUConjugateGradientMethod(thread_id, data_.domain_, data_.boundary_array_, data_.divergence_array_, data_.pressure_array_);
	UpdateVelocityByPressureGradient(thread_id, data_);
}
*/

/*
void ProjectionUniform3D::UpdateVelocityByPressureGradient(const int& thread_id, SIM_DATA_UNIFORM_3D& data_)
{
	SIM_DOMAIN_UNIFORM_3D &grid_(data_.domain_);
	T *pressure_array_ = data_.pressure_array_;
	int *boundary_array_ = data_.boundary_array_;
	TV3 *velocity_array_ = data_.velocity_array_;
	const int i_start = grid_.i_start_, i_end = grid_.i_end_;
	const int j_start = grid_.j_start_, j_end = grid_.j_end_;
	const int cell_num_i = grid_.cell_num_i_;

	const int i_res = grid_.cell_num_i_;
	const int ij_res = grid_.cell_num_i_*grid_.cell_num_j_;

	// update velocity
	T p_ijk, p_ip, p_im, p_jp, p_jm, p_kp, p_km;
	BEGIN_PARALLEL_FOR(k, grid_.k_start_, grid_.k_end_)
	{
		for(int j=j_start;j<=j_end;++j) for(int i=i_start,ix=grid_.Get1DIndexCell(i,j,k);i<=i_end;++i,++ix)
		{
			p_ijk = pressure_array_[ix];

			if(boundary_array_[ix+1] >= 0)						 p_ip = pressure_array_[ix+1];
			else if(boundary_array_[ix+1] <= BOUNDARY_NEUMANN)	 p_ip =  p_ijk;
			else if(boundary_array_[ix+1] == BOUNDARY_DIRICHLET) p_ip = -p_ijk;
			else {assert(false);}

			if(boundary_array_[ix-1] >= 0)						 p_im = pressure_array_[ix-1];
			else if(boundary_array_[ix-1] <= BOUNDARY_NEUMANN)   p_im =  p_ijk;
			else if(boundary_array_[ix-1] == BOUNDARY_DIRICHLET) p_im = -p_ijk;
			else {assert(false);}


			if(boundary_array_[ix+i_res] >= 0)						 p_jp = pressure_array_[ix+i_res];
			else if(boundary_array_[ix+i_res] <= BOUNDARY_NEUMANN)	 p_jp =  p_ijk;
			else if(boundary_array_[ix+i_res] == BOUNDARY_DIRICHLET) p_jp = -p_ijk;
			else {assert(false);}

			if(boundary_array_[ix-i_res] >= 0)						 p_jm = pressure_array_[ix-i_res];
			else if(boundary_array_[ix-i_res] <= BOUNDARY_NEUMANN)	 p_jm =  p_ijk;
			else if(boundary_array_[ix-i_res] == BOUNDARY_DIRICHLET) p_jm = -p_ijk;
			else {assert(false);}


			if(boundary_array_[ix+ij_res] >= 0)						  p_kp = pressure_array_[ix+ij_res];
			else if(boundary_array_[ix+ij_res] <= BOUNDARY_NEUMANN)	  p_kp =  p_ijk;
			else if(boundary_array_[ix+ij_res] == BOUNDARY_DIRICHLET) p_kp = -p_ijk;
			else {assert(false);}

			if(boundary_array_[ix-ij_res] >= 0)						  p_km = pressure_array_[ix-ij_res];
			else if(boundary_array_[ix-ij_res] <= BOUNDARY_NEUMANN)	  p_km =  p_ijk;
			else if(boundary_array_[ix-ij_res] == BOUNDARY_DIRICHLET) p_km = -p_ijk;
			else {assert(false);}

			TV3 acc = TV3((T)0.5*(p_ip-p_im), (T)0.5*(p_jp-p_jm), (T)0.5*(p_kp-p_km));

			velocity_array_[ix] -= acc;
//			data_.mpm_force_array_[ix] = TV2((T)-0.5*(p_ip-p_im), (T)-0.5*(p_jp-p_jm));			
		}
	}
	END_PARALLEL_FOR;
}
*/

/*
void ProjectionUniform3D::DetermineDivergence(const int& thread_id, const PROJECTION_MODE& mode_, SIM_DATA_UNIFORM_3D& data_)
{
	const SIM_DOMAIN_UNIFORM_3D &grid_(data_.domain_);

	int *boundary_array_	= data_.boundary_array_;
	TV3	*velocity_array_	= data_.velocity_array_;
	T	*divergence_array_	= data_.divergence_array_;
	T   *expansion_array_   = data_.expansion_array_;
	T	*mpm_density_array_ = data_.mpm_density_array_;

	const T coef = (T)1;
	const T rigid_body_coef = (T)0.1;
	
	const T rho0 = params_.rho0_;	//TODO: this value need to be considered when particles are seeded at first.
	const T inv_rho0 = (T)1/rho0;
	const T pt_rho_coef = mode_ == PROJECT_WATER ? params_.rho0_div_control_coef_ : 0;		// volume control with respect to the MPM density. 

	// calculate divergence
	T u_ijk, v_ijk, w_ijk;
	T u_p, u_m, v_p, v_m, w_p, w_m;

	int i_res = grid_.cell_num_i_;
	int ij_res = grid_.cell_num_i_*grid_.cell_num_j_;

	BEGIN_PARALLEL_FOR(k, grid_.k_start_, grid_.k_end_)
	{
		for(int j=grid_.j_start_; j<=grid_.j_end_; j++) for(int i=grid_.i_start_,ix=grid_.Get1DIndexCell(i,j,k);i<=grid_.i_end_;++i,++ix)
		{
			if(boundary_array_[ix] < 0) continue;

			u_ijk = velocity_array_[ix].x_;
			v_ijk = velocity_array_[ix].y_;
			w_ijk = velocity_array_[ix].z_;

			if     (boundary_array_[ix+1] >= 0)                  u_p = velocity_array_[ix+1].x_;
			else if(boundary_array_[ix+1] <= BOUNDARY_NEUMANN)   u_p = coef*(velocity_array_[ix+1].x_-u_ijk)+u_ijk;
			else if(boundary_array_[ix+1] <= BOUNDARY_DIRICHLET) u_p = u_ijk;
			else {assert(false);}

			if     (boundary_array_[ix-1] >= 0)                  u_m = velocity_array_[ix-1].x_;
			else if(boundary_array_[ix-1] <= BOUNDARY_NEUMANN)   u_m = coef*(velocity_array_[ix-1].x_-u_ijk)+u_ijk;
			else if(boundary_array_[ix-1] <= BOUNDARY_DIRICHLET) u_m = u_ijk;
			else {assert(false);}


			if     (boundary_array_[ix+i_res] >= 0)                  v_p = velocity_array_[ix+i_res].y_;
			else if(boundary_array_[ix+i_res] <= BOUNDARY_NEUMANN)   v_p = coef*(velocity_array_[ix+i_res].y_-v_ijk)+v_ijk;
			else if(boundary_array_[ix+i_res] <= BOUNDARY_DIRICHLET) v_p = v_ijk;
			else {assert(false);}

			if     (boundary_array_[ix-i_res] >= 0)                  v_m = velocity_array_[ix-i_res].y_;
			else if(boundary_array_[ix-i_res] <= BOUNDARY_NEUMANN)   v_m = coef*(velocity_array_[ix-i_res].y_-v_ijk)+v_ijk;
			else if(boundary_array_[ix-i_res] <= BOUNDARY_DIRICHLET) v_m = v_ijk;
			else {assert(false);}


			if     (boundary_array_[ix+ij_res] >= 0)                  w_p = velocity_array_[ix+ij_res].z_;
			else if(boundary_array_[ix+ij_res] <= BOUNDARY_NEUMANN)   w_p = coef*(velocity_array_[ix+ij_res].z_-w_ijk)+w_ijk;
			else if(boundary_array_[ix+ij_res] <= BOUNDARY_DIRICHLET) w_p = w_ijk;
			else {assert(false);}

			if     (boundary_array_[ix-ij_res] >= 0)                  w_m = velocity_array_[ix-ij_res].z_;
			else if(boundary_array_[ix-ij_res] <= BOUNDARY_NEUMANN)   w_m = coef*(velocity_array_[ix-ij_res].z_-w_ijk)+w_ijk;
			else if(boundary_array_[ix-ij_res] <= BOUNDARY_DIRICHLET) w_m = w_ijk;
			else {assert(false);}

			divergence_array_[ix] = (T)0.5*(u_p-u_m+v_p-v_m+w_p-w_m) + expansion_array_[ix] - pt_rho_coef*MAX((T)0, (mpm_density_array_[ix]*inv_rho0-(T)1));
//			divergence_array_[ix] = (T)0.5*(u_p-u_m+v_p-v_m+w_p-w_m) - pt_rho_coef*MAX((T)0, (mpm_density_array_[ix]*inv_rho0 - (T)1));
//			divergence_array_[ix] = (T)0.5*(u_p-u_m+v_p-v_m+w_p-w_m) - pt_rho_coef*(mpm_density_array_[ix]*inv_rho0 - (T)1);
//			divergence_array_[ix] = (T)0.5*(u_p-u_m+v_p-v_m+w_p-w_m);
		}
	};
	END_PARALLEL_FOR;
}
*/

void ProjectionUniform3D::InitializeLinearSystem(MT* mt, const int& thread_id, const GridUniform3D& grid_g, int* boundary_array, CSRMatrix<T>& A_matrix, VectorND<T>& x_vector, VectorND<T>& b_vector)
{
	int full_ix(0);
	BEGIN_PARALLEL_FOR(k, grid_g.k_start_, grid_g.k_end_)
	{
		for (int j = grid_g.j_start_; j <= grid_g.j_end_; ++j) for (int i = grid_g.i_start_, ix = grid_g.get1DIndex(i,j,k); i <= grid_g.i_end_; ++i, ++ix)
			if(boundary_array[ix] >= 0) full_ix++;
	}
	END_PARALLEL_FOR;

	// correct full cell index for multi-threading and count number of non-zero elements of A matrix
	int nnz(0);
	int i_res = grid_g.i_res_;
	int ij_res = grid_g.i_res_*grid_g.j_res_;
	int start_full_ix, end_full_ix;

	mt->syncDomainIndices1D(thread_id, full_ix, start_full_ix, end_full_ix);
	BEGIN_PARALLEL_FOR(k, grid_g.k_start_, grid_g.k_end_)
	{
		for(int j=grid_g.j_start_; j <= grid_g.j_end_; ++j) for(int i=grid_g.i_start_, ix=grid_g.get1DIndex(i,j,k); i <= grid_g.i_end_; ++i, ++ix)
		{
			if(boundary_array[ix] >= 0)
			{
				boundary_array[ix] = start_full_ix ++;

				nnz ++;
				if(boundary_array[ix+1]      >= 0) nnz ++;
				if(boundary_array[ix-1]      >= 0) nnz ++;
				if(boundary_array[ix+i_res]  >= 0) nnz ++;
				if(boundary_array[ix-i_res]  >= 0) nnz ++;
				if(boundary_array[ix+ij_res] >= 0) nnz ++;
				if(boundary_array[ix-ij_res] >= 0) nnz ++;
			}
		}
	}
	END_PARALLEL_FOR;

	assert(start_full_ix-1 == end_full_ix);

	full_ix = mt->syncSum(full_ix);
	nnz = mt->syncSum(nnz);

	ONE_THREAD_WORK(A_matrix.Initialize(mt, full_ix, nnz));
	ONE_THREAD_WORK(x_vector.Initialize(full_ix, true));
	ONE_THREAD_WORK(b_vector.Initialize(full_ix));
}

/*
void ProjectionUniform3D::DeterminePressureCPUConjugateGradientMethod(const int& thread_id, const SIM_DOMAIN_UNIFORM_3D& grid_, int* boundary_array_, T* divergence_array_, T* pressure_array_)
{
	InitializeLinearSystem(thread_id, grid_, boundary_array_, A_matrix_, x_vector_, b_vector_);

	const T dt = params_.dt_;
	const T dt_dt = dt*dt;

	const T safe_dist = grid_.dx_w_*(T)5;

	BEGIN_HEAD_THREAD_WORK
	{
		A_matrix_.PrepareForMultithreading();
	}
	END_HEAD_THREAD_WORK;

	int i_res = grid_.cell_num_i_;
	int ij_res = grid_.cell_num_i_*grid_.cell_num_j_;

	BEGIN_PARALLEL_FOR(k, grid_.k_start_g_, grid_.k_end_g_)
	{
		for(int j=grid_.j_start_g_;j<=grid_.j_end_g_;++j) for(int i=grid_.i_start_g_,ix=grid_.Get1DIndexCell(i,j,k);i<=grid_.i_end_g_;++i,++ix)
		{
			if(boundary_array_[ix] < 0) continue;

			const int bc_ijk = boundary_array_[ix];

			int coeff_ijk = 0; // for optimization, inv_dxdx is multiplied at the end.

			// if neighbor is full cell
			if(boundary_array_[ix+1] >= 0)
			{
				coeff_ijk += 1;
				A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array_[ix+1], -1);
			}
			if(boundary_array_[ix-1] >= 0)
			{
				coeff_ijk += 1;
				A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array_[ix-1], -1);
			}


			if(boundary_array_[ix+i_res] >= 0)
			{
				coeff_ijk += 1;
				A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array_[ix+i_res], -1);
			}
			if(boundary_array_[ix-i_res] >= 0)
			{
				coeff_ijk += 1;
				A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array_[ix-i_res], -1);
			}


			if(boundary_array_[ix+ij_res] >= 0)
			{
				coeff_ijk += 1;
				A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array_[ix+ij_res], -1);
			}
			if(boundary_array_[ix-ij_res] >= 0)
			{
				coeff_ijk += 1;
				A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array_[ix-ij_res], -1);
			}

			if(boundary_array_[ix+1]      == BOUNDARY_DIRICHLET) coeff_ijk += 2;
			if(boundary_array_[ix-1]      == BOUNDARY_DIRICHLET) coeff_ijk += 2;
			if(boundary_array_[ix+i_res]  == BOUNDARY_DIRICHLET) coeff_ijk += 2;
			if(boundary_array_[ix-i_res]  == BOUNDARY_DIRICHLET) coeff_ijk += 2;
			if(boundary_array_[ix+ij_res] == BOUNDARY_DIRICHLET) coeff_ijk += 2;
			if(boundary_array_[ix-ij_res] == BOUNDARY_DIRICHLET) coeff_ijk += 2;

			if(coeff_ijk == 0){assert(false);}// no connected FULL and DIR cells, all Neumann condition (null space)

			A_matrix_.AssignValue(thread_id, bc_ijk, bc_ijk, (T)coeff_ijk);

			b_vector_[bc_ijk] = -divergence_array_[ix];// Note: We solve -Ax = -b
//			x_vector_[bc_ijk] = pressure_array_[ix];
			x_vector_[bc_ijk] = 0;
		}
	}
	END_PARALLEL_FOR;

	if (use_cuda_)
	{
	 	BEGIN_HEAD_THREAD_WORK
	 	{
			GPUCGSolve(A_matrix_.row_ptr_, A_matrix_.column_index_, A_matrix_.values_, x_vector_.values_, b_vector_.values_, x_vector_.num_dimension_, x_vector_.num_dimension_, A_matrix_.nz_, 0);
		}
	 	END_HEAD_THREAD_WORK;
	}
	else
	{
		CPUCGSolve(thread_id, A_matrix_, x_vector_, b_vector_);
	}

	BEGIN_PARALLEL_FOR(k, grid_.k_start_, grid_.k_end_)
	{
		for(int j=grid_.j_start_;j<=grid_.j_end_;++j) for(int i=grid_.i_start_,ix=grid_.Get1DIndexCell(i,j,k);i<=grid_.i_end_;++i,++ix)
		{
			if(boundary_array_[ix] >= 0) pressure_array_[ix] = x_vector_[boundary_array_[ix]];

			if(pressure_array_[ix]*dt_dt > safe_dist) pressure_array_[ix] = (T)0;
		}
	}
	END_PARALLEL_FOR;
}
*/

void ProjectionUniform3D::CPUCGSolve(MT* mt, const int& thread_id, const CSRMatrix<T>& A, VectorND<T>& x, const VectorND<T>& b)
{
// 	BEGIN_ONE_THREAD_WORK
// 	{
// 		multithreading_->SplitDomainIndex1D(0, A.N_);
// 	}
// 	END_ONE_THREAD_WORK;

	const Vector2D<int> A_range = mt->getParallelRange(thread_id, 0, A.N_ - 1);

	const int N(x.num_dimension_), start_ix(A_range.t_min_), end_ix(A_range.t_max_);

	ONE_THREAD_WORK(res_.Initialize(N));
	ONE_THREAD_WORK(p_.Initialize(N));
	ONE_THREAD_WORK(Ap_.Initialize(N));
	ONE_THREAD_WORK(num_iteration_ = 0);

	T *rval(res_.values_), *pval(p_.values_), *Apval(Ap_.values_), *xval(x.values_);

	T alpha, res_old, res_new;

	A.ComputeResidual(mt, thread_id, x, b, res_);

	for (int i = start_ix; i <= end_ix; i++)
	{
		p_.values_[i] = res_.values_[i];
	}
	mt->sync();

	res_old = DotProduct(mt, thread_id, res_, p_);

	if (res_old < sqr_tolerance_)
	{
		num_iteration_ = 0;
		return;
	}

	num_iteration_ = 0;
	while (num_iteration_ < max_iteration_)
	{
		A.Multiply(mt, thread_id, p_, Ap_);

		alpha = res_old / DotProduct(mt, thread_id, p_, Ap_);

		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] += alpha * pval[i];
			rval[i] -= alpha * Apval[i];
		}
		mt->sync();

		res_new = DotProduct(mt, thread_id, res_, res_);

		if (res_new < sqr_tolerance_) break;			// L2 norm

		for (int i = start_ix; i <= end_ix; i++)
		{
			const T k = res_new / res_old;

			pval[i] = res_.values_[i] + k*p_.values_[i];
		}
		mt->sync();

		res_old = res_new;

		BEGIN_ONE_THREAD_WORK
		{
			num_iteration_++;
		}
		END_ONE_THREAD_WORK;
	}

	mt->sync();

	BEGIN_ONE_THREAD_WORK
	{
		residual_ = sqrt(res_new);
//		std::cout<<"Iteration = "<<num_iteration_<<" residual = "<<residual_<<std::endl;
	}
	END_ONE_THREAD_WORK;
}

/*
void ProjectionUniform3D::AccumulateSolidForcesByPressure(const int& thread_id, const PROJECTION_MODE mode)
{
	const T *pressure_array_ = data_.pressure_array_;
	const int *boundary_array_ = data_.boundary_array_;

	const int i_res = domain_.cell_num_i_;
	const int ij_res = domain_.cell_num_i_ * domain_.cell_num_j_;

	std::list<SIM_OBJECT_3D*>::iterator it = data_.object_list_.begin();

	for(;it != data_.object_list_.end(); it++)
	{
		if((*it)->object_type_ == SIM_OBJECT_TYPE_RIGID)
		{
			RIGID_OBJECT_3D* rigid_object = (RIGID_OBJECT_3D*)(*it);

			const BOX_3D<T>& aabb_box = rigid_object->sim_aabb_;
			BOX_3D<int> itr_box;
			itr_box.i_start_ = (int)aabb_box.i_start_; itr_box.i_end_ = (int)aabb_box.i_end_; 
			itr_box.j_start_ = (int)aabb_box.j_start_; itr_box.j_end_ = (int)aabb_box.j_end_; 
			itr_box.k_start_ = (int)aabb_box.k_start_; itr_box.k_end_ = (int)aabb_box.k_end_; 

			T pressure_coeff = (mode == PROJECT_WATER) ? rigid_object->options_.water_pressure_coef_ : rigid_object->options_.air_pressure_coef_;
			pressure_coeff = pressure_coeff * domain_.dx_w_;

			if(pressure_coeff <= (T)0) continue;

			BEGIN_BOX_ITERATION_ARR_3D(itr_box, domain_)
			{
				TV3 ijk_s = TV3((T)i + (T)0.5, (T)j + (T)0.5, (T)k + (T)0.5);

				const T phi = rigid_object->GetSignedDistanceS(ijk_s);

				if(phi <= (T)0) // cell inside object
				{
					TV3 acc_force = TV3();
					T density;

					density = data_.mpm_density_array_[arr_ix+1] + data_.smoke_density_array_[arr_ix+1];
					if(density > (T)0) acc_force += pressure_array_[arr_ix+1]*pressure_coeff*TV3((T)-1,(T)0,(T)0)*density;

					density = data_.mpm_density_array_[arr_ix-1] + data_.smoke_density_array_[arr_ix-1];
					if(density > (T)0 ) acc_force += pressure_array_[arr_ix-1]*pressure_coeff*TV3((T) 1,(T)0,(T)0)*density;

					density = data_.mpm_density_array_[arr_ix+i_res] + data_.smoke_density_array_[arr_ix+i_res];
					if(density > (T)0 )	acc_force += pressure_array_[arr_ix+i_res]*pressure_coeff*TV3((T)0,(T)-1,(T)0)*density;

					density = data_.mpm_density_array_[arr_ix-i_res] + data_.smoke_density_array_[arr_ix-i_res];
					if(density > (T)0 )	acc_force += pressure_array_[arr_ix-i_res]*pressure_coeff*TV3((T)0,(T) 1,(T)0)*density;

					density = data_.mpm_density_array_[arr_ix+ij_res] + data_.smoke_density_array_[arr_ix+ij_res];
					if(density > (T)0 ) acc_force += pressure_array_[arr_ix+ij_res]*pressure_coeff*TV3((T)0,(T)0,(T)-1)*density;

					density = data_.mpm_density_array_[arr_ix-ij_res] + data_.smoke_density_array_[arr_ix-ij_res];
					if(density > (T)0 ) acc_force += pressure_array_[arr_ix-ij_res]*pressure_coeff*TV3((T)0,(T)0,(T) 1)*density;

					rigid_object->force_  += acc_force;
					rigid_object->torque_ += CrossProduct(acc_force, rigid_object->position_-ijk_s);					
				}
			}
			END_GRID_ITERATION;
		}
	}
}
*/

/*
void ProjectionUniform3D::AccumulateDeformableForcesByPressure(const int& thread_id, const PROJECTION_MODE mode)
{
	const T *pressure_array_ = data_.pressure_array_;
	const int *boundary_array_ = data_.boundary_array_;

	const int i_res = domain_.cell_num_i_;
	const int ij_res = domain_.cell_num_i_ * domain_.cell_num_j_;

	std::list<SIM_OBJECT_3D*>::iterator it = data_.object_list_.begin();

	for (; it != data_.object_list_.end(); it++)
	{
		if ((*it)->object_type_ == SIM_OBJECT_TYPE_DEFORMABLE)
		{
			DEFORMABLE_OBJECT_3D* deformable_object = (DEFORMABLE_OBJECT_3D*)(*it);

			const BOX_3D<T>& aabb_box = deformable_object->sim_aabb_;
			BOX_3D<int> itr_box;
			itr_box.i_start_ = (int)aabb_box.i_start_; itr_box.i_end_ = (int)aabb_box.i_end_;
			itr_box.j_start_ = (int)aabb_box.j_start_; itr_box.j_end_ = (int)aabb_box.j_end_;
			itr_box.k_start_ = (int)aabb_box.k_start_; itr_box.k_end_ = (int)aabb_box.k_end_;

			const T pressure_coeff = (mode == PROJECT_WATER) ? deformable_object->options_.water_pressure_coef_ : deformable_object->options_.air_pressure_coef_ * domain_.dx_w_;

			int i_start(itr_box.x_min_), i_end(itr_box.x_max_), j_start(itr_box.y_min_), j_end(itr_box.y_max_), k_start(itr_box.z_min_), k_end(itr_box.z_max_);

			PREPARE_FOR_1D_ITERATION(deformable_object->vertex_array_.num_elements_)
			BEGIN_1D_ITERATION
			{								
				TV3 &pos_i = deformable_object->vertex_array_[p]->x_;
				TV3 &vel_i = deformable_object->vertex_array_[p]->velocity_;
				TV3 &force_i = deformable_object->vertex_array_[p]->force_;

				TV3 ix = TV3((int)pos_i.x_ + (T)0.5, (int)pos_i.y_ + (T)0.5, (int)pos_i.z_ + (T)0.5);
				int arr_ix = domain_.Get1DIndexCell((int)ix.x_, (int)ix.y_, (int)ix.z_);
				const T phi = data_.water_levelset_->GetSignedDistanceS(pos_i);

				if (phi <= (T)0 || mode == PROJECT_AIR) // cell inside object
				{
					TV3 acc_force = TV3();

					int num_fluid_cell = 0;
					if (boundary_array_[arr_ix + 1] >= 0) { acc_force += pressure_array_[arr_ix + 1] * pressure_coeff*TV3((T)-1, (T)0, (T)0); num_fluid_cell++; }

					if (boundary_array_[arr_ix - 1] >= 0) { acc_force += pressure_array_[arr_ix - 1] * pressure_coeff*TV3((T)1, (T)0, (T)0); num_fluid_cell++; }

					if (boundary_array_[arr_ix + i_res] >= 0)	{ acc_force += pressure_array_[arr_ix + i_res] * pressure_coeff*TV3((T)0, (T)-1, (T)0); num_fluid_cell++; }

					if (boundary_array_[arr_ix - i_res] >= 0)	{ acc_force += pressure_array_[arr_ix - i_res] * pressure_coeff*TV3((T)0, (T)1, (T)0); num_fluid_cell++; }

					if (boundary_array_[arr_ix + ij_res] >= 0)	{ acc_force += pressure_array_[arr_ix + ij_res] * pressure_coeff*TV3((T)0, (T)0, (T)-1); num_fluid_cell++; }

					if (boundary_array_[arr_ix - ij_res] >= 0)	{ acc_force += pressure_array_[arr_ix - ij_res] * pressure_coeff*TV3((T)0, (T)0, (T)1); num_fluid_cell++; }

					if (num_fluid_cell > 0)
					{
				//		acc_force /= (T)num_fluid_cell;

						TV related_velocity = acc_force;
						vel_i += related_velocity;
					}
				}
			}
			END_1D_ITERATION;
		}
	}
}
*/

/*
int ProjectionUniform3D::CheckStatus(cusparseStatus_t status, char *msg)
{
	if (status != CUSPARSE_STATUS_SUCCESS)
	{
		fprintf(stderr, "!!!! CUSPARSE %s ERROR \n", msg);
		return 1;
	}

	return 0;
}
*/

/*
int ProjectionUniform3D::GPUCGSolve(int* row_pointers, int* column_indices, float* values, float* x_vector, float* b_vector, int M, int N, int nz, const int device)
{
	//	LOG::Begin("CUDACG_LINEAR_SOLVER::CGSolve()");
	boost::chrono::system_clock::time_point start_time;
	static T accum_time1((T)0), elapsed_time1;
	static T accum_time2((T)0), elapsed_time2;
	static T accum_time3((T)0), elapsed_time3;
	//static int frame(0);
	//frame++;

	//	if(use_detailed_log_) LOG::Begin(start_time, "¡Ù CUDACG_LINEAR_SOLVER::CGSolve()");

	static T euclidean_res_accum((T)0);
	cusparseHandle_t handle = 0;
	cusparseStatus_t status;

	// Get handle to the CUSPARSE context
	status = cusparseCreate(&handle);
	if (CheckStatus(status, (char*)"initialization")) return EXIT_FAILURE;

	//	cudaSetDevice(device);
	cudaSetDevice(dev_id_);

	// Description of the A matrix
	cusparseMatDescr_t descr = 0;
	status = cusparseCreateMatDescr(&descr);
	if (CheckStatus(status, (char*)"cusparseCreateMatDescr")) return EXIT_FAILURE;
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	int *d_col, *d_row;
	float *d_val, *d_x, *d_y;
	float *d_r, *d_p, *d_Ap;

	cudaMalloc((void**)&d_col, nz*sizeof(int));
	cudaMalloc((void**)&d_row, (N + 1)*sizeof(int));
	cudaMalloc((void**)&d_val, nz*sizeof(float));
	cudaMalloc((void**)&d_x, N*sizeof(float));
	cudaMalloc((void**)&d_y, N*sizeof(float));
	cudaMalloc((void**)&d_r, N*sizeof(float));
	cudaMalloc((void**)&d_p, N*sizeof(float));
	cudaMalloc((void**)&d_Ap, N*sizeof(float));

	cudaMemcpy(d_col, column_indices, nz*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, row_pointers, (N + 1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_val, values, nz*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_x, x_vector, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, b_vector, N*sizeof(float), cudaMemcpyHostToDevice);

	// 	if(use_detailed_log_) LOG::End(start_time, elapsed_time1, "¡Ú CUDACG_LINEAR_SOLVER::CGSolve()", accum_time1, frame);
	// 	if(use_detailed_log_) LOG::Begin(start_time, "¡Ù CUDACG_LINEAR_SOLVER::CGSolve()");

	// Example code in GNU Octave http://en.wikipedia.org/wiki/Conjugate_gradient_method
	float res_old, res_new, alpha;
	float sqr_res;

	cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, -1.0, descr, d_val, d_row, d_col, d_x, 1.0, d_r);		// r = -A*x + r;
	cublasScopy(N, d_r, 1, d_p, 1);																						// p=r;
	res_old = cublasSdot(N, d_r, 1, d_r, 1);																			// res_old=r'*r;

	if (std::fabs(res_old) > FLT_EPSILON)
	{
		num_iteration_ = 0;
		while (num_iteration_ < max_iteration_)
		{
			cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, 1.0, descr, d_val, d_row, d_col, d_p, 0.0, d_Ap); // Ap=A*p;
			alpha = res_old / cublasSdot(N, d_p, 1, d_Ap, 1);																	 // alpha=res_old/(p'*Ap);
			cublasSaxpy(N, alpha, d_p, 1, d_x, 1);																			 // x=x+alpha*p;
			cublasSaxpy(N, -alpha, d_Ap, 1, d_r, 1);																			 // r=r-alpha*Ap;
			res_new = cublasSdot(N, d_r, 1, d_r, 1);																		 // res_new=r'*r;

			//		LOG::cout<<"itr:"<<num_iteration<<" res:"<<sqrt(res_new)<<" N:"<<N<<" nz:"<<nz<<" sqrt(res_new/N):"<<sqrt(res_new/(T)N)<<std::endl;

			//		if(res_new < sqr_tolerance_) break;
			sqr_res = res_new / (float)N;
			if (sqr_res < sqr_tolerance_) break;

			cublasSscal(N, res_new / res_old, d_p, 1);	// p=r+ress_new/rsold*p;
			cublasSaxpy(N, 1.0, d_r, 1, d_p, 1);
			res_old = res_new;							// res_old=res_new;

			// 		if(use_detailed_log_ && num_iteration_%100==0)
			// 		{
			// 			T Euclidean_res = sqrt(sqr_residual);
			// //LOG::cout<<"[CUDACG] res:"<<Euclidean_res<<" itr:"<<num_iteration_<<std::endl;
			// 			LOG::cout<<"[CUDACG(float)] res:"<<Euclidean_res<<" itr:"<<num_iteration_<<" sqr_res:"<<sqr_residual<<" res_new:"<<res_new<<" N:"<<N<<std::endl;
			// 		}

			num_iteration_++;
		}
	}

	//if (use_detailed_log_)
	//{
	//		T euclidean_res = sqrt(sqr_res);
	//		euclidean_res_accum += euclidean_res;
	//		LOG::cout << "[CUDACG(float)] res:" << euclidean_res << " itr:" << num_iteration_ << " accm:" << euclidean_res_accum << " n_frame:" << frame << " ave:" << euclidean_res_accum / (T)frame << " n_full:" << N << std::endl;
	//}
	//else
	//{
	//	LOG::cout << "Iteration = " << num_iteration_ << " residual = " << sqrt(res_new) << std::endl;
	//}

	//	if(use_detailed_log_) LOG::End(start_time, elapsed_time2, "¡Ú CUDACG_LINEAR_SOLVER::CGSolve()", accum_time2, frame);
	//	if(use_detailed_log_) LOG::Begin(start_time, "¡Ù CUDACG_LINEAR_SOLVER::CGSolve()");

	cudaMemcpy(x_vector, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);
	cusparseDestroy(handle);
	cudaFree(d_col);
	cudaFree(d_row);
	cudaFree(d_val);
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_Ap);

	//	LOG::End();
	//	if(use_detailed_log_) LOG::End(start_time, elapsed_time3, "¡Ú CUDACG_LINEAR_SOLVER::CGSolve()", accum_time3, frame);
	//MULTIGRID_PROJECTION_UNIFORM_3D::coarsest_elapsed_time_ = elapsed_time1 + elapsed_time2 + elapsed_time3;

	return 1;
}
*/

/*
int ProjectionUniform3D::GPUCGSolve(int* row_pointers, int* column_indices, double* values, double* x_vector, double* b_vector, int M, int N, int nz, const int device)
{
	// LOG::Begin("CUDACG(Double)Solve");

	//Copy to double arrays
	// TODO : make this to use multithreading functionality.
	double* values_d = new double[nz];
	double* x_vector_d = new double[N];
	double* b_vector_d = new double[N];
	int i;
	for (i = 0; i<nz; i++)
	{
		values_d[i] = (double)values[i];
	}
	for (i = 0; i<N; i++)
	{
		x_vector_d[i] = (double)x_vector[i];
		b_vector_d[i] = (double)b_vector[i];
	}

	static double Euclidean_res_accum((T)0), abs_residual_accum((T)0);
	static int cnt(0);

	cusparseHandle_t handle = 0;
	cusparseStatus_t status;

	// Get handle to the CUSPARSE context
	status = cusparseCreate(&handle);
	if (CheckStatus(status, (char*)"initialization")) return EXIT_FAILURE;

	cudaSetDevice(dev_id_);

	// Description of the A matrix
	cusparseMatDescr_t descr = 0;
	status = cusparseCreateMatDescr(&descr);
	if (CheckStatus(status, (char*)"cusparseCreateMatDescr")) return EXIT_FAILURE;
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	int *d_col, *d_row;
	double *d_val, *d_x;
	double *d_r, *d_p, *d_Ap;

	cudaMalloc((void**)&d_col, nz*sizeof(int));
	cudaMalloc((void**)&d_row, (N + 1)*sizeof(int));

	cudaMalloc((void**)&d_val, nz*sizeof(double));//
	cudaMalloc((void**)&d_x, N*sizeof(double));//	
	cudaMalloc((void**)&d_r, N*sizeof(double));//

	cudaMalloc((void**)&d_p, N*sizeof(double));
	cudaMalloc((void**)&d_Ap, N*sizeof(double));

	cudaMemcpy(d_col, column_indices, nz*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_row, row_pointers, (N + 1)*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(d_val, values_d, nz*sizeof(double), cudaMemcpyHostToDevice);//
	cudaMemcpy(d_x, x_vector_d, N*sizeof(double), cudaMemcpyHostToDevice);//
	cudaMemcpy(d_r, b_vector_d, N*sizeof(double), cudaMemcpyHostToDevice);//

	// Example code in GNU Octave http://en.wikipedia.org/wiki/Conjugate_gradient_method
	double res_old, res_new, alpha;
	double sqr_residual;

	cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, -1.0, descr, d_val, d_row, d_col, d_x, 1.0, d_r);		    // r = -A*x + r;
	cublasDcopy(N, d_r, 1, d_p, 1);																							// p=r;
	res_old = cublasDdot(N, d_r, 1, d_r, 1);																				// res_old=r'*r;

	num_iteration_ = 0;
	while (num_iteration_ < max_iteration_)
	{
		cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, 1.0, descr, d_val, d_row, d_col, d_p, 0.0, d_Ap);	// Ap=A*p;
		alpha = res_old / cublasDdot(N, d_p, 1, d_Ap, 1);																		// alpha=res_old/(p'*Ap);
		cublasDaxpy(N, alpha, d_p, 1, d_x, 1);																				// x=x+alpha*p;
		cublasDaxpy(N, -alpha, d_Ap, 1, d_r, 1);																			    // r=r-alpha*Ap;
		res_new = cublasDdot(N, d_r, 1, d_r, 1);																			// res_new=r'*r;

		//		LOG::cout<<"itr:"<<num_iteration<<" res:"<<sqrt(res_new)<<" N:"<<N<<" nz:"<<nz<<" sqrt(res_new/N):"<<sqrt(res_new/(T)N)<<std::endl;

		//		if(res_new < sqr_tolerance_) break;
		sqr_residual = res_new / (T)N;
		if (sqr_residual < sqr_tolerance_) break;

		cublasDscal(N, res_new / res_old, d_p, 1);//d_p+=res_new/res_old														// p=r+ress_new/rsold*p;
		cublasDaxpy(N, 1.0, d_r, 1, d_p, 1);    //d_p += d_r

		res_old = res_new;																									// res_old=res_new;

		//if (use_detailed_log_ && num_iteration_ % 20 == 0)
		//{
		//	double Euclidean_res = sqrt(sqr_residual);
		//	LOG::cout << "[CUDACG(double)] res:" << Euclidean_res << " itr:" << num_iteration_ << std::endl;
		//}

		num_iteration_++;
	}

	//Euclicean residual
	cnt++;

	//if (use_detailed_log_)
	//{
	//	double Euclidean_res = sqrt(sqr_residual);
	//	Euclidean_res_accum += Euclidean_res;
	//	LOG::cout << "[CUDACG(double)] res:" << Euclidean_res << " itr:" << num_iteration_ << " accm:" << Euclidean_res_accum << " n_frame:" << cnt
	//		<< " ave:" << Euclidean_res_accum / (T)cnt << " n_full:" << N << std::endl;
	//}
	//else
	//{
	//	LOG::cout << "Iteration = " << num_iteration_ << " residual = " << sqrt(res_new) << std::endl;
	//}

	cudaMemcpy(x_vector_d, d_x, N*sizeof(double), cudaMemcpyDeviceToHost);

	for (i = 0; i<N; i++)
	{
		x_vector[i] = (float)x_vector_d[i];
	}

	// finalize
	cusparseDestroy(handle);
	cudaFree(d_col);
	cudaFree(d_row);
	cudaFree(d_val);
	cudaFree(d_x);
	cudaFree(d_r);
	cudaFree(d_p);
	cudaFree(d_Ap);

	SAFE_DELETE_ARRAY(values_d);
	SAFE_DELETE_ARRAY(x_vector_d);
	SAFE_DELETE_ARRAY(b_vector_d);

	//LOG::End();
	return 1;
}
*/
