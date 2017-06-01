#pragma once

#include "DataStructure/CSRMatrix.h"
#include "DataStructure/GridUniform3D.h"
#include "Parallelism/MultiThreading.h"

/*
#include <helper_functions.h>
#include <cusparse.h>
#include <cublas.h>
*/

enum PROJECTION_MODE {PROJECT_WATER, PROJECT_AIR};

class ProjectionUniform3D
{
public:
	CSRMatrix<T> A_matrix_;
	VectorND<T>  x_vector_;
	VectorND<T>  b_vector_;

	int num_iteration_;
	int max_iteration_;

	T tolerance_;
	T sqr_tolerance_;
	T residual_;

	bool use_cuda_;
	int dev_id_;

	VectorND<T> res_, p_, Ap_;	// res_: residual vector;

public:
	ProjectionUniform3D(void);
	~ProjectionUniform3D(void);

	void Initialize()
	{
		// Poisson solver parameters
		num_iteration_ = 0;
		max_iteration_ = 20;
		tolerance_ = (T)0.0001;
		sqr_tolerance_ = tolerance_*tolerance_;

		use_cuda_ = false;		//TODO: script option
	}

// 	void Update(const int& thread_id, const PROJECTION_MODE& mode, SIM_DATA_UNIFORM_3D& data_);
// 	void DetermineDivergence(const int& thread_id, const PROJECTION_MODE& mode_, SIM_DATA_UNIFORM_3D& data_);
// 	void DeterminePressureCPUConjugateGradientMethod(const int& thread_id, const SIM_DOMAIN_UNIFORM_3D& grid_, int* boundary_array_, T* divergence_array_, T* pressure_array_);
	void CPUCGSolve(MT* mt, const int& thread_id, const CSRMatrix<T>& A, VectorND<T>& x, const VectorND<T>& b);
 	void InitializeLinearSystem(MT* mt, const int& thread_id, const GridUniform3D& grid_g, int* boundary_array_, CSRMatrix<T>& A_matrix_, VectorND<T>& x_vector_, VectorND<T>& b_vector_);
// 	void UpdateVelocityByPressureGradient(const int& thread_id, SIM_DATA_UNIFORM_3D& data_);
// 	void AccumulateSolidForcesByPressure(const int& thread_id, const PROJECTION_MODE mode = PROJECT_WATER);
// 	void AccumulateDeformableForcesByPressure(const int& thread_id, const PROJECTION_MODE mode = PROJECT_WATER);

public: // CUDA
/*
	int CheckStatus(cusparseStatus_t status, char *msg);
	int GPUCGSolve(int* row_pointers, int* column_indices, float* values, float* x_vector, float* b_vector, int M, int N, int nz, const int device);
	int GPUCGSolve(int* row_pointers, int* column_indices, double* values, double* x_vector, double* b_vector, int M, int N, int nz, const int device);
*/
};
