// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "VectorND.h"
#include "Parallelism/MultiThreading.h"

// compressed sparse row (CSR or CRS) http://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR_or_CRS.29
// see http://developer.download.nvidia.com/compute/DevZone/docs/html/CUDALibraries/doc/CUSPARSE_Library.pdf for CUDA compatable descriptions.

template <class T>
class CSRMatrix
{
public:
	int N_;
	int nz_;
	T *values_;				// the values of the nonzero elements
	int *row_ptr_;			// the locations in the val vector that start a row
	int *column_index_;		// the column indexes of the elements in the val vector
	
	int values_ix_;
	int row_count_;			//row_count is not used.
	int prev_row_;

	int *start_ix_, *end_ix_;
	int *prev_row_array_;
	int *values_ix_array_;

public:
	CSRMatrix(void)
		: values_(0)
		, column_index_(0)
		, row_ptr_(0)
		, start_ix_(0)
		, end_ix_(0)
		, prev_row_array_(0)
		, values_ix_array_(0)
	{}

	~CSRMatrix(void)
	{
		DeleteMemory();
	}

	void DeleteMemory()
	{
		SAFE_DELETE_ARRAY(values_);
		SAFE_DELETE_ARRAY(row_ptr_);
		SAFE_DELETE_ARRAY(column_index_);
		SAFE_DELETE_ARRAY(start_ix_);
		SAFE_DELETE_ARRAY(end_ix_);
		SAFE_DELETE_ARRAY(prev_row_array_);
		SAFE_DELETE_ARRAY(values_ix_array_);
	}

public:
	void Initialize(const int num_threads, const int &N_input, const int &nz_input)// N by N square matrix. nz: number of non zero elements
	{
		DeleteMemory();

		start_ix_ = new int[num_threads];
		end_ix_ = new int[num_threads];
		prev_row_array_ = new int[num_threads];
		values_ix_array_ = new int[num_threads];

		N_ = N_input;
		nz_ = nz_input;

		values_ = new T [nz_];
		row_ptr_ = new int [N_+1];				// +1 is for final row iteration
		column_index_ = new int [nz_];

		values_ix_ = 0;
		row_count_ = 0;
		prev_row_ = -1;

		row_ptr_[N_input] = nz_input;
	}

	void PrepareForMultithreading(MT* mt)
	{
		mt->splitRange(0, N_ - 1);

		start_ix_[0] = 0;
		end_ix_[0] = mt->sync_value_int_[0] - 1;
		
		prev_row_array_[0] = -1;
		values_ix_array_[0] = 0;

		for(int thread_id = 1; thread_id < (int)mt->num_threads_; thread_id++)
		{
			start_ix_[thread_id] = end_ix_[thread_id-1] + 1;
			end_ix_[thread_id] = end_ix_[thread_id-1] + mt->sync_value_int_[thread_id];
			prev_row_array_[thread_id] = -1;
			values_ix_array_[thread_id] = start_ix_[thread_id];
		}
	}

	void AssignValue(const int& row_input, const int& column_input, const T& values_input)
	{
		values_[values_ix_] = values_input;

		if(row_input != prev_row_)
		{
//			assert(row_input == prev_row_+1);

			row_ptr_[row_input] = values_ix_;
			prev_row_ = row_input;
		}

		column_index_[values_ix_] = column_input;

		values_ix_ ++;
	}

	void AssignValue(const int& thread_id, const int& row_input, const int& column_input, const T& values_input)
	{

		values_[values_ix_array_[thread_id]] = values_input;

		if(row_input != prev_row_array_[thread_id])
		{
			//check whether the matrix is well-made or not.
			//it can be omitted in simulation stage.
//			assert(row_input == prev_row_array_[thread_id]+1);

			row_ptr_[row_input] = values_ix_array_[thread_id];
			prev_row_array_[thread_id] = row_input;
		}

		column_index_[values_ix_array_[thread_id]] = column_input;

		values_ix_array_[thread_id] ++;
	}
	
	inline void Multiply(const VectorND<T>& x, VectorND<T>& b) const // this_matrix * x -> b
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);

		T *bval(b.values_), *xval(x.values_);

		for(int row = 0; row < N_; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			T v=0;
			for(int vix = row_ptr_[row]; vix < row_ptr_[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values_[vix]*xval[column_index_[vix]];
			}

			bval[row] = v;
		}
	}

	inline void Multiply(const VectorND<T>& x, VectorND<T>& b, const int& k_start, const int& k_end) const // this_matrix * x -> b
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);
		T *bval(b.values_), *xval(x.values_);

		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			T v=0;
			for(int vix = row_ptr_[row]; vix < row_ptr_[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values_[vix]*xval[column_index_[vix]];
			}

			bval[row] = v;
		}
	}

	void Multiply(MT* mt, const int& thread_id, const VectorND<float>& x, VectorND<float>& b) const // this_matrix * x -> b
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);

		float *bval(b.values_), *xval(x.values_);

		const Vector2D<int> k_range = mt->getParallelRange(thread_id, 0, b.num_dimension_ - 1);
		const int k_start(k_range.t_min_), k_end(k_range.t_max_);
		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			float v = 0;
			const int vix_start = row_ptr_[row];	assert(vix_start >= 0);
			const int vix_end = row_ptr_[row+1];	assert(vix_start < nz_);
			for(int vix = vix_start; vix < vix_end; vix ++) // iterate all components of 'row'th row of this matrix
			{
				assert(column_index_[vix] < N_);

				v += values_[vix]*xval[column_index_[vix]];
			}

			bval[row] = v;
		}

		mt->sync();
	}

	void Multiply(MT* mt, const int& thread_id, const VectorND<double>& x, VectorND<double>& b) const // this_matrix * x -> b
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);

		double *bval(b.values_), *xval(x.values_);

		const Vector2D<int> k_range = mt->getParallelRange(thread_id, 0, b.num_dimension_ - 1);
		const int k_start(k_range.t_min_), k_end(k_range.t_max_);
		for (int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			double v = 0;
			const int vix_start = row_ptr_[row];	assert(vix_start >= 0);
			const int vix_end = row_ptr_[row + 1];	assert(vix_start < nz_);
			for (int vix = vix_start; vix < vix_end; vix++) // iterate all components of 'row'th row of this matrix
			{
				assert(column_index_[vix] < N_);

				v += values_[vix] * xval[column_index_[vix]];
			}

			bval[row] = v;
		}

		mt->sync();
	}

	inline void ComputeResidual(const VectorND<T>& x, const VectorND<T>& b, VectorND<T>& residual) const // residual = b - this_matrix*x
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);
		assert(residual.num_dimension_ == N_);

		T *bval(b.values_), *xval(x.values_), *rval(residual.values_);

		for(int row = 0; row < N_; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			// compute A*x
			T v=0;
			for(int vix = row_ptr_[row]; vix < row_ptr_[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values_[vix]*xval[column_index_[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}
	}

	inline void ComputeResidual(const VectorND<T>& x, const VectorND<T>& b, VectorND<T>& residual, const int& k_start, const int& k_end) const // residual = b - this_matrix*x
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);
		assert(residual.num_dimension_ == N_);

		// speed-up pointers
		T *bval(b.values_), *xval(x.values_), *rval(residual.values_);

		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			// compute A*x
			// TODO: we may optimize row_ptr_[row] and row_ptr_[row+1] access
			T v=0;
			for(int vix = row_ptr_[row]; vix < row_ptr_[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values_[vix]*xval[column_index_[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}
	}

	inline void ComputeResidual(MT* mt, const int& thread_id, const VectorND<T>& x, const VectorND<T>& b, VectorND<T>& residual) const // residual = b - this_matrix*x
	{
		assert(N_ == x.num_dimension_);
		assert(x.num_dimension_ == b.num_dimension_);
		assert(residual.num_dimension_ == N_);

		// speed-up pointers
		T *bval(b.values_), *xval(x.values_), *rval(residual.values_);

		const Vector2D<int> k_range = mt->getParallelRange(thread_id, 0, b.num_dimension_ - 1);
		const int k_start(k_range.t_min_), k_end(k_range.t_max_);
		for(int row = k_start; row <= k_end; row++) // multiply 'row'th row of this matrix and vector x to compute b[row]
		{
			// compute A*x
			// TODO: we may optimize row_ptr_[row] and row_ptr_[row+1] access
			T v=0;
			for(int vix = row_ptr_[row]; vix < row_ptr_[row+1]; vix ++) // iterate all components of 'row'th row of this matrix
			{
				v += values_[vix]*xval[column_index_[vix]];
			}

			// residual = b - A*x
			rval[row] = bval[row] - v;
		}

		mt->sync();
	}

	void operator*=(const T& s)
	{
		for(int i = 0; i < nz_; i++)
		{
			values_[i] *= s;
		}
	}

	T& GetValue(const int& row, const int& column)
	{
		for(int vix = row_ptr_[row]; vix < row_ptr_[row+1]; vix ++)
		{
			if(column_index_[vix] == column) return &values_[vix];
		}

		std::cout<<"CSRMatrix error "<<std::endl;
		exit(-1);
	}
};


