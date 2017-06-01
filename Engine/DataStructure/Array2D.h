// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <assert.h>
#include <list>

#include "GENERIC_DEFINITIONS.h"

template<class TT>
class Array2D
{
public:
	union{
		struct{int i_start_, j_start_, i_end_, j_end_;};
		struct{int ix_start_[2], ix_end_[2];};};

public:// data Array1D
	TT *values_;

public:// speedup variables
	int i_res_, j_res_;
	int ij_res_;

public:
	Array2D(void)
		: values_(0)
	{}

	Array2D(const TV2_INT& ij_start_input, const TV2_INT& ij_res_input, const bool& _initialize=false)
		: values_(0)
	{
		initialize(ij_start_input.i_, ij_start_input.j_, ij_res_input.i_, ij_res_input.j_, _initialize);
	}

	Array2D(const int& i_start_input, const int& j_start_input, const int& i_res_input, const int& j_res_input, const bool& _initialize = false)
		: values_(0)
	{
		initialize(i_start_input, j_start_input, i_res_input, j_res_input, _initialize);
	}

	Array2D(const Array2D<TT>& input_arr)
		: values_(0)
	{
		initialize(input_arr.i_start_, input_arr.j_start_, input_arr.i_res_, input_arr.j_res_, false);

		for (int p = 0; p < ij_res_; ++p)
			values_[p] = input_arr.values_[p];
	}

	~Array2D(void)
	{
		if(values_!=0) delete [] values_;
	};

public:// initialization
	void initialize(const int& i_start_input, const int& j_start_input, const int& i_res_input, const int& j_res_input, const bool& _initialize = false)
	{
		if(values_!=0) delete [] values_;

		i_start_ = i_start_input;
		j_start_ = j_start_input;

		i_res_ = i_res_input;
		j_res_ = j_res_input;

		i_end_ = i_start_ + i_res_ - 1;// for(int i=i_start;i<=i_end_;i++){...}
		j_end_ = j_start_ + j_res_ - 1;
		
		// speed up constants
		ij_res_ = i_res_*j_res_;

		assert(i_res_>0 && j_res_>0);
		values_ = new TT [ij_res_];

		if (_initialize == true) assignAllValues(TT());
	}

	void swapIJ()
	{
		SWAP(i_start_, j_start_, int);
		SWAP(i_res_, j_res_, int);
		SWAP(i_end_, j_end_, int);
	}

	int getNumAllValues() const
	{
		return i_res_ * j_res_;
	}

	void freeMemory()
	{
		initialize(0, 0, 0, 0, false);

		SAFE_DELETE_ARRAY(values_);
	}

	void allocateMemory()
	{
		SAFE_DELETE_ARRAY(values_);
		values_ = new TT [ij_res_];
	}

public:// indexing
	const int get1DIndex(const int& i, const int& j) const // 2D index to 1D Array1D index
	{
		assert(i>=i_start_ && i<=i_end_);
		assert(j>=j_start_ && j<=j_end_);

		return (i-i_start_) + (j-j_start_)*i_res_;//TODO: try pointer operation for optimization (check performance)
	}

	const int get1DIndex(const TV2_INT& index) const // 2D index to 1D Array1D index
	{
		assert(index.i_>=i_start_ && index.i_<=i_end_);
		assert(index.j_>=j_start_ && index.j_<=j_end_);

		return (index.i_-i_start_) + (index.j_-j_start_)*i_res_;//TODO: try pointer operation for optimization (check performance)
	}

	const bool isValid(const int& i, const int& j) const
	{
		// don't need to check the rage of i and j because this function is made for it
		// assert(i >= i_start_ && i <= i_end_);
		// assert(j >= j_start_ && j <= j_end_);

		if (i < i_start_) return false;
		else if (i > i_end_) return false;
		else if (j < j_start_) return false;
		else if (j > j_end_) return false;
		else return true;
	}


	inline const TV2_INT getClampedIndex(const int& i, const int& j) const
	{
		return TV2_INT(clamp(i,i_start_,i_end_), clamp(j,j_start_,j_end_));
	}

	inline const TV2_INT getClampedIndex(const TV2_INT& ix) const
	{
		return getClampedIndex(ix.i_, ix.j_);//TODO: check performance compared to the commented line
	}

public:// overloaded operators
	TT& operator () (const int& ix) const
	{
		assert(ix >= 0 && ix <= ij_res_);

		return values_[ix];
	}

	TT& operator () (const int& i, const int& j) const
	{
		assert(i>=i_start_ && i<=i_end_);
		assert(j>=j_start_ && j<=j_end_);

		//TODO: check performance 1. use Index function call, 2. use pointer operation for indexing.

		return *(values_+(i-i_start_) + (j-j_start_)*i_res_);
	}

	TT& getClamped(int i, int j) const
	{
		i = MAX2(i, i_start_);
		i = MIN2(i, i_end_);

		j = MAX2(j, j_start_);
		j = MIN2(j, j_end_);

		return (*this)(i, j);		
	}

	TT& getIRepeatedJClamped(int i, int j) const
	{
		j = MAX2(j, j_start_);
		j = MIN2(j, j_end_);

		if (i < i_start_) i = i_end_;
		else if (i > i_end_) i = i_start_;

		return (*this)(i, j);
	}

	TT& getIClampedJRepeated(int i, int j) const
	{
		i = MAX2(i, i_start_);
		i = MIN2(i, i_end_);

		if (j < j_start_) j = j_end_;
		else if (j > j_end_) j = j_start_;

		return (*this)(i, j);
	}

	TT& operator () (const TV2_INT& ix) const
	{
		return (*this)(ix.i_, ix.j_);
	}

	TT& getDeviatedX (const int& base, const int& i_deviation) const
	{
		assert((base + i_deviation) >= 0 && (base + i_deviation) <= ij_res_);

		return values_[base + i_deviation];
	}

	TT& getDeviatedY (const int& base, const int& j_deviation) const
	{
		assert((base + j_deviation*i_res_) >= 0 && (base + j_deviation*i_res_) <= ij_res_);

		return values_[base + j_deviation*i_res_];
	}

public:// manipulation functions
	void assignAllValues(const TT& constant)
	{
		for(int w=0;w<ij_res_;w++) values_[w] = constant;//TODO: multi-threading
	}

// 	void AssignAllValues(MULTITHREADING* multithreading_, const int& thread_id, const TT& constant)
// 	{
// 		PREPARE_FOR_1D_ITERATION(ij_res_);
// 		std::memset(values_+multithreading_->start_ix_1D_[thread_id], constant, (multithreading_->end_ix_1D_[thread_id]-multithreading_->start_ix_1D_[thread_id]+1)*sizeof(TT));
// 		multithreading_->Sync(thread_id);
// 	}
// 
// 	void AssignAllValuesZero(MULTITHREADING* multithreading_, const int& thread_id)
// 	{
// 		PREPARE_FOR_1D_ITERATION(ij_res_);
// 		std::memset(values_+multithreading_->start_ix_1D_[thread_id], 0, (multithreading_->end_ix_1D_[thread_id]-multithreading_->start_ix_1D_[thread_id]+1)*sizeof(TT));
// 		multithreading_->Sync(thread_id);
// 	}

	void assignRegionalValues(const TT& constant, const int& i_start, const int& j_start, const int& i_end_, const int& j_end_)
	{
		for(int j=j_start;j<=j_end_;j++)
		{
			for(int i=i_start;i<=i_end_;i++)
			{
				values_[get1DIndex(i,j)] = constant;
			}
		}
	}

	void operator *= (const T& constant)
	{
		for(int w=0;w<ij_res_;w++) values_[w] *= constant;
	}

	void operator += (const T& constant)
	{
		for(int w=0;w<ij_res_;w++) values_[w] += constant;
	}

	void operator -= (const T& constant)
	{
		for(int w=0;w<ij_res_;w++) values_[w] -= constant;
	}

	void copyFrom(const Array2D<TT>& input)
	{
		assert(this->ij_res_ == input.ij_res_);

		for (int w = 0; w < ij_res_; w++) values_[w] = input.values_[w];
	}

	void extendOneColumn(const int column, const int extend_width)
	{
		Array2D<TT> temp;
		temp.initialize(i_start_, j_start_, i_res_ + extend_width, j_res_, false);

		for (int j = j_start_; j <= j_end_; ++j)
		for (int i = i_start_; i < column; ++i)
		{
			temp(i, j) = (*this)(i, j);
		}

		for (int j = j_start_; j <= j_end_; ++j)
		for (int i = column; i <= column + extend_width; ++i)
		{
			temp(i, j) = (*this)(column, j);
		}

		for (int j = j_start_; j <= j_end_; ++j)
		for (int i = column+1; i <= i_end_; ++i)
		{
			temp(i + extend_width, j) = (*this)(i, j);
		}

		initialize(temp.i_start_, temp.j_start_, temp.i_res_, temp.j_res_, false);

		copyFrom(temp);
	}

	void FloodFill(const int i_start, const int j_start, const TT stop_flag, const TT new_flag)
	{
		std::list<TV2_INT> ij_stack;

		ij_stack.push_back(TV2_INT(i_start, j_start));

		while (!ij_stack.empty())
		{
			const TV2_INT ij = ij_stack.back();

			ij_stack.pop_back();

			const int &i = ij.i_, &j = ij.j_;

			if ((*this)(ij.i_, ij.j_) == stop_flag || (*this)(ij.i_, ij.j_) == new_flag) continue;

			(*this)(i, j) = new_flag;

			if (i != i_start_) ij_stack.push_back(TV2_INT(i - 1, j));

			if (i != i_end_) ij_stack.push_back(TV2_INT(i + 1, j));

			if (j != j_start_) ij_stack.push_back(TV2_INT(i, j - 1));

			if (j != j_end_)  ij_stack.push_back(TV2_INT(i, j + 1));
		}		
	}

	void smoothLaplacian(const T& alpha, const int& repeat)
	{
		for (int r = 0; r < repeat; r++)
		{
			for (int j = j_start_; j <= j_end_; j++)
				for (int i = i_start_; i <= i_end_; i++)
				{
					const T average = (getClamped(i - 1, j) + getClamped(i + 1, j) + getClamped(i, j - 1) + getClamped(i, j + 1)) * (T)0.25;

					(*this)(i, j) = (*this)(i, j) * ((T)1 - alpha) + average * alpha;
				}
		}
	}

	void smoothMaxClamped(const T& alpha, const T& max, const int& repeat)
	{
		for (int r = 0; r < repeat; r++)
		{
			for (int j = j_start_; j <= j_end_; j++)
				for (int i = i_start_; i <= i_end_; i++)
				{
					if ((*this)(i, j) >= max) continue;

					const T average = (getClamped(i - 1, j) + getClamped(i + 1, j) + getClamped(i, j - 1) + getClamped(i, j + 1)) * (T)0.25;

					(*this)(i, j) = (*this)(i, j) * ((T)1 - alpha) + average * alpha;
				}
		}
	}

	void clampMinMaxZeroOne(const TT& min, const TT& max)
	{
		for (int j = j_start_; j <= j_end_; j++)
			for (int i = i_start_; i <= i_end_; i++)
			{
				T height = (*this)(i, j);

				if (height <= min) height = 0.0f;
				if (height > max) height = 1.0f;

				(*this)(i, j) = height;
			}
	}

public:
	void read()
	{}

	void write()
	{}
};

