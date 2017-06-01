// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../GENERIC_DEFINITIONS.h"
#include <iostream>
#include <fstream>

template<class TT>
class Array3D
{
public:// essential variables
	union{
		struct{int i_start_, j_start_, k_start_, i_end_, j_end_, k_end_;};
		struct{int ix_start_[3], ix_end_[3];};};

public:// data Array1D
	TT *values_;

public:// speedup variables
	int i_res_, j_res_, k_res_;
	int ij_res_;
	int ijk_res_;

public:
	Array3D(void)
		: values_(0)
	{}

	Array3D(const TV_INT& ijk_start_input, const TV_INT& ijk_res_input, const bool& initialize=false)
		: values_(0)
	{
		Initialize(ijk_start_input.i_, ijk_start_input.j_, ijk_start_input.k_, ijk_res_input.i_, ijk_res_input.j_, ijk_res_input.k_, initialize);
	}

	Array3D(const int& i_start_input, const int& j_start_input, const int &k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input,const bool& initialize=false)
		: values_(0)
	{
		Initialize(i_start_input, j_start_input, k_start_input, i_res_input, j_res_input, k_res_input, initialize);
	}

	~Array3D(void)
	{
		if(values_!=0) delete [] values_;
	};

public:// initialization
	void Initialize(const int& i_start_input, const int& j_start_input, const int &k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const bool& initialize=false)
	{
		if(values_!=0)
		{
			delete [] values_;
			values_= 0;
		}

		i_start_ = i_start_input;
		j_start_ = j_start_input;
		k_start_ = k_start_input;

		i_res_ = i_res_input;
		j_res_ = j_res_input;
		k_res_ = k_res_input;
		
		i_end_ = i_start_ + i_res_ - 1;// for(int i=i_start_;i<=i_end_;i++){...}
		j_end_ = j_start_ + j_res_ - 1;
		k_end_ = k_start_ + k_res_ - 1;
		
		// speed up constants
		ij_res_ = i_res_*j_res_;
		ijk_res_ = i_res_*j_res_*k_res_;

		assert(i_res_>0 && j_res_>0 && k_res_>0);
		values_ = new TT [ijk_res_];

		if(initialize==true) AssignAllValues(TT());
	}

	void reset()
	{
		if (values_ != 0)
		{
			delete[] values_;
			values_ = 0;
		}
	}

public:// indexing
	inline const int Index1D(const int& i, const int& j, const int& k) const // 3D index to 1D Array1D index
	{
		assert(i>=i_start_ && i<=i_end_);
		assert(j>=j_start_ && j<=j_end_);
		assert(k>=k_start_ && k<=k_end_);

		return (i-i_start_) + (j-j_start_)*i_res_ + (k-k_start_)*ij_res_;//TODO: try pointer operation for optimization (check performance)
	}

	inline const int Index1D(const TV_INT& index) const // 3D index to 1D Array1D index
	{
		assert(index.i_>=i_start_ && index.i_<=i_end_);
		assert(index.j_>=j_start_ && index.j_<=j_end_);
		assert(index.k_>=k_start_ && index.k_<=k_end_);

		return (index.i_-i_start_) + (index.j_-j_start_)*i_res_ + (index.k_-k_start_)*ij_res_;//TODO: try pointer operation for optimization (check performance)
	}

	TT& getClamped(int i, int j, int k) const
	{
		i = MAX2(i, i_start_);
		i = MIN2(i, i_end_);

		j = MAX2(j, j_start_);
		j = MIN2(j, j_end_);

		k = MAX2(k, k_start_);
		k = MIN2(k, k_end_);

		return (*this)(i, j, k);
	}

	int Get1DIndex(const int& i, const int& j, const int& k) const
	{
		return Index1D(i,j,k);
	}

	int Get1DIndex(const TV_INT& index) const
	{
		return Index1D(index);
	}

	const TV_INT Get3DIndex(const int& index_1d) const
	{
		const int k = index_1d/ij_res_;
		const int j = (index_1d - k*ij_res_)/j_res_;
		const int i = index_1d - k*ij_res_ - j*i_res_;

		return TV_INT(i,j,k);
	}

	const TV_INT ClampedIndex(const int& i, const int& j, const int& k) const
	{
		return TV_INT(clamp(i,i_start_,i_end_), clamp(j,j_start_,j_end_), clamp(k,k_start_,k_end_));
	}

	const TV_INT ClampedIndex(const TV_INT& ix) const
	{
//		return TV_INT(clamp(ix.i,i_start_,i_end_), clamp(ix.j,j_start_,j_end_), clamp(ix.k,k_start_,k_end_));
		return ClampedIndex(ix.i_, ix.j_, ix.k_);//TODO: check performance compared to the commented line
	}

	inline bool Inside(const TV_INT& ix) const
	{
		if(ix.i_ < i_start_) return false;
		else if(ix.i_ > i_end_) return false;
		else if(ix.j_ < j_start_) return false;
		else if(ix.j_ > j_end_) return false;
		else if(ix.k_ < k_start_) return false;
		else if(ix.k_ > k_end_) return false;
		else return true;
	}

	inline bool Inside(const int& i, const int& j, const int& k) const
	{
		if(i < i_start_) return false;
		else if(i > i_end_) return false;
		else if(j < j_start_) return false;
		else if(j > j_end_) return false;
		else if(k < k_start_) return false;
		else if(k > k_end_) return false;
		else return true;
	}

public:// operator overloadings
	TT& operator () (const int& ix) const
	{
		assert(ix >= 0 && ix <= ijk_res_);

		return values_[ix];
	}

	TT& operator () (const int& i, const int& j, const int& k) const
	{
		assert(i>=i_start_ && i<=i_end_);
		assert(j>=j_start_ && j<=j_end_);
		assert(k>=k_start_ && k<=k_end_);

		//TODO: check performance 1. use Index function call, 2. use pointer operation for indexing.
//		return values_[(i-i_start_) + (j-j_start_)*i_res_ + (k-k_start_)*ij_res_];
		return *(values_+(i-i_start_) + (j-j_start_)*i_res_ + (k-k_start_)*ij_res_);
	}

	inline void Assign(const int& i, const int& j, const int& k, const TT& value)
	{
		*(values_+(i-i_start_) + (j-j_start_)*i_res_ + (k-k_start_)*ij_res_) = value;
	}

	TT& operator()(const TV_INT& ix) const
	{
		//TODO: check performance compared to using values_[(i-i_start_) + (j-j_start_)*i_res_ + (k-k_start_)*ij_res_]
		return (*this)(ix.i_, ix.j_, ix.k_);
	}

	TT& DeviatedX (const int& base, const int& i_deviation) const
	{
		assert((base + i_deviation) >= 0 && (base + i_deviation) <= ijk_res_);

		return values_[base + i_deviation];
	}

	TT& DeviatedY (const int& base, const int& j_deviation) const
	{
		assert((base + j_deviation*i_res_) >= 0 && (base + j_deviation*i_res_) <= ijk_res_);

		return values_[base + j_deviation*i_res_];
	}

	TT& DeviatedZ (const int& base, const int& k_deviation) const
	{
		assert((base + k_deviation*ij_res_) >= 0 && (base + k_deviation*ij_res_) <= ijk_res_);

		return values_[base + k_deviation*ij_res_];
	}

public:// manipulation functions
	void AssignAllValues(const TT& constant)
	{
		for(int w=0;w<ijk_res_;w++) values_[w] = constant;//TODO: multi-threading
	}

	void AssignAllValuesZeroGhost()
	{
		std::memset(values_, 0, (ijk_res_)*sizeof(TT));
	}

	void AssignAllValuesZeroGhost(const int& thread_id)
	{
//		PREPARE_FOR_1D_ITERATION(array_.ijk_res_);
//		std::memset(array_.values_+multithreading_->start_ix_1D_[thread_id], 0, (multithreading_->end_ix_1D_[thread_id]-multithreading_->start_ix_1D_[thread_id]+1)*sizeof(TT));
	}

	void AssignRegionalValues(const TT& constant, const int& i_start_, const int& j_start_, const int& k_start_, const int& i_end_, const int& j_end_, const int& k_end_)
	{
		/*
		for(int k=k_start_;k<=k_end_;k++)
			for(int j=j_start_;j<=j_end_;j++)
				for(int i=i_start_;i<=i_end_;i++)
				{
					values_[index] = Indexl(i,j,k);
				}
		*/		
		for(int k=k_start_;k<=k_end_;k++)
		{
			for(int j=j_start_;j<=j_end_;j++)
			{
				for(int i=i_start_;i<=i_end_;i++)
				{
					values_[Get1DIndex(i,j,k)] = constant;
				}
			}
		}
	}

	void operator *= (const T& constant)
	{
		for(int w=0;w<ijk_res_;w++) values_[w] *= constant;
	}

	void operator += (const T& constant)
	{
		for(int w=0;w<ijk_res_;w++) values_[w] += constant;
	}

	void operator -= (const T& constant)
	{
		for(int w=0;w<ijk_res_;w++) values_[w] -= constant;
	}

public:
	std::ofstream& Write(std::ofstream& os)
	{
		for(int i=0; i<ijk_res_; i++)
		{
			TT& value = values_[i];
			os.write((char*)&value, sizeof(TT));			
		}
		return os;
	}

	std::ifstream& Read(std::ifstream& is)
	{
		for(int i=0; i<ijk_res_; i++)
		{
			TT& value = values_[i];
			is.read((char*)&value, sizeof(TT));	
		}
		return is;
	}
};

