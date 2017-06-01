// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include <assert.h>
#include "Parallelism/MultiThreading.h"

template<class TT>
class VectorND
{
public:
	int num_dimension_;
	TT *values_;

public:
	VectorND(void)
		: num_dimension_(0), values_(0)
	{};

	VectorND(const int& num)
	{
		values_ = 0;

		Initialize(num);
	}

	VectorND(const VectorND<TT>& vector_)
	{
		values_ = 0;

		Initialize(vector_.num_dimension_, false);

		for(int i=0; i < num_dimension_; i ++)	values_[i] = vector_[i];
	}

	~VectorND(void)
	{
		if(values_ != 0) delete [] values_;
		num_dimension_ = 0;
	}

public:
	void Initialize(const int& num_, const bool initialize = false)
	{
		num_dimension_ = num_;

		SAFE_DELETE(values_);

		if(num_dimension_>0)
		{
			values_ = new TT [num_dimension_];
			if(initialize==true) for(int i=0; i<num_dimension_; i++) values_[i] = TT();
		}
	}

public:
	inline void operator = (const VectorND<TT>& from)
	{
		num_dimension_ = from.num_dimension_;

		SAFE_DELETE(values_);

		values_ = new TT[num_dimension_];

		for(int i = 0; i < num_dimension_; i++) values_[i] = from[i];
	}

	inline TT& operator[](const int& i) const
	{
		return values_[i];
	}

	inline VectorND<TT> operator+(const VectorND & vector) const
	{
		assert(num_dimension_ == vector.num_dimension_);

		VectorND<TT> result(num_dimension_);

		for(int i=0;i<num_dimension_;i++) result[i] = values_[i] + vector.values_[i];

		return result;
	}

	inline VectorND<TT> operator-(const VectorND & vector)const
	{
		assert(num_dimension_ == vector.num_dimension_);

		VectorND<TT> result(num_dimension_);

		for(int i=0;i<num_dimension_;i++) result[i] = values_[i] - vector.values_[i];

		return result;
	}

	T SqrMagnitude() const
    {
		T sum(0);
		for(int i=0;i<num_dimension_;i++)
			sum+=SQUARE(values_[i]);
		return sum;
	}

	T MaxAbs() const
	{
		T max_abs(0);
		for(int i=0;i<num_dimension_;i++) max_abs = max(abs(values_[i]), max_abs);
		return max_abs;
	}
	void operator += (const TT& s)
	{
		for(int i = 0; i < num_dimension_; i++) values_[i] += s;
	}

	void operator -= (const TT& s)
	{
		for(int i = 0; i < num_dimension_; i++) values_[i] -= s;
	}

	void operator *= (const int& s)
	{
		for(int i = 0; i < num_dimension_; i++) values_[i] *= s;
	}

	void operator *= (const T& s)
	{
		for(int i = 0; i < num_dimension_; i++) values_[i] *= s;
	}

	void operator /= (const TT& s)
	{
		for(int i = 0; i < num_dimension_; i++) values_[i] /= s;
	}

	void operator += (const VectorND<TT>& s)
	{
		assert(num_dimension_ == s.num_dimension_);

		for(int i = 0; i < num_dimension_; i++) values_[i] += s.values_[i];
	}

	void operator -= (const VectorND<TT>& s)
	{
		assert(num_dimension_ == s.num_dimension_);
		
		for(int i = 0; i < num_dimension_; i++) values_[i] -= s.values_[i];
	}

	void operator *= (const VectorND<TT>& s)
	{
		assert(num_dimension_ == s.num_dimension_);

		for(int i = 0; i < num_dimension_; i++) values_[i] *= s.values_[i];
	}

	void operator /= (const VectorND<TT>& s)
	{
		assert(num_dimension_ == s.num_dimension_);

		for(int i = 0; i < num_dimension_; i++) values_[i] /= s.values_[i];
	}

	VectorND <TT> operator*(const T& s) const
	{
		VectorND <TT> V(num_dimension_);

		for(int i=0; i < num_dimension_; i++) V.values_[i]=values_[i]*s;

		return V;
	}
};

static D dotProduct(const VectorND<D>& v1, const VectorND<D>& v2)
{
	assert(v1.num_dimension_ == v2.num_dimension_);

	T sum = 0;

	for (int i = 0; i < v1.num_dimension_; i++)
		sum += v1.values_[i] * v2.values_[i];

	return sum;
}


inline static T dotProduct(const VectorND<T>& v1,const VectorND<T>& v2)
{
	assert(v1.num_dimension_ == v2.num_dimension_);

	T sum=0;

	for(int i=0; i<v1.num_dimension_; i++)
		sum += v1.values_[i]*v2.values_[i];

	return sum;
}

inline static T dotProduct(const VectorND<T>& v1, const VectorND<T>& v2, const int& start_ix, const int& end_ix)
{
	assert(v1.num_dimension_ == v2.num_dimension_);
	
	T *v1_values = v1.values_, *v2_values = v2.values_;

	T sum(0);

	for(int i = start_ix; i <= end_ix; i++) sum += v1_values[i] * v2_values[i];

	return sum;
}

template<class TT>
TT dotProduct(MT* mt, const int& thread_id, const VectorND<TT>& v1, const VectorND<TT>& v2)
{
	assert(v1.num_dimension_ == v2.num_dimension_);

	TT *v1_values = v1.values_, *v2_values = v2.values_;

	const int start_ix(mt->start_ix_1D_[thread_id]), end_ix(mt->end_ix_1D_[thread_id]);

	TT dot_product = (TT)0;

	for(int i = start_ix; i <= end_ix; i++) dot_product += v1_values[i] * v2_values[i];

	dot_product = mt->syncSum(thread_id, dot_product);

	return dot_product;
}

template<class T> inline std::ostream&
operator<<(std::ostream& output,const VectorND<T>& v)
{
	for(int i=0; i<v.num_dimension_; i++) output<<v.values_[i]<<" ";
	output << std::endl;
	return output;
}
