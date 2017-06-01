// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <assert.h>

#define MAX(a, b)							((a) > (b) ? (a) : (b))

template<class TT>
class DynamicArray
{
public:
	int num_elements_;// number of all elements, i.e. last index is num_element - 1.
	int array_size_;// size of val Array1D. array_size >= num_element

	TT *values_;

public:
	int num_resize_;// number of adding memory size when reinitialized due to the lack of memory size.

public:
// 	DYNAMIC_ARRAY(void)
// 		: num_elements_(0), array_size_(0), values_(0), num_resize_(1)
// 	{}

	DynamicArray(const int& array_size = 0, const int& num_resize = 1)
		: num_elements_(0), array_size_(array_size), values_(0), num_resize_(num_resize)
	{
		Initialize(array_size_, num_resize_);
	}

	~DynamicArray(void)
	{
		SAFE_DELETE_ARRAY(values_);
	}

public:
	void Initialize(const int& array_size_input = 0, const int& num_resize_input = 1)
	{
		assert(array_size_input >= 0 && num_resize_input > 0);	// allows the number of Array1D elements to be zero

		array_size_ = array_size_input;
		num_resize_ = num_resize_input;

		if(values_ != 0)
		{
			delete [] values_;
			values_ = 0;
		}

		if(array_size_ > 0)
			values_ = new TT [array_size_];

		num_elements_ = 0;
	}

	void Realocate(const int& array_size_input)
	{
		assert(array_size_input >= num_elements_);	// not to lose existing data

		array_size_ = MAX(array_size_input, 1);		// Array1D size should be 0 at least.

		TT *new_array = new TT [array_size_];		// allocate memory for new Array1D

		for(int i = 0; i < num_elements_; i++) new_array[i] = values_[i];// copy existing data

		delete [] values_;
		values_ = new_array;
	}

	void ResizeAsReservedArray(const int& array_size_input)
	{
		if(array_size_ < array_size_input) Initialize(array_size_input);

		num_elements_ = array_size_input;
	}

	void ReserveAndEmptify(const int& array_size_input)
	{
		if(array_size_ < array_size_input) Initialize(array_size_input);

		num_elements_ = 0;
	}

	void Minimize()// reallocate array_size to num_elements (or + 1)
	{
		assert(num_elements_ <= array_size_);

		if(num_elements_ == array_size_) return;

		TT *new_array = new TT [num_elements_];// reallocate memory for larger Array1D
		for(int i=0; i < num_elements_; i++) new_array[i] = values_[i];// copy existing data
		array_size_ = num_elements_;
		delete [] values_;
		values_ = new_array;
	}

	void Push(const TT& data_input)
	{
		if(array_size_ > num_elements_)
		{
			values_[num_elements_] = data_input;
			num_elements_ ++;
		}
		else
		{
			Realocate(array_size_ + num_resize_);

			values_[num_elements_] = data_input;
			num_elements_ ++;
		}
	}

	TT& Push()
	{
		if(array_size_ > num_elements_) num_elements_ ++;
		else
		{
			Realocate(array_size_ + num_resize_);
			num_elements_ ++;
		}

		return values_[num_elements_-1];
	}

	bool Pop(TT& data_output)
	{
		if(num_elements_ > 0)
		{
			data_output = values_[num_elements_-1];
			num_elements_ --;

			return true;
		}
		else
			return false;
	}

	TT Pop()						//NOTE: check num_elements before use this Pop() function.
	{
		assert(num_elements_ > 0);
		
		return values_[--num_elements_];
	}

	void Emptify()				//NOTE: this function does not return memory to OS.
	{
		num_elements_ = 0;
	}

	inline void Reset()			//NOTE: this function returns memory to OS.
	{
		if(values_ == 0) return;
		
		delete [] values_;

		values_ = 0;

		num_elements_ = 0;// do not delete or reinitialize Array1D for efficiency
		array_size_ = 0;// allows the number of elements to be zero
	}

	inline TT& operator[](const int& i) const
	{
		assert(i >= 0 && i < num_elements_);

		return values_[i];
	}
};

