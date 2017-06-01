// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "DataStructure/Array1D.h"

class MergeSort
{
public:
	Array1D<int> ix_array_;
	Array1D<int> sorted_;

	// text book sorting
	void sort(Array1D<int>& arr);
	void sortArray(int list[], int left, int right);
	void merge(int list[], int left, int mid, int right);

	// key sorting
	template<class TT>
	void sortWithKey(Array1D<TT>& key_arr)
	{
		ix_array_.initialize(key_arr.num_elements_);
		sorted_.initialize(key_arr.num_elements_);

		for (int i = 0; i < ix_array_.num_elements_; i++) ix_array_.values_[i] = i;			// initialize index

		sortArrayWithKey(key_arr.values_, ix_array_.values_, 0, key_arr.num_elements_ - 1);
//		sortArrayWithKeyNonrecursive(key_arr, ix_array_.values_);
	}

	template<class TT>
	void mergeWithKey(const TT* key_arr, int* list, int left, int mid, int right)
	{
		int i = left, j = mid + 1, k = left;

		// sort while merging
		while (i <= mid && j <= right)
		{
			if (key_arr[list[i]] <= key_arr[list[j]])		// use key to compare
				sorted_.values_[k++] = list[i++];
			else
				sorted_.values_[k++] = list[j++];
		}

		// copy rest 
		if (i > mid)
		{
			for (int l = j; l <= right; l++)
				sorted_.values_[k++] = list[l];
		}
		else
		{
			for (int l = i; l <= mid; l++)
				sorted_.values_[k++] = list[l];
		}

		// copy sorted back to the input array
		for (int l = left; l <= right; l++)
			list[l] = sorted_.values_[l];
	}

	template<class TT>
	void sortArrayWithKey(const TT* key_arr, int* list, int left, int right)
	{
		int mid;
		if (left < right)
		{
			mid = (left + right) / 2;

//			std::cout<<"("<<left<<" "<<right<<")"<<std::endl;

			sortArrayWithKey(key_arr, list, left, mid);			// recursive merge sort of left array
			sortArrayWithKey(key_arr, list, mid + 1, right);		// recursive merge sort of right array
			mergeWithKey(key_arr, list, left, mid, right);	// merge
		}
	}

	template<class TT>
	void sortArrayWithKeyNonrecursive(const Array1D<TT> &key_arr, int* list)
	{
		int spacing = 2;

		while (true)
		{
			for (int i = 0; i < key_arr.num_elements_; i += spacing)
			{			
				const int left = i;
				int right = MIN2(i + spacing - 1, key_arr.num_elements_ - 1);

				if (left < right)
				{
					const int mid = (left + right) / 2;

//					std::cout << "(" << left << " " << right << ")";

					mergeWithKey(key_arr.values_, list, i, mid, right);
				}
			}

//			std::cout << std::endl;

			if ((key_arr.num_elements_ / spacing) == 0) break;

			spacing *= 2;
		}

// 		
// 		for (int i = 0; i < key_arr.num_elements_; i+= spacing)
// 		{
// 			mergeWithKey(key_arr.values_, list, i, spacing*i);
// 		}
// 
// 		int mid;
// 		if (left < right)
// 		{
// 			mid = (left + right) / 2;
// 
// 			std::cout << "(" << left << " " << right << ")" << std::endl;
// 
// 			sortArrayWithKey(key_arr, list, left, mid);			// recursive merge sort of left array
// 			sortArrayWithKey(key_arr, list, mid + 1, right);		// recursive merge sort of right array
// 			mergeWithKey(key_arr, list, left, mid, right);	// merge
// 		}
	}
};