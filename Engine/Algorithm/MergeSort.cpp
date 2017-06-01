// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "MergeSort.h"

void MergeSort::merge(int* list, int left, int mid, int right)
{
	int i = left, j = mid + 1, k = left;

	// sort while merging
	while (i <= mid && j <= right)
	{
		if (list[i] <= list[j])
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

void MergeSort::sortArray(int* list, int left, int right)
{
	int mid;
	if (left < right)
	{
		mid = (left + right) / 2;
		sortArray(list, left, mid);			// recursive merge sort of left array
		sortArray(list, mid + 1, right);		// recursive merge sort of right array
		merge(list, left, mid, right);	// merge
	}
}

void MergeSort::sort(Array1D<int>& arr)
{
	sorted_.initialize(arr.num_elements_);

	sortArray(arr.values_, 0, arr.num_elements_ - 1);
}