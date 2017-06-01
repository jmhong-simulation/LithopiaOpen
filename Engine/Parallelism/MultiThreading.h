// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <atomic>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>

#include "DataStructure/Vector2D.h"
#include "DataStructure/Array1D.h"
#include "Utilities/RandomNumberGenerator.h"
#include <boost/chrono.hpp>

class MultiThreading
{
public:
	unsigned int num_threads_;

	std::condition_variable cond_;
	std::mutex mtx_lock_;

	std::once_flag o_flag;

	Array1D<std::thread*> thread_list_;

	Array1D<int> start_ix_1D_, end_ix_1D_;

	std::atomic<int> num_waiting_threads_;

	std::atomic<int> num_of_synced_threads_;

	float *sync_value_float_;
	int *sync_value_int_;
	double *sync_value_double_;

	float sync_value_float_temp_;
	double sync_value_double_temp_;

	Array1D<RandomNumberGenerator*> random_;

public:
	MultiThreading(const unsigned int _num_threads = std::thread::hardware_concurrency())
		: num_threads_(0), num_waiting_threads_(0), sync_value_float_(nullptr), sync_value_int_(nullptr), sync_value_double_(nullptr), num_of_synced_threads_(0)
	{
		initialize(_num_threads);
	}

	~MultiThreading()
	{
		deleteAllThreads();

		for (int i = 0; i < random_.num_elements_; i++)
		{
			SAFE_DELETE(random_[i]);
		}
	}

	void initialize(const unsigned int _num_threads = std::thread::hardware_concurrency())
	{
		deleteAllThreads();

		num_threads_ = _num_threads;
		
		thread_list_.initialize(num_threads_, nullptr);

		start_ix_1D_.initialize(num_threads_);
		end_ix_1D_.initialize(num_threads_);

		if (sync_value_float_ != nullptr) delete[] sync_value_float_;
		sync_value_float_ = new float[num_threads_];

		if (sync_value_int_ != nullptr) delete[] sync_value_int_;
		sync_value_int_ = new int[num_threads_];

		if (sync_value_double_ != nullptr) delete[] sync_value_double_;
		sync_value_double_ = new double[num_threads_];

		random_.initialize(num_threads_);

		for (int i = 0; i < random_.num_elements_; i++)
			random_[i] = new RandomNumberGenerator(i);
	}

	RandomNumberGenerator& getRandom(const int thread_id)
	{
		return *random_[thread_id];
	}

	void lock()
	{
		mtx_lock_.lock();
	}

	void unlock()
	{
		mtx_lock_.unlock();
	}

	void sync()
	{
		std::unique_lock<std::mutex> lock(mtx_lock_);

		num_of_synced_threads_++;

		if (num_of_synced_threads_ == num_threads_)
		{
			num_of_synced_threads_ = 0;

			cond_.notify_all();
		}
		else
			cond_.wait(lock);
	}

	void wait()	// called by the threads to sleep
	{
		std::unique_lock<std::mutex> lock(mtx_lock_);
		
		num_waiting_threads_++;

		cond_.wait(lock);

		num_waiting_threads_--;
	}

	bool isAllWaiting()	// not tested yet
	{
		if (num_waiting_threads_ == num_threads_) return true;
		else return false;
	}

	void notifyAll()	// called by a different thread outside of this class
	{
		cond_.notify_all();		// be careful not to notify the threads in sync.
	}

	void deleteAllThreads()
	{
		for (unsigned int i = 0; i < num_threads_; i++) SAFE_DELETE(thread_list_[i]);
	}

	void joinAll()
	{
		for (unsigned int i = 0; i < num_threads_; i++) thread_list_[i]->join();

		deleteAllThreads();
	}

	void splitRange(const int& k_start, const int& k_end)
	{
		const int k_res = k_end - k_start + 1;
		const int quotient = k_res / num_threads_;
		const int remainder = k_res % num_threads_;

		int k_start_p = k_start;

		for (int i = 0; i < (int)num_threads_; i++)
		{
			int k_depth = i < remainder ? (quotient + 1) : quotient;

			start_ix_1D_[i] = k_start_p;
			end_ix_1D_[i] = k_start_p + k_depth - 1;

			k_start_p += k_depth;
		}
	}

	Vector2D<int> getParallelRange(const int thread_id, const int k_start, const int k_end)
	{
		if (thread_id == 0) splitRange(k_start, k_end);

		sync();

		Vector2D<int> range(start_ix_1D_[thread_id], end_ix_1D_[thread_id]);

		return range;
	}

	template <class A0, class ...Args>
	void printf_once(const A0& a0, const Args& ...args)
	{
		std::call_once(o_flag, printf, a0, args...);
	}

	double syncSum(const int& thread_id, double value)	//Note: the variable 'value' should not be accessed by many threads.
	{
		sync_value_double_[thread_id] = value;

		sync();

		double sync_value_double_temp = sync_value_double_[0];
		for (int i = 1; i < (int)num_threads_; i++) sync_value_double_temp += sync_value_double_[i];

		return sync_value_double_temp;
	}

	/*
	double syncSum(const int& thread_id, double value)	//Note: the variable 'value' should not be accessed by many threads.
	{
		sync_value_double_[thread_id] = value;

		sync();

		if (thread_id == 0)
		{
			sync_value_double_temp_ = sync_value_double_[0];// to remove one assignment
			for (int i = 1; i < (int)num_threads_; i++) sync_value_double_temp_ += sync_value_double_[i];
		}
		//Note: we can let head thread perform the adding operation in order to save memory bandwidth.

		sync();

		return sync_value_double_temp_;
	}*/

	float syncSum(const int& thread_id, float value)	//Note: the variable 'value' should not be accessed by many threads.
	{
		sync_value_float_[thread_id] = value;

		sync();

		if (thread_id == 0)
		{
			sync_value_float_temp_ = sync_value_float_[0];// to remove one assignment
			for (int i = 1; i < (int)num_threads_; i++) sync_value_float_temp_ += sync_value_float_[i];
		}
		//Note: we can let head thread perform the adding operation in order to save memory bandwidth.

		sync();

		return sync_value_float_temp_;
	}

	int syncSum(const int& thread_id, int value)
	{
		sync();

		sync_value_int_[thread_id] = value;

		sync();

		value = sync_value_int_[0];// to remove one assignment
		for (int i = 1; i < (int)num_threads_; i++) value += sync_value_int_[i];

		sync();

		return value;
	}

	void syncMax(const int& thread_id, float& value)
	{
		sync();

		sync_value_float_[thread_id] = value;

		sync();

		for (unsigned int i = 0; i < num_threads_; i++) if (value < sync_value_float_[i]) value = sync_value_float_[i];

		sync();
	}

	void syncMax(const int& thread_id, int& value)
	{
		sync();

		sync_value_int_[thread_id] = value;

		sync();

		for (unsigned int i = 0; i < num_threads_; i++) if (value < sync_value_int_[i]) value = sync_value_int_[i];

		sync();
	}

	void syncMin(const int& thread_id, float& value)
	{
		sync();

		sync_value_float_[thread_id] = value;

		sync();

		for (unsigned int i = 0; i < num_threads_; i++) if (value > sync_value_float_[i]) value = sync_value_float_[i];

		sync();
	}

	template <class F, class ...Args>
	void run(F f, const Args& ...args)
	{
		for (unsigned int thread_id = 0; thread_id < num_threads_; thread_id++)
			thread_list_[thread_id] = new std::thread(f, args...);
	}

	// works with class member functions only due to this pointer
	template <class F, class Athis, class ...Aetc>
	void runWithID(F f, Athis at, const Aetc& ...args)
	{
		for (unsigned int thread_id = 0; thread_id < num_threads_; thread_id++)
			thread_list_[thread_id] = new std::thread(f, at, this, thread_id, args...);
//			thread_list_[thread_id] = new std::thread(f, at, *this, thread_id, args...);
	}

	template<class T>
	void copyArray(const int thread_id, const Array1D<T>& from_array, Array1D<T>& to_array)
	{
		if (thread_id == 0) to_array.initialize(from_array.num_elements_);

		sync();

		const TV2_INT range = getParallelRange(thread_id, 0, from_array.num_elements_ - 1);

		memcpy(&to_array.values_[range.t_min_], &from_array.values_[range.t_min_], (range.t_max_ - range.t_min_ + 1)*sizeof(T));

		sync();
	}

	template<class TT>
	void copyArray(const int thread_id, const int& num_values, const TT* from_array, TT* to_array)
	{
		if (thread_id == 0)
		{
			SAFE_DELETE_ARRAY(to_array);
			to_array = new TT [num_values];
		}

		sync();

		const TV2_INT range = getParallelRange(thread_id, 0, num_values - 1);

		memcpy(&to_array[range.t_min_], &from_array[range.t_min_], (range.t_max_ - range.t_min_ + 1)*sizeof(TT));

		sync();
	}

	template<class TT>
	void setArray(const int thread_id, const int& num_values, const TT value, TT* to_array)
	{
		const TV2_INT range = getParallelRange(thread_id, 0, num_values - 1);

		memset(&to_array[range.t_min_], (int)value, (range.t_max_ - range.t_min_ + 1)*sizeof(TT));

		sync();
	}

	// ex) thread 0's res = 3, 1's = 5, 2's 2 -> thread 0's range = [0, 2], 1's range = [3, 7], 2's range = [8, 9]
	Vector2D<int> getIncrementalRange(const int& thread_id, const int& res)
	{
		sync_value_int_[thread_id] = res;

		sync();

		int i_start = 0;
		for (int i = 0; i <= thread_id - 1; i++)
		{
			i_start += sync_value_int_[i];
		}

		int i_end = i_start + res - 1;

		sync();

		return Vector2D<int>(i_start, i_end);
	}
};

typedef MultiThreading MT;

#define ONE_THREAD_WORK(expression) if(thread_id == 0){expression;};mt->sync();

#define BEGIN_BOX_ITERATION_ARR_3D_SYNC(_grid)	{const TV2_INT k_range = mt->getParallelRange(thread_id, _grid.k_start_, _grid.k_end_);\
												int i_start(_grid.i_start_),i_end(_grid.i_end_),j_start(_grid.j_start_),j_end(_grid.j_end_),k_start(k_range.t_min_),k_end(k_range.t_max_);\
												int i, j, k, k_itr, k_size=k_end-k_start+1, arr_ix; mt->syncMax(thread_id, k_size);\
												for(k=k_start,k_itr=0;k_itr<k_size;++k_itr,++k){mt->sync();if(k > k_end)continue;\
												for(j=j_start; j<=j_end; ++j) for(i=i_start,arr_ix=_grid.get1DIndex(i,j,k);i<=i_end;++i,++arr_ix)
#define END_GRID_ITERATION_Z_SYNC				} mt->sync(); }

#define BEGIN_1D_ITERATION(num)					{const TV2_INT k_range = mt->getParallelRange(thread_id, 0, num-1);\
												const int _p_start(k_range.t_min_), _p_end(k_range.t_max_); \
												for(int p = _p_start; p <= _p_end; p++)
#define END_1D_ITERATION						mt->sync();}

#define BEGIN_ONE_THREAD_WORK					if(thread_id == 0)
#define END_ONE_THREAD_WORK						mt->sync();

#define BEGIN_PARALLEL_FOR(_itr, _start, _end)	{const Vector2D<int> range = mt->getParallelRange(thread_id, _start, _end);\
												const int _itr_start(range.t_min_), _itr_end(range.t_max_);\
												for(int _itr=_itr_start;_itr<=_itr_end;++_itr)
#define END_PARALLEL_FOR						mt->sync();}


#define BEGIN_GRID_ITERATION_2D(_grid)		{const TV2_INT j_range = mt->getParallelRange(thread_id, _grid.j_start_, _grid.j_end_);\
												GridUniform2D &grid_itr(_grid);	\
												int i, j;						\
												const int j_start = j_range.t_min_, j_end = j_range.t_max_, i_start = grid_itr.i_start_, i_end = grid_itr.i_end_;	\
												for(j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)
#define END_GRID_ITERATION_2D					mt->sync();}

#define BEGIN_ARRAY_ITERATION_2D(_arr)		{const TV2_INT j_range = mt->getParallelRange(thread_id, _arr.j_start_, _arr.j_end_);\
											int i, j;						\
											const int j_start = j_range.t_min_, j_end = j_range.t_max_, i_start = _arr.i_start_, i_end = _arr.i_end_;	\
											for(j = j_start; j <= j_end; ++j) for(i = i_start; i <= i_end; ++i)
#define END_ARRAY_ITERATION_2D				mt->sync();}

#define BEGIN_GRID_ITERATION_3D(_grid)			{const TV2_INT k_range = mt->getParallelRange(thread_id, _grid.k_start_, _grid.k_end_);\
												int i_start(_grid.i_start_),i_end(_grid.i_end_),j_start(_grid.j_start_),j_end(_grid.j_end_),k_start(k_range.t_min_),k_end(k_range.t_max_);\
												int i, j, k, arr_ix;\
												for(k=k_start;k<=k_end;++k)for(j=j_start; j<=j_end; ++j) for(i=i_start,arr_ix=_grid.get1DIndex(i,j,k);i<=i_end;++i,++arr_ix)
#define END_GRID_ITERATION_3D					mt->sync(); }