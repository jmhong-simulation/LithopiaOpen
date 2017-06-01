// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include <chrono>

#include "TaskManager.h"
#include "DataStructure/LinkedArray.h"
#include "Parallelism/MultiThreading.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/Timer.h"

TaskManager::TaskManager(DataDepot* data)
	: run_flag_(false), exit_(false), write_files_(true), data_depot_(data), run_repeat_(false)
{
	data_depot_ = data;

	manager_thread_.initialize(1);		// number of manager thread should be 1.
	work_threads_.initialize();

	runManagerThread();		// run manager thread when this instance is created. it will sleep until resume is called.
}

TaskManager::~TaskManager()
{
	exit_ = true;

	// don't delete task_object or data_depot here
}

void TaskManager::pause()
{
	run_flag_ = false;

	printf("Pause simulation\n");
}

void TaskManager::resume()
{
	if (run_flag_ == false)
	{
		run_flag_ = true;

		manager_thread_.notifyAll();
	}

	printf("Resume simulation\n");
}

void TaskManager::changeRunStatus()
{
	if (run_flag_ == true) pause();
	else resume();
}

void TaskManager::runManagerThread()
{
	// starts manager thread
	manager_thread_.run(&TaskManager::updateParallel, this);
}

void TaskManager::updateParallel()
{
	Timer frame_timer;

	while (true)
	{	
		if (run_flag_ == false)	
		{
			manager_thread_.wait();		// start waiting when simulation is paused by UI. resume() call by UI can resume the simulation
			continue;					//TODO: necessary?
		}

//		frame_timer.start();

		work_threads_.runWithID(&TaskManager::updateOneStep, this);		// adavanceOneFrame takes care of simulation
		work_threads_.joinAll();

//		frame_timer.end();

		// make copies of simulation data for rendering and writing
		// stop simulation during copying
		// TODO: use atomic operator

		data_depot_->lock();

		work_threads_.runWithID(&TaskManager::updateDataDepot, this, data_depot_);
		work_threads_.joinAll();

		if (write_files_ == true)
			std::thread *file_write_thread = new std::thread(&TaskManager::writeFiles, this); //TODO: check if this pointer needs to be deleted after writing is finished
		else data_depot_->unlock();	// writeFiles makes data_lock_ = false at the end

		data_depot_->updated_ = true;			// UI knows whether data need to be updated or not.

		if (run_repeat_ == false) run_flag_ = false;
		if (exit_ == true) break;
	}
}

/*
void TaskManager::writeFiles()
{
	// write data files here
//	data_lock_ = false;	// free data lock when file writing is done.
	data_depot_->writeFiles();

	data_depot_->unlock();
}
*/

/*
void TaskManager::updateOneStep(MT* mt, const int thread_id)
{
	task_object_->update(mt, thread_id);
}
*/

/*
void TaskManager::updateDataDepot(MT* mt, const int thread_id)
{
//	fluid_display_data_.copyData(mt, thread_id, fluids_);

//	data_depot_->update(mt, thread_id);

//	task_object_->updateIOData(mt, thread_id, *data_depot_);
}
*/

void TaskManager::updateDataDepotParallel()
{
	work_threads_.runWithID(&TaskManager::updateDataDepot, this, data_depot_);
	work_threads_.joinAll();
}