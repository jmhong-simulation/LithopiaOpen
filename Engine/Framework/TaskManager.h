// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"

#include "Parallelism/MultiThreading.h"
#include "DataDepot.h"
#include "Geometry/BOX_3D.h"
#include "Utilities/ScriptParameterList.h"

class TaskManager
{
public:
	TaskManager(){};		// temporary
	TaskManager(DataDepot* data);
	~TaskManager();

public:
	bool write_files_;
	bool run_repeat_;

	std::atomic<bool> run_flag_;		//TODO: need to be atomic if this is controlled by one controller class?
	std::atomic<bool> exit_;

	MultiThreading manager_thread_;		// manager_thread controls work_threads
	MultiThreading work_threads_;		// these threads

public:
	DataDepot	*data_depot_;

public:	
	void pause();
	void resume();
	void changeRunStatus();

	void runManagerThread();
	void updateParallel();
	void updateDataDepotParallel();
	
	virtual void updateOneStep(MT* mt, const int thread_id) {std::cout << "TaskManager::updateOneStep(MT* mt, const int thread_id)" << std::endl;}
	virtual void updateDataDepot(MT* mt, const int thread_id, DataDepot* data) { std::cout << "TaskManager::updateDataDepot(MT* mt, const int thread_id)" << std::endl; }
	virtual void writeFiles(){ std::cout << "TaskManager::writeFiles()" << std::endl; data_depot_->unlock(); }
	virtual BOX_3D<T> getAABB(){ std::cout << "TaskManager::getAABB()" << std::endl; return BOX_3D<T>(0, 0, 0, 1, 1, 1); }
	virtual void receiveMessage(ScriptParameterList& message_parser){ std::cout << "TaskManager::receiveMessage(const ScriptParameterList& message_parser)" << std::endl; }
};

