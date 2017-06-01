// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "Object.h"
#include "Simulation/FluidSimulation.h"
#include "Simulation/VolumeEditor.h"
#include "Simulation/LevelsetEditor2D.h"
#include "Simulation/CookieCutterMaker.h"
#include "Simulation/ShrinkWrappingExample2D.h"
#include "Simulation/Photo2_5D.h"
#include "Simulation/LithophaneMaker.h"
#include "Simulation/SoundTo3D.h"

Object::Object(const std::string& message)
	: Object()
{
	ScriptParameterList message_parser;

	message_parser.initializeFromMessage(message);

	message_parser.printAll();

	if (message_parser.getStringValue("type") == "FLUID_SIMULATION")
	{
		FluidSimulation *new_task = new FluidSimulation(&data_depot_);
		new_task->initialize();		//TODO: add script name

		name_ = std::string("fluid_simulation_object");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;
		
		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();

	}
	else if (message_parser.getStringValue("type") == "VOLUME_EDITOR")
	{
		VolumeEditor *new_task = new VolumeEditor(&data_depot_);
		new_task->initializeFromUltrasonography("./Ultrasonography/3D_0002.Raw");		//TODO: add script name
		new_task->write_files_ = false;

		name_ = std::string("volume_editor");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();
	}
	else if (message_parser.getStringValue("type") == "LEVELSET_EDITOR_2D")
	{
		LevelsetEditor2D *new_task = new LevelsetEditor2D(&data_depot_);
		new_task->initialize();
		new_task->write_files_ = false;

		name_ = std::string("levelset_editor_2d");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();
	}
	else if (message_parser.getStringValue("type") == "COOKIE_CUTTER_MAKER")
	{
		CookieCutterMaker *new_task = new CookieCutterMaker(&data_depot_);
		new_task->initializeFromScript(message_parser.getStringValue("script_filename").c_str());
		new_task->write_files_ = false;

		name_ = std::string("COOKIE_CUTTER_MAKER");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();

		updateAABB(task_manager_->getAABB());
	}
	else if (message_parser.getStringValue("type") == "SHRINK_WRAPPING_EXAMPLE_2D")
	{
		ShrinkWrappingExample2D *new_task = new ShrinkWrappingExample2D(&data_depot_);
		new_task->initializeFromScript(message_parser.getStringValue("script_filename").c_str());
		new_task->write_files_ = false;

		name_ = std::string("shrink_wrapping_exmaple_2d");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();

		updateAABB(task_manager_->getAABB());
	}
	else if (message_parser.getStringValue("type") == "PHOTO_2.5D")
	{
		Photo2_5D *new_task = new Photo2_5D(&data_depot_);
		new_task->initializeFromScript(message_parser.getStringValue("script_filename").c_str());
		new_task->write_files_ = false;

		name_ = std::string("photo_2_5D");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();

		updateAABB(task_manager_->getAABB());
	}
	else if (message_parser.getStringValue("type") == "LITHOPHANE_MAKER")
	{
		LithophaneMaker *new_task = new LithophaneMaker(&data_depot_);
		new_task->initializeFromScript(message_parser.getStringValue("script_filename").c_str());
		new_task->write_files_ = false;

		name_ = std::string("LITHOPHANE_MAKER");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = true;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();

		updateAABB(task_manager_->getAABB());
	}
	else if (message_parser.getStringValue("type") == "SOUND_TO_3D")
	{
		SoundTo3D *new_task = new SoundTo3D(&data_depot_);
		new_task->initializeFromScript(message_parser.getStringValue("script_filename").c_str());
		new_task->write_files_ = false;

		name_ = std::string("SOUND_TO_3D");
		task_manager_ = new_task;
		task_manager_->run_repeat_ = false;

		// update first data
		data_depot_.lock();
		task_manager_->updateDataDepotParallel();
		data_depot_.updated_ = true;
		data_depot_.unlock();

		updateAABB(task_manager_->getAABB());
	}
	else
	{
		std::cout << "Object::Object(const std::string& message) Unknown object type" << std::endl;
		exit(1);
	}
}

Object::Object()
	: parent_(nullptr), name_(std::string("noname")), task_manager_(nullptr), aabb_(BOX_3D<T>(0, 0, 0, 1, 1, 1))
{
	updateAABB(aabb_);
}
