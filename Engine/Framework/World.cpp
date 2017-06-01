// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "World.h"

World* World::world_ = nullptr;

World& World::GetInstance()
{
	if (world_ == nullptr)
		world_ = new World;

	return *world_;
}

World* World::GetPointer()
{
	if (world_ == nullptr)
		world_ = new World;

	return world_;
}

Object* World::GetObjectPtr(const int& ix)
{
	object_list_.begin();

	int i = 0;
	while (true)
	{
		if (object_list_.valid() == false) return nullptr;

		if (i == ix) return object_list_.itr_->value_;

		object_list_.next();
		i++;
	}
}

Object* World::getObjectPtr(const std::string& name)
{
	object_list_.begin();

	int i = 0;
	while (true)
	{
		if (object_list_.valid() == false) return nullptr;

		if (object_list_.itr_->value_->name_ == name) return object_list_.itr_->value_;

		object_list_.next();
		i++;
	}
}

void World::ResetObjectList()
{
	for (object_list_.begin(); object_list_.valid();)
	{
		Object *temp = object_list_.getItr();

		object_list_.next();

		SAFE_DELETE(temp);
	}

	object_list_.reset();
}