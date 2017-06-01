// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "Object.h"
#include "DataStructure/SinglyLinkedList.h"

class World
{
private:
	static World *world_;

	World()
		: sim_obj_(nullptr)
	{}

	~World()
	{
		SAFE_DELETE(sim_obj_);
		ResetObjectList();
	}

public:
	static World& GetInstance();
	static World* GetPointer();

public:
	Object *sim_obj_;	// TODO remove

	SinglyLinkedList<Object*> object_list_;

	Object* GetObjectPtr(const int& ix);
	Object* getObjectPtr(const std::string& name);

	void ResetObjectList();

	void update()
	{

	}
};