// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <atomic>
#include <iostream>
#include <glm/glm.hpp>
#include "Parallelism/MultiThreading.h"
#include "DataStructure/SinglyLinkedList.h"
#include "DataStructure/LinkedArray.h"

class ColoredParticlesData
{
public:
	ColoredParticlesData()
		: name_(std::string("noname")), point_size_(1.0f)
	{}

	std::string name_;
	float point_size_;

	Array1D<TV>			position_;
	Array1D<glm::vec4>	color_;
};

class LinesData
{
public:
	LinesData()
		: name_(std::string("noname")), line_width_(1.0f), color_(glm::vec4(0,0,0,0))
	{}

	float line_width_;
	glm::vec4 color_;
	std::string name_;

	LinkedArray<TV> vertices_;
};

class PhongTrianglesData
{
public:
	PhongTrianglesData()
		: name_(std::string("noname")), rasterize_mode_(0x1B02)	// 0x1B02 is GL_FILL, 0x1B01 is GL_LINE
	{}

	std::string name_;

	int rasterize_mode_;

	LinkedArray<TV>	 positions_;
	LinkedArray<TV>  normals_;
};

class DataDepot
{
public:
	std::atomic<bool> is_locked_;		// check if display data is locked
	std::atomic<bool> updated_;

	SinglyLinkedList<ColoredParticlesData*> colored_particles_list_;
	SinglyLinkedList<LinesData*>			lines_list_;
	SinglyLinkedList<PhongTrianglesData*>   phong_triangles_list_;

	DataDepot()
		: is_locked_(false), updated_(false)
	{}

	~DataDepot()
	{
		reset();
	}

//	virtual void update(MT* mt, const int thread_id){ std::cout << "DataDepot::void update(MT* mt, const int thread_id)" << std::endl; }
	virtual void writeFiles(){ std::cout << "DataDepot::writeFiles()" << std::endl; }

	void lock()
	{
		while (is_locked_ == true) std::this_thread::sleep_for(std::chrono::milliseconds(1));

		// try compare_exchange_strong
		
		is_locked_ = true;		//TODO: wait when is_locked_ is true
	}
	
	void unlock()
	{
		is_locked_ = false;
	}

	void reset()
	{
		colored_particles_list_.resetPointers();
		lines_list_.resetPointers();
		phong_triangles_list_.resetPointers();
	}
};