// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "GENERIC_DEFINITIONS.h"
#include "DataDepot.h"
#include "TaskManager.h"
#include "DataStructure/GridAdaptive3D.h"
#include "DataStructure/GridUniform3D.h"
#include "DataStructure/SinglyLinkedList.h"
#include "Geometry/StaticTriangularSurface.h"
#include "Geometry/SLICER.h"
#include "Geometry/BOX_3D.h"
#include "Geometry/DynamicParticles.h"
#include "GL/GL_MATERIAL_PARAMETERS.h"
#include "Utilities/ScriptParameterList.h"

class Object	// Object receives all messages (commands) from UI and distributes them to task_object, data_depot, task_manager, etc. properly.
{
public:
	std::string name_;

	TaskManager *task_manager_;
	DataDepot	data_depot_;

	// rendering options (what to draw) so that GL_Object passively update rendering data. Object determines what to draw.
	// UI options

	BOX_3D<T>	aabb_;					// TODO: move AABB and m_matrix to UI_Object because they can be determined by DataDepot
	glm::mat4	m_matrix_;				// for rendering, move to UI_Object (?)

	StaticTriangularSurface surface_;
	SLICER slicer_;
	GridAdaptive3D amr_grid_;
	//DynamicParticles particles_; use SimulationWorld::GetInstance().particles_. Don't Simulation include Object.

//	LEVELSET volume_;

	// material properties
	GL_MATERIAL_PARAMETERS material_;	//TODO: move to UI_Object

	// scene graph
	Object *parent_;
	SinglyLinkedList<Object*> children_;

	Object();
	Object(const std::string& message);

	void updateAABB(const BOX_3D<T>& _aabb)
	{
		aabb_ = _aabb;

		// update m_matrix from aabb_
		float scale = MAX3(aabb_.x_max_ - aabb_.x_min_, aabb_.y_max_ - aabb_.y_min_, aabb_.z_max_ - aabb_.z_min_);
		if (scale != 0.0f) scale = 1.0f / scale;

		glm::vec3 scalevec(scale, scale, scale);
		const Vector3D<float> aabb_center = aabb_.GetCenter();
		const glm::vec3 center(aabb_center.x_, aabb_center.y_, aabb_center.z_);

		m_matrix_ = glm::scale(scalevec) * glm::translate(-center);
	}

	void receiveMessage(const std::string& message)
	{
		ScriptParameterList message_parser;

		message_parser.initializeFromMessage(message);

		task_manager_->receiveMessage(message_parser);
	}
};
