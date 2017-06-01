#pragma once

#include "GL_OBJECT.h"
#include "Framework/FluidObject.h"

class GL_FluidObject : public GL_OBJECT
{
public:
// 	GL_Fluid()
// 	{}

	GL_FluidObject(Object* object_input)
		: GL_OBJECT(object_input)
	{}

	void initialize(const FluidObject& fluid)
	{
		{GL_ELEMENT *gl_element = new GL_ELEMENT;
		gl_element->name_ = std::string("flip_particles");
		gl_element->setPrimitiveType(GL_POINTS);				// GL_POINTS doesn't work!
		gl_element->point_size_ = 3.0f;
// 		gl_element->setPositionArray(&SimulationWorld::getInstance().fluid_display_data_.flip_particles_pos_);
// 		gl_element->setColorArray(&SimulationWorld::getInstance().fluid_display_data_.flip_particles_color_);
		gl_element_list_.push_back(gl_element); }

		{GL_ELEMENT *gl_element = new GL_ELEMENT;
		gl_element->name_ = std::string("fluid_aabb");
		gl_element->setPrimitiveType(GL_LINES);				// GL_POINTS doesn't work!
		gl_element->line_width_ = 1.0f;
//		gl_element->setPositionArray(&SimulationWorld::getInstance().fluid_display_data_.aabb_vertices_);
		gl_element_list_.push_back(gl_element); }

		{GL_ELEMENT *gl_element = new GL_ELEMENT;
		gl_element->name_ = std::string("water_surface");
		gl_element->setPrimitiveType(GL_TRIANGLES);				// GL_POINTS doesn't work!
		gl_element->setFacingMode(GL_FRONT_AND_BACK);
		gl_element->setRasterizeMode(GL_FILL);
// 		gl_element->setPositionArray(&SimulationWorld::getInstance().fluid_display_data_.water_surface_vertex_position_);
// 		gl_element->setNormalArray(&SimulationWorld::getInstance().fluid_display_data_.water_surface_vertex_normal_);
		gl_element_list_.push_back(gl_element); }
	}
};