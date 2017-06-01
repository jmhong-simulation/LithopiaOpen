#pragma once

#include "GL_OBJECT.h"
#include "GL_VIEW.h"
#include "GL_SHADER_MANAGER.h"

#include "Geometry/SLICER.h"
#include "DataStructure/SinglyLinkedList.h"

class GL_WORLD
{
private:
	static GL_WORLD *gl_world_;

	GL_WORLD()
		: axis_element_(nullptr), draw_axis_(true), gl_sim_obj_(nullptr), capture_(false)
	{
		InitializeCamera();
	}

	~GL_WORLD();

public:
	bool draw_axis_;
	bool capture_;

	GL_VIEW camera_;
	GL_ELEMENT *axis_element_;
	GL_OBJECT *gl_sim_obj_;
	SinglyLinkedList<GL_OBJECT*> object_list_;

	GL_SHADER_MANAGER shader_programs_;

public:
	static GL_WORLD& GetInstance();
	static GL_WORLD* GetPointer();

public:
	void InitializeObjectList();
	void InitializeRenderOptions();
	void InitializeCamera();
	void InitializeAxis();

	void Draw();
	void DrawAxis();
	
	void ResetObjectList();

	GL_OBJECT* GL_WORLD::GetObjectPtr(const int& ix);
};