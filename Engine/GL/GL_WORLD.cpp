#include "GL_WORLD.h"
#include "MAIN_WINDOW.h"
#include "Simulation/SoundTo3D.h"
//#include "Simulation/SimulationWorld.h"

GL_WORLD* GL_WORLD::gl_world_ = nullptr;

GL_WORLD& GL_WORLD::GetInstance()
{
	if (gl_world_ == nullptr) gl_world_ = new GL_WORLD;

	return *gl_world_;
}

GL_WORLD* GL_WORLD::GetPointer()
{
	if (gl_world_ == nullptr) gl_world_ = new GL_WORLD;

	return gl_world_;
}

GL_WORLD::~GL_WORLD()
{
	ResetObjectList();

	SAFE_DELETE(axis_element_);
	SAFE_DELETE(gl_sim_obj_);
}

void GL_WORLD::InitializeCamera()
{
	camera_.SetPanScale(0.01f);
	camera_.SetDollyStartPosition(-5.0f);
	camera_.SetDollyScale(0.001f);
	camera_.SetTrackballScale(100.0f);
	camera_.SetCenterOfRotation(glm::vec3(0, 0, 0));
}

void GL_WORLD::ResetObjectList()
{
	for (object_list_.begin(); object_list_.valid(); object_list_.next())
	{
		delete object_list_.getItr();
	}	

	object_list_.reset();
}

void GL_WORLD::InitializeRenderOptions()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glEnable(GL_DEPTH_TEST);		// Enable depth buffer
	glEnable(GL_CULL_FACE);			// Enable back face culling
	glFrontFace(GL_CCW);
//	glDisable(GL_CULL_FACE);
}

void GL_WORLD::Draw()
{
	if (draw_axis_) DrawAxis();

// 	if (gl_sim_obj_)
// 	{
//  		if (SimulationWorld::getInstance().updated_ == true)
//  		{
// 			while (SimulationWorld::getInstance().data_lock_ == true) std::this_thread::sleep_for(std::chrono::milliseconds(1));		// wait if simulation data is being copied
// 
// 			SimulationWorld::getInstance().data_lock_ = true;
// 
// 			gl_sim_obj_->update();
// 
// 			SimulationWorld::getInstance().updated_ = false;
// 
// 			SimulationWorld::getInstance().data_lock_ = false;	
// 
// 			capture_ = true;
// 		}
// 
// 		gl_sim_obj_->Draw(camera_);
// 	}

	//TODO: disable picking for performance reason
//	glInitNames();

	for (object_list_.begin(); object_list_.valid(); object_list_.next())
	{
		GL_OBJECT *gl_obj = object_list_.getItr();

//		if (gl_obj->object_->data_depot_.updated_ == true) capture_ = true;

		if (MAIN_WINDOW::getInstance().audio_input_test_ != nullptr)
		{
			if (gl_obj->object_->name_ == std::string("SOUND_TO_3D"))
			{
				SoundTo3D *task_manager = dynamic_cast<SoundTo3D*>(gl_obj->object_->task_manager_);

				const int num_samples = MAIN_WINDOW::getInstance().audio_input_test_->m_audioInfo->audio_samples_.num_elements_;

				if (task_manager->audio_spectrum_.num_elements_ != num_samples)
				{
					task_manager->audio_spectrum_.initialize(num_samples);					
				}

//				task_manager->audio_spectrum_.copyFrom(MAIN_WINDOW::getInstance().audio_input_test_->m_audioInfo->audio_samples_);
				task_manager->updateHeight(MAIN_WINDOW::getInstance().audio_input_test_->m_audioInfo->audio_samples_);

				task_manager->updateDataDepotParallel();
				task_manager->data_depot_->updated_ = true;
			}
		}

		gl_obj->update();
	}

	for (object_list_.begin(); object_list_.valid(); object_list_.next())
	{
//		std::cout << "Drawing GL_Object " << object_list_.getItr()->object_->name_ << std::endl;

		object_list_.getItr()->Draw(camera_);
	}

	// Note that you can push a dummy name (some unused value) right at the start, and afterwards just use glLoadName, instead of Pushing and Popping names. However it is faster to disregard objects that don't have a name than to check if the name is the dummy name. On the other hand it is probably faster to call glLoadName than it is to Pop followed by Push.
	// http://www.lighthouse3d.com/opengl/picking/
}

void GL_WORLD::DrawAxis()
{
	axis_element_->program_->bind();
	axis_element_->setUniformMVPMatrix(glm::translate(glm::vec3(-0.87f, -0.8f, 0.0f)) * camera_.GetScreenViewMatrix());
	axis_element_->drawArrays();
	axis_element_->program_->release();
}

void GL_WORLD::InitializeObjectList()
{
	InitializeAxis();
}

void GL_WORLD::InitializeAxis()
{
	axis_element_ = new GL_ELEMENT;
	axis_element_->initializeQTGL();
	axis_element_->setProgram(&shader_programs_.program_pc_);

	Array1D<TV> positions(6);
	{
		int ix = 0;
		positions[ix++] = TV(0.0f, 0.0f, 0.0f) * 0.1f;
		positions[ix++] = TV(1.0f, 0.0f, 0.0f) * 0.1f;
		positions[ix++] = TV(0.0f, 0.0f, 0.0f) * 0.1f;
		positions[ix++] = TV(0.0f, 1.0f, 0.0f) * 0.1f;
		positions[ix++] = TV(0.0f, 0.0f, 0.0f) * 0.1f;
		positions[ix++] = TV(0.0f, 0.0f, 1.0f) * 0.1f;
	}

	Array1D<glm::vec4> colors(6);
	{
		int ix = 0;
		colors[ix++] = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
		colors[ix++] = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
		colors[ix++] = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
		colors[ix++] = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
		colors[ix++] = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);
		colors[ix++] = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);
	}

	axis_element_->begin(GL_LINES, USE_POSITION | USE_COLOR);
	axis_element_->bindArray(USE_POSITION, positions);
	axis_element_->bindArray(USE_COLOR, colors);
	axis_element_->end();
}

GL_OBJECT* GL_WORLD::GetObjectPtr(const int& ix)
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