#pragma once

#include <glm/glm.hpp>

class GL_MATERIAL_PARAMETERS
{
public:
	glm::vec4 Ka_;
	glm::vec4 Kd_;
	glm::vec4 Ks_;
	glm::vec4 Ke_;
	float Se_;

	GL_MATERIAL_PARAMETERS()
	{
		Ka_ = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
		Kd_ = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
		Ks_ = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
		Ke_ = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f);

		Se_ = 20.0f;
	}

	~GL_MATERIAL_PARAMETERS()
	{}
};