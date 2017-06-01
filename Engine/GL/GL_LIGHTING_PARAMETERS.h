#pragma once

#include <glm/glm.hpp>

class GL_LIGHT_PARAMETERS
{
public:
	glm::vec3 spot_direction_;
	int spot_exponent_;
	int spot_cutoff_;

	float Kc_;
	float Kl_;
	float Kq_;

	glm::vec4 light_pos_;
	glm::vec4 light_Ka_;
	glm::vec4 light_Kd_;
	glm::vec4 light_Ks_;

	GL_LIGHT_PARAMETERS()
	{
		spot_direction_ = glm::vec3(1.0f, -1.0f, -1.0f);
		spot_exponent_ = 10;
		spot_cutoff_ = 180;

		Kc_ = 1.0f;
		Kl_ = 0.0f;
		Kq_ = 0.0f;

		light_pos_ = glm::vec4(0.0f, 2.0f, -3.0f, 1.0f);
		light_Ka_ = glm::vec4(0.1f, 0.1f, 0.1f, 1.0f);
		light_Kd_ = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
		light_Ks_ = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
	}
};

class GL_LIGHTING_PARAMETERS
{
public:
	bool use_lighting_;
	bool use_light_[8];
	bool normalize_normals_;

	glm::vec4 lmKa_;

	GL_LIGHT_PARAMETERS light_params_[8];

	GL_LIGHTING_PARAMETERS()
	{
		use_lighting_ = true;
		use_light_[0] = true;
		normalize_normals_ = true;

		lmKa_ = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f);
	}
};
