// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "LithophaneMaker.h"
#include "GL/GL_TOOLS.h"
#include "Geometry/LITHOPHANE.h"

LithophaneMaker::LithophaneMaker(DataDepot* data_)
	: TaskManager(data_)
{}

void LithophaneMaker::initializeFromScript(const std::string script_filename)
{
	ScriptParameterList par_reader;

	par_reader.initializeFromFile(script_filename.c_str());

	using namespace std;

	DEFINE_PARAMETER(string, model_type);
	DEFINE_PARAMETER(int, target_width);
	DEFINE_PARAMETER(int, target_height);

	DEFINE_PARAMETER(string, input_path);
	DEFINE_PARAMETER(string, input_filename);
	DEFINE_PARAMETER_WITH_DEFAULT(string, output_path, input_path);
	DEFINE_PARAMETER(string, output_filename_prefix);

	output_path_ = output_path;
	output_filename_prefix_ = output_filename_prefix;

	DEFINE_PARAMETER(float, inner_radius);
	DEFINE_PARAMETER(float, base_thickness);
	DEFINE_PARAMETER(float, front_thickness);
	DEFINE_PARAMETER(float, height_scale);
	DEFINE_PARAMETER(int, num_smoothing);

	DEFINE_PARAMETER(bool, write_stl);
	DEFINE_PARAMETER(bool, write_ctm);
	DEFINE_PARAMETER(bool, write_outer_base_ctm);		// this contains texture uv

	DEFINE_PARAMETER(bool, write_texture);
	DEFINE_PARAMETER(bool, write_bumpmap);

	IMAGE_2D height_map;	// current write as BMP24. Use gray le BMP for networking
	IMAGE_2D texture;

	if (model_type.compare(std::string("CURVED_CYLINDERICAL_LITHOPHANE")) == 0)
	{
		const string image_filename = input_path + input_filename;

		IMAGE_2D image;
		image.ReadFileAndFitHeight(image_filename, 2 * target_height);
		image.ScaleTo(target_width, target_height);

		DEFINE_PARAMETER(int, mapping_mode);

		ParametricCurve curve;
		curve.initializeFromScript(par_reader, inner_radius, (T)1, height_scale, mapping_mode);

		const Vector3D<float> temp_v = curve.getPosition(1);
		const T length = curve.getLength();

		ShellOfRevolution shell_maker;
		shell_maker.generateLithophane(image, surface_, curve, inner_radius, base_thickness, front_thickness, num_smoothing, height_scale);
		shell_maker.generateLithophaneHeightMap(image, num_smoothing, height_map);
		shell_maker.generateLithophaneTexture(image, num_smoothing, texture);

		//image.ScaleTo(target_width/2, target_height/2);
		//shell_maker.GenerateLithophane(image, surface_low, curve, inner_radius, base_thickness, front_thickness, num_smoothing);

		image.ScaleTo(target_width / 4, target_height / 4);
		shell_maker.generateOuterSurface(image, outer_base_surface_, curve, inner_radius, base_thickness, 0.0f, num_smoothing);
		//shell_maker.GenerateInnerSurface(image, inner_base_surface, curve, inner_radius, 0.0f, num_smoothing);
	}
	else if (model_type.compare(std::string("FREEFORM_CURVED_CYLINDERICAL_LITHOPHANE")) == 0)
	{
		DEFINE_PARAMETER(string, input_shape_filename);

		const string shape_image_filename = input_path + input_shape_filename;
		const string litho_image_filename = input_path + input_filename;

		IMAGE_2D shape_image;
		shape_image.ReadFileAndFitHeight(shape_image_filename, 2 * target_height);
		shape_image.ScaleTo(target_width, target_height);

		IMAGE_2D litho_image;
		litho_image.ReadFileAndFitHeight(litho_image_filename, 2 * target_height);
		litho_image.ScaleTo(target_width, target_height);

		DEFINE_PARAMETER(int, mapping_mode);

		ParametricCurve curve;
		curve.initializeFromScript(par_reader, inner_radius, (T)1, height_scale, mapping_mode);

		ShellOfRevolution shell_maker;
		shell_maker.generateFreefromCurvedCylindericalLithophane(shape_image, litho_image, surface_, curve, base_thickness, front_thickness, base_thickness, inner_radius, num_smoothing, height_scale);

		BOX_3D<T> box = surface_.getAABB();

		std::cout << box.x_max_ - box.x_min_ << " x " << box.y_max_ - box.y_min_ << " x " << box.z_max_ - box.z_min_ << std::endl;
	}
	else if (model_type.compare(std::string("CIRCULAR_TRUNCATED_CONE")) == 0)
	{
//		CircularTruncatedCone

		const string image_filename = input_path + input_filename;

		IMAGE_2D image;
		image.ReadFileAndFitHeight(image_filename, 2 * target_height);
		image.ScaleTo(target_width, target_height);

		DEFINE_PARAMETER(float, upper_inner_radius);
		DEFINE_PARAMETER(float, bottom_inner_radius);
		DEFINE_PARAMETER(int, ceiling_thickness);

		LITHOPHANE lithophane_maker;
		image.ReflectLeftRight();
		image.Rotate90();
		lithophane_maker.InitializeCircularTruncatedCone(image, surface_, 0, upper_inner_radius, bottom_inner_radius, base_thickness, front_thickness, num_smoothing, ceiling_thickness);
	}

	surface_.findAdjacentTrianglesOfVertices();
	surface_.determineFaceAveragedVertexNormals();

	// Write STL file
// 	if (write_stl) surface_.WriteSTL((output_path + "/" + output_filename_prefix + ".stl").c_str());
// 	if (write_ctm) surface_.WriteCTM((output_path + "/" + output_filename_prefix + ".ctm").c_str());
// 
// 	if (write_outer_base_ctm) outer_base_surface_.WriteCTM((output_path + "/" + output_filename_prefix + "_outer.ctm").c_str());// includes texture uv
// 
// 	if (write_texture) texture.WriteBMP24((output_path + "/" + output_filename_prefix + "_texture.bmp").c_str());
// 	texture.ReverseColors();
// 	if (write_bumpmap) texture.WriteBMP24((output_path + "/" + output_filename_prefix + "_bump.bmp").c_str());
}


void LithophaneMaker::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		data_depot->reset();
	}
	END_ONE_THREAD_WORK;
	
	static PhongTrianglesData *phong_triangles_temp = nullptr;

	//	if (step_ >= 3)
	{
		BEGIN_ONE_THREAD_WORK
		{
			phong_triangles_temp = new PhongTrianglesData;

			phong_triangles_temp->rasterize_mode_ = 0x1B02;

			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);

			// 		phong_triangles_temp->positions_.Reset();
			// 		phong_triangles_temp->normals_.Reset();

			//	LevelsetUniform2D &levelset_(levelset_image_);

			surface_.use_face_normal_ = true;

			surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);
		}
		END_ONE_THREAD_WORK;

		// mesh view
// 		BEGIN_ONE_THREAD_WORK
// 		{
// 			const BOX_3D<T> aabb = surface_.getAABB();
// 
// 			surface_.translate(TV3(aabb.x_max_ - aabb.x_min_, (T)0, (T)0));
// 
// 			phong_triangles_temp = new PhongTrianglesData;
// 
// 			phong_triangles_temp->rasterize_mode_ = 0x1B01;
// 
// 			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);
// 
// 			// 		phong_triangles_temp->positions_.Reset();
// 			// 		phong_triangles_temp->normals_.Reset();
// 
// 			//	LevelsetUniform2D &levelset_(levelset_image_);
// 
// 			surface_.use_face_normal_ = true;
// 
// 			surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);
// 		}
// 		END_ONE_THREAD_WORK;
	}
}

BOX_3D<T> LithophaneMaker::getAABB()
{
	return surface_.getAABB();
}

void LithophaneMaker::writeFiles()
{
//	if (write_stl) surface.WriteSTL((output_path_ + "/" + output_filename_prefix_ + ".stl").c_str());
	surface_.writeSTL((output_path_ + "/" + output_filename_prefix_ + ".stl").c_str());
// 	if (write_ctm) surface.WriteCTM((output_path + "/" + output_filename_prefix + ".ctm").c_str());
// 
// 	if (write_outer_base_ctm) outer_base_surface.WriteCTM((output_path + "/" + output_filename_prefix + "_outer.ctm").c_str());// includes texture uv
// 
// 	if (write_texture) texture.WriteBMP24((output_path + "/" + output_filename_prefix + "_texture.bmp").c_str());
// 	texture.ReverseColors();
// 	if (write_bumpmap) texture.WriteBMP24((output_path + "/" + output_filename_prefix + "_bump.bmp").c_str());
}