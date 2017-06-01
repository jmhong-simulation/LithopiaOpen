// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.

#include "../Engine/Framework/WORLD.h"
#include "../Engine/Geometry/LITHOPHANE.h"
#include "../Engine/Geometry/SHELL_OF_REVOLUTION.h"
#include "../Engine/DataStructure/GridUniform2D.h"
#include "../Engine/CONVENTIONAL_MACROS.h"
#include <boost/date_time/gregorian/gregorian.hpp>

#include <IL/il.h>
#include <IL/ilu.h>

#define CYLINDERICAL_LITHOPHANE 0
#define SPHERICAL_PROJECTION 1
#define PLANAR_LITHOPHANE 2
#define TRIANGULAR_PRISM 3
#define RECTANGULAR_PRISM 4
#define PENTAGONAL_PRISM 5
#define HEXAGONARL_PRISM 6
#define PLANE_WITH_BLACK_HOLES 7
#define PLANE_WITH_WHITE_HOLES 8
#define PLANE_KEY_RING 9
#define PLANE_LITHOPHANE_FREE 10
#define CIRCULAR_TRUNCATED_CONE_LITHOPHANE 11
#define CURVED_CYLINDERICAL_LITHOPHANE 12
#define FREEFORM_CURVED_CYLINDERICAL_LITHOPHANE 13

//#define KEYRING_RELEASE 1

using namespace std;

string image_file_name;

void ilResizeImage(const int& ground_width, const int& ground_height);
void InitializeImageFromFile(IMAGE_2D& image, const string image_file_path, const int& ground_width, const int& ground_height);
void ilReadImageFromFile(const string image_file_path);
void CopyImageFromILImage(IMAGE_2D& image, const string image_file_path);
void PrepareForWriting(StaticTriangularSurface& surface);

int main(int argc, char *argv[])
{
	ilInit();

	const int object_type = std::atoi(argv[1]);

	string output_directory(argv[2]);

	StaticTriangularSurface surface;

#ifdef KEYRING_RELEASE
	std::cout << "This keyring maker is developed by Jeong-Mo Hong." << std::endl;

	std::string user_name = "marantz2000";
	boost::gregorian::date expiration_day(boost::gregorian::from_simple_string(std::string("2015-09-01")));

	if (object_type != PLANE_KEY_RING)
	{
		std::cout << "The first parameter should be 9. Terminating." << std::endl;
		return -1;
	}
	else
	{
		using namespace boost::gregorian;

		date today = day_clock::local_day();
//		date today(from_simple_string(std::string("2015-09-03")));

		std::cout << "This release is allowed for " << user_name << " until " << to_simple_string(expiration_day) << std::endl;

		days days_left = expiration_day - today;

		if (days_left.is_negative())
		{
			std::cout << "License has been expired. Terminating." << std::endl;

			return -1;
		}
		else
		{
			std::cout << "License is valid. OK to proceed." << std::endl;
		}

// 		date::ymd_type ymd = today.year_month_day();
// 		greg_weekday wd = today.day_of_week();
// 		std::cout << wd.as_long_string() << " " << ymd.month.as_long_string() << " " << ymd.day << ", " << ymd.year << std::endl;
	}
#endif

#ifndef KEYRING_RELEASE
	if (object_type == CYLINDERICAL_LITHOPHANE)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);

		const float inner_radius = std::atof(argv[6]);
		const float base_thickness = std::atof(argv[7]);
		const float front_thickness = std::atof(argv[8]);
		const int smooth = std::atoi(argv[9]);

		LITHOPHANE lithophane_maker;

		image.ReflectLeftRight();
		image.Rotate90();
		lithophane_maker.InitializeCylinder(image, surface, inner_radius, base_thickness, front_thickness, smooth);
	}
	else if (object_type == SPHERICAL_PROJECTION)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);

		const float print_thickness = std::atof(argv[6]);
		const float print_width = std::atof(argv[7]);
		const int smooth = std::atoi(argv[8]);

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializeExtrusionAndSphericalMapping(image, surface, print_thickness, print_width, smooth);
	}
	else if (object_type == PLANAR_LITHOPHANE)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);

		const float base_depth = std::atof(argv[6]);
		const float print_depth = std::atof(argv[7]);
		const float print_width = std::atof(argv[8]);
		const int smooth = std::atoi(argv[9]);

		LITHOPHANE lithophane_maker;
		image.ReflectLeftRight();
		lithophane_maker.InitializePlane(image, surface, base_depth, print_depth, print_width, smooth);
	}
	else if (object_type >= TRIANGULAR_PRISM && object_type <= HEXAGONARL_PRISM)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);

		const int num_sides = object_type;				// TRIANGULAR_PRISM = 3, ...
		const float inner_radius = std::atof(argv[6]);
		const float base_thickness = std::atof(argv[7]);
		const float print_thickness = std::atof(argv[8]);
		const float corner_angle = std::atof(argv[9]);		// degree
		const int smoothing = std::atoi(argv[10]);

		image.ReflectLeftRight();
		image.Rotate90();

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializePrism(image, surface, num_sides, inner_radius, base_thickness, print_thickness, corner_angle, smoothing);
	}
	else if (object_type == PLANE_WITH_BLACK_HOLES)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);
		image.ReflectLeftRight();
		image.ReflectUpDown();

		const float thickness = std::atof(argv[6]);
		const float width = std::atof(argv[7]);
		const int smoothing = std::atoi(argv[8]);

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializeExtrusion(image, surface, thickness, width, smoothing);
	}
	else if (object_type == PLANE_WITH_WHITE_HOLES)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);
		image.ReverseColors();
		image.ReflectLeftRight();
		image.ReflectUpDown();

		const float thickness = std::atof(argv[6]);
		const float width = std::atof(argv[7]);
		const int smoothing = std::atoi(argv[8]);

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializeExtrusion(image, surface, thickness, width, smoothing);
	}
	else if (object_type == PLANE_KEY_RING)
	{
		string shape_image_file_path(argv[3]);
		string word_image_file_path(argv[4]);
		string litho_image_file_path(argv[5]);

		const int ground_width = std::atoi(argv[6]);
		const int ground_height = std::atoi(argv[7]);

		const float base_thickness = std::atof(argv[8]);
		const float front_print_thickness = std::atof(argv[9]);
		const float back_print_thickness = std::atof(argv[10]);
		const float width = std::atof(argv[11]);
		const int smoothing = std::atoi(argv[12]);

		// Read Image
		IMAGE_2D shape_image, word_image, litho_image;

		InitializeImageFromFile(shape_image, shape_image_file_path, ground_width, ground_height);
		InitializeImageFromFile(word_image, word_image_file_path, ground_width, ground_height);
		InitializeImageFromFile(litho_image, litho_image_file_path, ground_width, ground_height);

		const BOX_2D<int> blue_box = shape_image.GetColorBox(BMP_RGB(0, 0, 255));
		const T scale = MIN2((T)(blue_box.i_end_ - blue_box.i_start_ + 1) / (T)litho_image.res_x_, (T)(blue_box.j_end_ - blue_box.j_start_ + 1) / (T)litho_image.res_y_);
		litho_image.ScaleBy(scale, scale);
		litho_image.ReverseColors();
		litho_image.SetMultipleAllColor(0.5f);
		litho_image.ReverseColors();

		// litho_image + word_image
		const int i_offset = MAX(0, ((blue_box.i_end_ - blue_box.i_start_ + 1) - litho_image.res_x_)) / 2;
		const int j_offset = MAX(0, ((blue_box.j_end_ - blue_box.j_start_ + 1) - litho_image.res_y_)) / 2;

		for (int j = litho_image.data_.j_start_; j <= litho_image.data_.j_end_; ++j)
		for (int i = litho_image.data_.i_start_; i <= litho_image.data_.i_end_; ++i)
		{
			int ii = i + blue_box.i_start_ + i_offset;
			int jj = j + blue_box.j_start_ + j_offset;

			word_image.data_(ii, jj) = litho_image.data_.getClamped(i, j);
		}

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializePlanarKeyring(shape_image, word_image, surface, base_thickness, front_print_thickness, back_print_thickness, width, smoothing);
	}
	else if (object_type == PLANE_LITHOPHANE_FREE)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);
		//		image.ReverseColors();	// be careful

		const float base_thickness = std::atof(argv[6]);
		const float print_thickness = std::atof(argv[7]);
		const float width = std::atof(argv[8]);
		const int smoothing = std::atoi(argv[9]);

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializeFreeFormPlanarLithophane(image, surface, base_thickness, print_thickness, width, smoothing);
	}
	else if (object_type == CIRCULAR_TRUNCATED_CONE_LITHOPHANE)
	{
		string image_file_path(argv[3]);
		const int ground_width = std::atoi(argv[4]);
		const int ground_height = std::atoi(argv[5]);

		IMAGE_2D image;
		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);

		const float upper_inner_radius = std::atof(argv[6]);
		const float bottom_inner_radius = std::atof(argv[7]);
		const float base_thickness = std::atof(argv[8]);
		const float front_thickness = std::atof(argv[9]);
		const int smooth = std::atoi(argv[10]);

		LITHOPHANE lithophane_maker;

		image.ReflectLeftRight();
		image.Rotate90();
		lithophane_maker.InitializeCircularTruncatedCone(image, surface, 0, upper_inner_radius, bottom_inner_radius, base_thickness, front_thickness, smooth);
	}
	else if (object_type == CURVED_CYLINDERICAL_LITHOPHANE)
	{
		// disabled.
// 		string image_file_path(argv[3]);
// 		const int ground_width = std::atoi(argv[4]);
// 		const int ground_height = std::atoi(argv[5]);
// 
// 		IMAGE_2D image;
// 		InitializeImageFromFile(image, image_file_path, ground_width, ground_height);
// 
// 		const float inner_radius = std::atof(argv[6]);
// 		const float base_thickness = std::atof(argv[7]);
// 		const float front_thickness = std::atof(argv[8]);
// 		const int smooth = std::atoi(argv[9]);
// 
// 		//TODO: use input
// 		ParametricCurve curve;
// 		curve.
// 		curve.x_.InitializeZeroOneDomain(1, 0, (T)(34 * 2) / (T)(46 * 2), -0.8); // radius direction
// 		curve.y_.InitializeZeroOneDomain(0, 0, 0, 0);
// 		curve.z_.InitializeZeroOneDomain(0, 1, 1, 1);  // height
// 
// 		SHELL_OF_REVOLUTION shell_maker;
// 		shell_maker.GenerateLithophane(image, surface, curve, inner_radius, base_thickness, front_thickness, smooth);
	}
	else if (object_type == FREEFORM_CURVED_CYLINDERICAL_LITHOPHANE)
	{
/*
		string shape_image_file_path(argv[3]);
		string litho_image_file_path(argv[4]);
		const int ground_width = std::atoi(argv[5]);
		const int ground_height = std::atoi(argv[6]);

		ilReadImageFromFile(shape_image_file_path);

		IMAGE_2D image;
		CopyImageFromILImage(image, shape_image_file_path);

		const BOX_2D<int> blue_box = image.GetColorBox(BMP_RGB(0, 0, 255));

//		image.ExtendOneColumn(image.res_x_ / 2, 20);

		IMAGE_2D input_image;
		ilReadImageFromFile(litho_image_file_path);
		CopyImageFromILImage(input_image, litho_image_file_path);

		const T scale = MIN2((T)(blue_box.i_end_ - blue_box.i_start_ + 1) / (T)input_image.res_x_, (T)(blue_box.j_end_ - blue_box.j_start_ + 1) / (T)input_image.res_y_);

		input_image.ScaleBy(scale, scale);

		IMAGE_2D litho_image;
		litho_image.Initialize(image.res_x_, image.res_y_);
		litho_image.SetAllColor(BMP_RGB(0, 0, 0));

		// copy image to litho height map
		
		// centering
		const int i_offset = MAX(0, ((blue_box.i_end_ - blue_box.i_start_ + 1) - input_image.res_x_)) / 2;
		const int j_offset = MAX(0, ((blue_box.j_end_ - blue_box.j_start_ + 1) - input_image.res_y_)) / 2;

		for (int j = input_image.data_.j_start_; j <= input_image.data_.j_end_; ++j)
		for (int i = input_image.data_.i_start_; i <= input_image.data_.i_end_; ++i)
		{
			litho_image.data_(i + blue_box.i_start_ + i_offset, j + blue_box.j_start_ + j_offset) = input_image.data_.getClamped(i, j);
		}

//		litho_image.ExtendOneColumn(litho_image.res_x_ / 2, 20);

		//TODO: auto extending
		//TODO: scaling

		const float base_thickness = std::atof(argv[7]);
		const float front_print_thickness = std::atof(argv[8]);
		const float back_print_thickness = std::atof(argv[9]);
		const float radius = std::atof(argv[10]);
		const int smoothing = std::atoi(argv[11]);

		//TODO: use input
		PARAMETRIC_CURVE_SEGMENT curve;
		curve.x_.InitializeZeroOneDomain(1, 0, 1, 0);	// radius direction
		curve.y_.InitializeZeroOneDomain(0, 0, 0, 0);
		curve.z_.InitializeZeroOneDomain(0, 1, 1, 1);		// height

		SHELL_OF_REVOLUTION shell_maker;
		shell_maker.GenerateFreefromCurvedCylindericalLithophane(image, litho_image, surface, curve, base_thickness, front_print_thickness, back_print_thickness, radius, smoothing);
*/
	}
	else
	{
		std::cout << "Object IX is invalid" << std::endl;
		exit(1);
	}

#else if
	if (object_type == PLANE_KEY_RING)
	{
		string shape_image_file_path(argv[3]);
		string word_image_file_path(argv[4]);
		string litho_image_file_path(argv[5]);

		const int ground_width = std::atoi(argv[6]);
		const int ground_height = std::atoi(argv[7]);

		const float base_thickness = std::atof(argv[8]);
		const float front_print_thickness = std::atof(argv[9]);
		const float back_print_thickness = std::atof(argv[10]);
		const float width = std::atof(argv[11]);
		const int smoothing = std::atoi(argv[12]);

		// Read Image
		IMAGE_2D shape_image, word_image, litho_image;

		InitializeImageFromFile(shape_image, shape_image_file_path, ground_width, ground_height);
		InitializeImageFromFile(word_image, word_image_file_path, ground_width, ground_height);
		InitializeImageFromFile(litho_image, litho_image_file_path, ground_width, ground_height);

		const BOX_2D<int> blue_box = shape_image.GetColorBox(BMP_RGB(0, 0, 255));
		const T scale = MIN2((T)(blue_box.i_end_ - blue_box.i_start_ + 1) / (T)litho_image.res_x_, (T)(blue_box.j_end_ - blue_box.j_start_ + 1) / (T)litho_image.res_y_);
		litho_image.ScaleBy(scale, scale);
		litho_image.ReverseColors();
		litho_image.SetMultipleAllColor(0.5f);
		litho_image.ReverseColors();

		// litho_image + word_image
		const int i_offset = MAX(0, ((blue_box.i_end_ - blue_box.i_start_ + 1) - litho_image.res_x_)) / 2;
		const int j_offset = MAX(0, ((blue_box.j_end_ - blue_box.j_start_ + 1) - litho_image.res_y_)) / 2;

		for (int j = litho_image.data_.j_start_; j <= litho_image.data_.j_end_; ++j)
		for (int i = litho_image.data_.i_start_; i <= litho_image.data_.i_end_; ++i)
		{
			int ii = i + blue_box.i_start_ + i_offset;
			int jj = j + blue_box.j_start_ + j_offset;

			word_image.data_(ii, jj) = litho_image.data_.GetClamped(i, j);
		}

		LITHOPHANE lithophane_maker;
		lithophane_maker.InitializePlanarKeyring(shape_image, word_image, surface, base_thickness, front_print_thickness, back_print_thickness, width, smoothing);
	}
#endif

	PrepareForWriting(surface);

	// Write STL file
	const string output_file_path_and_name_stl = output_directory + "/" + image_file_name + ".stl";
	surface.writeSTL(output_file_path_and_name_stl.c_str());

	// Write CTM file
	const string output_file_path_and_name_ctm = output_directory + "/" + image_file_name + ".ctm";
	surface.writeCTM(output_file_path_and_name_ctm.c_str());
	
	return 0;
}

void ilResizeImage(const int& ground_width, const int& ground_height)
{
	// image resizing
	int depth = ilGetInteger(IL_IMAGE_DEPTH);

	const int image_width = ilGetInteger(IL_IMAGE_WIDTH);
	const int image_height = ilGetInteger(IL_IMAGE_HEIGHT);

	const float dw = (float)ground_width / (float)image_width;
	const float dh = (float)ground_height / (float)image_height;

	const float dd = MIN2(dw, dh);

	const int final_width = (int)((float)image_width*dd);
	const int final_height = (int)((float)image_height*dd);

	iluScale(final_width, final_height, depth);
}

void ilReadImageFromFile(const string image_file_path)
{
	// Find file name from path
	if ((int)image_file_path.find("\\") > -1 || (int)image_file_path.find("/") > -1)
	{
		int ix = MAX2((int)image_file_path.find_last_of("\\"), (int)image_file_path.find_last_of("/"));
		image_file_name = image_file_path.substr(ix + 1);
		image_file_name;
	}
	else
	{
		image_file_name = image_file_path;
	}

	if (ilLoadImage(image_file_path.c_str()) == false)
	{
		std::cout << "Cannot read image file " << image_file_path.c_str() << std::endl;
		exit(1);
	}
}

void CopyImageFromILImage(IMAGE_2D& image, const string image_file_path)
{
	int width = ilGetInteger(IL_IMAGE_WIDTH);
	int height = ilGetInteger(IL_IMAGE_HEIGHT);
	int size = ilGetInteger(IL_IMAGE_SIZE_OF_DATA);
	int ch = ilGetInteger(IL_IMAGE_CHANNELS);

	if (ch < 3)
	{
		std::cout << "Image channel < 3" << std::endl;
		exit(1);
	}

	ILubyte * bytes = ilGetData();

	image.res_x_ = width;
	image.res_y_ = height;
	image.data_.initialize(0, 0, width, height, false);

	if (IL_BMP == ilTypeFromExt(image_file_path.c_str()))
	{
		// for bmp
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int ix = (j*width + i)*ch;

				image.data_(i, j).r_ = bytes[ix + 2];
				image.data_(i, j).g_ = bytes[ix + 1];
				image.data_(i, j).b_ = bytes[ix + 0];
			}
		}
	}
	else
	{
		// for others
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int ix = ((height - 1 - j)*width + i)*ch;

				image.data_(i, j).r_ = bytes[ix + 0];
				image.data_(i, j).g_ = bytes[ix + 1];
				image.data_(i, j).b_ = bytes[ix + 2];
			}
		}
	}
}

void InitializeImageFromFile(IMAGE_2D& image, const string image_file_path, const int& ground_width, const int& ground_height)
{
	// Find file name from path
	if ((int)image_file_path.find("\\") > -1 || (int)image_file_path.find("/") > -1)
	{
		int ix = MAX2((int)image_file_path.find_last_of("\\"), (int)image_file_path.find_last_of("/"));
		image_file_name = image_file_path.substr(ix + 1);
	}
	else
	{
		image_file_name = image_file_path;
	}

	if (ilLoadImage(image_file_path.c_str()) == false)
	{
		std::cout << "Cannot read image file " << image_file_path.c_str() << std::endl;
		exit(1);
	}

	ilResizeImage(ground_width, ground_height);

	int width = ilGetInteger(IL_IMAGE_WIDTH);
	int height = ilGetInteger(IL_IMAGE_HEIGHT);
	int size = ilGetInteger(IL_IMAGE_SIZE_OF_DATA);
	int ch = ilGetInteger(IL_IMAGE_CHANNELS);

	if (ch < 3)
	{
		std::cout << "Image channel < 3" << std::endl;
		exit(1);
	}

	ILubyte * bytes = ilGetData();

	image.res_x_ = width;
	image.res_y_ = height;
	image.data_.initialize(0, 0, width, height, false);

	if (IL_BMP == ilTypeFromExt(image_file_path.c_str()))
	{
		// for bmp
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int ix = (j*width + i)*ch;

				image.data_(i, j).r_ = bytes[ix + 2];
				image.data_(i, j).g_ = bytes[ix + 1];
				image.data_(i, j).b_ = bytes[ix + 0];
			}
		}
	}
	else
	{
		// for others
		for (int j = 0; j < height; j++)
		{
			for (int i = 0; i < width; i++)
			{
				int ix = ((height - 1 - j)*width + i)*ch;

				image.data_(i, j).r_ = bytes[ix + 0];
				image.data_(i, j).g_ = bytes[ix + 1];
				image.data_(i, j).b_ = bytes[ix + 2];
			}
		}
	}
}

void PrepareForWriting(StaticTriangularSurface& surface)
{
	BOX_3D<T> aabb = surface.getAABB();

	float scale = MAX3(aabb.x_max_ - aabb.x_min_, aabb.y_max_ - aabb.y_min_, aabb.z_max_ - aabb.z_min_);
	if (scale != 0.0f) scale = 1.0f / scale;

	glm::vec3 scalevec(scale, scale, scale);
	const Vector3D<float> aabb_center = aabb.GetCenter();
	const glm::vec3 center(aabb_center.x_, aabb_center.y_, aabb_center.z_);

	glm::mat4 m_matrix;
	m_matrix = glm::scale(scalevec) * glm::translate(-center);

	surface.OOBB_.Initialize(aabb);

//	std::cout << surface.vertex_positions_.num_elements_ << std::endl;

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();

	surface.translate(-aabb_center);
}