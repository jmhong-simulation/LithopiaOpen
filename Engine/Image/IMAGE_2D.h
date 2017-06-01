// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "BMP_RGB.h"
#include "CONVENTIONAL_MACROS.h"
#include "GENERIC_DEFINITIONS.h"

#include "DataStructure/Array2D.h"
#include "Geometry/BOX_2D.h"
#include "Geometry/StaticTriangularSurface.h"

#include <string>

class IMAGE_2D
{
public:
	union{
		struct{ int res_x_; int res_y_; };
		struct{ int width_; int height_; };
	};

	Array2D<BMP_RGB> data_;

	IMAGE_2D()
	{}

	~IMAGE_2D()
	{}

public:
	void Initialize(const int res_x_input, const int res_y_input);

	void Initialize(const IMAGE_2D& image);
	void Initialize(const Array2D<T>& height_map);

	void ReadFileAndFitHeight(const std::string image_file_path, const int max_height);

	const bool ReadBMP24(const char * imagepath);
	const bool WriteBMP24(const char * imagepath);
	
	void Rotate90();
	void ReflectLeftRight();
	void ReflectUpDown();
	void ReverseColors();
	
	void ExtendOneColumn(const int column, const int width);
	void CutLeftAndRight(const int cut_width);
	void CutTopAndBottom(const int cut_width);

	void SetAllColor(const BMP_RGB& input_color);
	void SetMultipleAllColor(T val);

	void ScaleBy(const T scale_x, const T scale_y);
	void ScaleTo(const int new_width, const int new_height);
	void ScaleIso(const T scale);	// scale_x = scale_y
	void ScaleToFitHeight(const int new_height);
	void ScaleToFitWidth(const int new_width);

	const int GetResizeWidth(const int new_height) const;	// get expected width when this is image is forced to have new_height without ratio change.
	const int GetResizeHeight(const int new_width) const;

	const BOX_2D<int> GetColorBox(const BMP_RGB& color_mask) const;

	const T GetAspectRatio() const;
};
