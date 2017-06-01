// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "ParametricCurve.h"
#include "Image/IMAGE_2D.h"

class ShellOfRevolution
{
public:
	ShellOfRevolution()
	{}

	void generateLithophane(const IMAGE_2D& image, StaticTriangularSurface& surface,
		ParametricCurve& curve,		
		const float radius = 45.0f, const float base_thickness = 0.75f, const float front_thickness = 1.0f, const int smoothing_repeat = 5, const float height_scale = 1.0f);

	void generateOuterSurface(const IMAGE_2D& image, StaticTriangularSurface& surface,
		ParametricCurve& curve,
		const float radius = 45.0f, const float base_thickness = 0.75f, const float front_thickness = 1.0f, const int smoothing_repeat = 5);

	void generateInnerSurface(const IMAGE_2D& image, StaticTriangularSurface& surface,
		ParametricCurve& curve,
		const float radius = 45.0f, const float base_thickness = 0.75f, const int smoothing_repeat = 5);

	void generateLithophaneHeightMap(const IMAGE_2D& image, const int smoothing_repeat, IMAGE_2D& heightmap);
	void generateLithophaneTexture(const IMAGE_2D& image, const int smoothing_repeat, IMAGE_2D& heightmap);

	void generateFreefromCurvedCylindericalLithophane(const IMAGE_2D& shape_image, const IMAGE_2D& height_image, 
		StaticTriangularSurface& surface,
		ParametricCurve& curve,
		const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float radius, const int smoothing_repeat, const float height_scale_input);

	const TV getShellX(const TV2& uv, ParametricCurve& curve, const T radial_thickness) const;

	// move these functions
	const TV ShellOfRevolution::cast(const TV2 v, const T z) const;
	const bool ShellOfRevolution::isInterfacial(const T& phi0, const T& phi1, const T& th) const;
	const bool ShellOfRevolution::isInterfacial(const Array2D<T>& field, const int i, const int j, const T th) const;
};