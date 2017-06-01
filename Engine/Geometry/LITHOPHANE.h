// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "Image/IMAGE_2D.h"
#include "StaticTriangularSurface.h"
#include "DataStructure/GridUniform2D.h"
#include "DataStructure/Array2D.h"

class LITHOPHANE
{
public:
	void InitializeCylinder(const IMAGE_2D& image, StaticTriangularSurface& surface, const float inner_radius = 45.0f, const float base_thickness = 0.75f, const float front_thickness = 1.0f, const int smoothing = 5);
	void InitializePlane(const IMAGE_2D& image, StaticTriangularSurface& surface, const float base_depth = 5.0f, const float print_depth = 5.0f, const float print_width = 100.0f, const int smoothing = 5);
	void InitializePlane(const Array2D<T>& height_field, StaticTriangularSurface& surface, const float base_depth, const float print_depth, const float print_width);
	void InitializeExtrusion(const IMAGE_2D& image, StaticTriangularSurface& surface, const float thickness = 0.05f, const float width = 3.0f, const int smoothing = 5);
	void InitializePlanarKeyring(const IMAGE_2D& image, const IMAGE_2D& litho_image, StaticTriangularSurface& surface,
		const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float width = 3.0f, const int smoothing = 5);
	void InitializeFreeFormPlanarLithophane(const IMAGE_2D& image, StaticTriangularSurface& surface, const float base_thickness, const float front_print_thickness, const float width = 3.0f, const int smoothing = 5);
	void InitializeExtrusionAndSphericalMapping(const IMAGE_2D& image, StaticTriangularSurface& surface, const float thickness = 0.02f, const float width = 3.0f, const int smoothing = 5);
	void InitializePrism(const IMAGE_2D& image, StaticTriangularSurface& surface, const int num_sides = 4, const float inner_radius = 45.0f, const float base_thickness = 5.0f, const float print_thickness = 5.0f, const float corner_angle = 0.0f, const int smoothing = 5);
	void InitializeCircularTruncatedCone(const IMAGE_2D& image, StaticTriangularSurface& surface, 
		const int num_sides_ = 0,
		const float upper_radius = 45.0f, const float bottom_radius = 90.0f, const float base_thickness = 0.75f, const float front_thickness = 1.0f, const int smoothing = 5, const int thickness_ceiling = 0, const int thickness_floor = 0);

	void makeFreeformStamp(const Array2D<T>& shape_field_node, const Array2D<T>& thickness_field_node, StaticTriangularSurface& surface,
		const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float width = 3.0f, const int smoothing = 5);

	void makeFreeformStampDoubleSided(const Array2D<T>& shape_field_node, const Array2D<T>& thickness_field_node, StaticTriangularSurface& surface,
		const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float width = 3.0f, const int smoothing = 5);

	void makeFreePlane(const GridUniform2D& grid, const Array2D<T>& shape_field_node, const Array2D<T>& height_field_node, const T iso_threshold, const T min_edge_length, StaticTriangularSurface& surface);
};