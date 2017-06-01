// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/LinkedArray.h"
#include "DataStructure/GridUniform2D.h"
#include "Geometry/StaticTriangularSurface.h"

#include "LINE_SEGMENT.h"
#include "StaticTriangle.h"
#include "PLANE.h"

class SLICER
{
public:
	Array1D<LinkedArray<LINE_SEGMENT>> z_sorted_contours_;

	GridUniform2D grid_;
	Array1D<Array2D<T>> signed_distance_layers_;

	T z_min_, dz_;
	T num_slices_;

	SLICER() : num_slices_(20)
	{}

	const TV GetSignedDistanceWeightedAverage(const T& phi0, const T& phi1, const TV& v0, const TV& v1);
	void SliceTriangle(const StaticTriangle triangle, const PLANE plane, LinkedArray<LINE_SEGMENT>& contour);
	void SortContours(LinkedArray<LINE_SEGMENT>& contour);
	const int FindConnectedContour(const int i, Array1D<LINE_SEGMENT>& contour_temp, const Array1D<bool>& contour_check, LinkedArray<LINE_SEGMENT>& contour, const double ep_sqr);
	void SliceOutlines(const StaticTriangularSurface& surface);
	void SliceDLP(const StaticTriangularSurface& surface);
	void CopyAllContoursToArray(Array1D<LINE_SEGMENT>& all_contours) const;
	void GenerateRaft(std::ofstream& file, T& epsilon) const;
	void GenerateGCode() const;
};