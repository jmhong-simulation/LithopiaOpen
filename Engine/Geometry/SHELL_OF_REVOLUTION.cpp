#include "SHELL_OF_REVOLUTION.h"
#include "UV_TO_XYZ.h"
#include "../DataStructure/Vector2D.h"
#include "../DataStructure/GridUniform2D.h"
#include "../Geometry/LINE_SEGMENT.h"
#include "../Operations/SMOOTHING_UNIFORM_2D.h"

#include <iostream>

const TV ShellOfRevolution::getShellX(const TV2& uv, ParametricCurve& curve, const T radial_thickness) const
{
	const T theta = 2.0f*PI*uv.u_;
	const TV p_on_curve = curve.getPosition(uv.v_);

	const T radius = p_on_curve.x_ + radial_thickness;
	const T height = p_on_curve.z_;

	return TV(radius*cos(theta), radius*sin(theta), height);
}

void ShellOfRevolution::generateLithophaneTexture(const IMAGE_2D& image, const int smoothing_repeat, IMAGE_2D& texture)
{
	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			// average 4 image cell colors to determine 1 node height
			height_field_node(i, j) = (image.data_.getIRepeatedJClamped(i - 1, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray()
				+ image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j).GetReversedGray()) * 0.25f;
		}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIRepeatJClamp(grid, height_field_node, smoothing_repeat);

	texture.Initialize(height_field_node);
	texture.ReverseColors();
}

void ShellOfRevolution::generateLithophaneHeightMap(const IMAGE_2D& image, const int smoothing_repeat, IMAGE_2D& height_map)
{
	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			// average 4 image cell colors to determine 1 node height
			height_field_node(i, j) = (image.data_.getIRepeatedJClamped(i - 1, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray()
				+ image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j).GetReversedGray()) * 0.25f;
		}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIRepeatJClamp(grid, height_field_node, smoothing_repeat);

	height_map.Initialize(height_field_node);
}

void ShellOfRevolution::generateLithophane(const IMAGE_2D& image, StaticTriangularSurface& surface, ParametricCurve& curve, const float inner_radius, const float base_thickness, const float front_print_thickness, const int smoothing_repeat, const float height_scale_input)
{
	const T outer_radius = inner_radius + base_thickness + front_print_thickness;
	const T circumference = 2.0f*PI*outer_radius;
//	const T height_scale = circumference / (T)(image.res_x_ + 1) * (T)image.res_y_ / curve.getLength() * height_scale_input;	// one more x cell
	const T height_scale = circumference / (T)(image.res_x_ + 1) * (T)image.res_y_ * height_scale_input;

	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		// average 4 image cell colors to determine 1 node height
		height_field_node(i, j) = (image.data_.getIRepeatedJClamped(i - 1, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray()
			+ image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j).GetReversedGray()) * 0.25f;
	}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIRepeatJClamp(grid, height_field_node, smoothing_repeat);

	surface.vertex_positions_.initialize(grid.GetNumAllNodes() * 2);
	//surface.vertex_uv_.initialize(grid.GetNumAllNodes() * 2);

	// define vertices
	int ix = 0;

	// vertices of outer cylinder
	for (int j = 0; j < height_field_node.j_res_; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		const T outer_radius = inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;
		const TV2 uv = TV2((T)i / (T)height_field_node.i_res_, (T)j / (height_field_node.j_res_ - 1));			// v = 1 when j = j_res_ - 1
		//surface.vertex_uv_[ix] = uv;
		surface.vertex_positions_[ix++] = getShellX(uv, curve, base_thickness + height_field_node(i, j)*front_print_thickness);
	}

	// all vertices of inner cylinder plane
	for (int j = 0; j < height_field_node.j_res_; ++j)
	for (int i = 0; i < height_field_node.i_res_; ++i)
	{
		const TV2 uv = TV2((T)i / (T)height_field_node.i_res_, (T)j / (height_field_node.j_res_ - 1));			// v = 1 when j = j_res_ - 1
		//surface.vertex_uv_[ix] = uv;
		surface.vertex_positions_[ix++] = getShellX(uv, curve, (T)0);
	}

	LinkedArray<TV_INT> tri;

	// triangles of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
	for (int i = 0; i < height_field_node.i_res_; i++)
	{
		const int i_plus1 = (i == height_field_node.i_res_ - 1) ? 0 : i + 1;		// repeated i + 1

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i_plus1, j), height_field_node.get1DIndex(i, j + 1));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i_plus1, j), height_field_node.get1DIndex(i_plus1, j + 1));
	}

	const int offset = grid.GetNumAllNodes();		// one j row is shared

	// triangles of inner cylinder
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
	for (int i = 0; i < height_field_node.i_res_; i++)
	{
		const int i_plus1 = (i == height_field_node.i_res_ - 1) ? 0 : i + 1;		// repeated i + 1

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i_plus1, j) + offset, height_field_node.get1DIndex(i, j + 1) + offset).getReversedCW();
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1) + offset, height_field_node.get1DIndex(i_plus1, j) + offset, height_field_node.get1DIndex(i_plus1, j + 1) + offset).getReversedCW();
	}
	
	//	bottom holed disk
	int j = 0;
	for (int i = 0; i < height_field_node.i_res_; i++)
	{
		const int i_plus1 = (i == height_field_node.i_res_ - 1) ? 0 : i + 1;		// repeated i + 1

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i_plus1, j) + offset, height_field_node.get1DIndex(i, j));
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i_plus1, j) + offset, height_field_node.get1DIndex(i_plus1, j));
	}

	// up holed disk

	j = height_field_node.j_res_ - 1;
	for (int i = 0; i < height_field_node.i_res_; i++)
	{
		const int i_plus1 = (i == height_field_node.i_res_ - 1) ? 0 : i + 1;		// repeated i + 1

		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j) + offset, height_field_node.get1DIndex(i_plus1, j) + offset, height_field_node.get1DIndex(i, j)).getReversedCW();
		tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i_plus1, j) + offset, height_field_node.get1DIndex(i_plus1, j)).getReversedCW();
	}

	tri.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();

	// aspect ratio sampling
	T aspect_ratio_sum = (T)0;	
	for (int j = 0; j < height_field_node.j_res_-1; ++j)
	{
		int i = 0;
		const TV2 uv_00 = TV2((T)i / (T)height_field_node.i_res_, (T)j / (height_field_node.j_res_ - 1));
		const TV2 uv_10 = TV2((T)(i + 1) / (T)height_field_node.i_res_, (T)j / (height_field_node.j_res_ - 1));
		const TV2 uv_01 = TV2((T)i / (T)height_field_node.i_res_, (T)(j+1) / (height_field_node.j_res_ - 1));

		const T dx = (getShellX(uv_00, curve, (T)0) - getShellX(uv_10, curve, (T)0)).getMagnitude();
		const T dy = (getShellX(uv_00, curve, (T)0) - getShellX(uv_01, curve, (T)0)).getMagnitude();

		aspect_ratio_sum += dx / dy;
	}

	std::cout << "Average aspecr ration (dx / dy) " << aspect_ratio_sum / (T)(height_field_node.j_res_-1);
}

void ShellOfRevolution::generateOuterSurface(const IMAGE_2D& image, StaticTriangularSurface& surface, ParametricCurve& curve, const float inner_radius, const float base_thickness, const float front_print_thickness, const int smoothing_repeat)
{
	const T outer_radius = inner_radius + base_thickness + front_print_thickness;
	const T circumference = 2.0f*PI*outer_radius;
	const T height_scale = circumference / (T)(image.res_x_ + 1) * (T)image.res_y_ / curve.getLength();	// one more x cell

	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			// average 4 image cell colors to determine 1 node height
			height_field_node(i, j) = (image.data_.getIRepeatedJClamped(i - 1, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray()
				+ image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j).GetReversedGray()) * 0.25f;
		}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIRepeatJClamp(grid, height_field_node, smoothing_repeat);		//TODO: implement

	surface.vertex_positions_.initialize(grid.GetNumAllNodes() * 2);
	surface.vertex_uv_.initialize(surface.vertex_positions_.num_elements_);

	// define vertices
	int ix = 0;

	// vertices of outer cylinder
	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			const T outer_radius = inner_radius + base_thickness + height_field_node(i, j)*front_print_thickness;
			const TV2 uv = TV2((T)i / (T)height_field_node.i_res_, (T)j / (height_field_node.j_res_ - 1));			// v = 1 when j = j_res_ - 1

			surface.vertex_uv_[ix] = uv;
			surface.vertex_positions_[ix++] = getShellX(uv, curve, base_thickness + height_field_node(i, j)*front_print_thickness);
		}

	LinkedArray<TV_INT> tri;

	// triangles of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
		for (int i = 0; i < height_field_node.i_res_; i++)
		{
			const int i_plus1 = (i == height_field_node.i_res_ - 1) ? 0 : i + 1;		// repeated i + 1

			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i_plus1, j), height_field_node.get1DIndex(i, j + 1));
			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i_plus1, j), height_field_node.get1DIndex(i_plus1, j + 1));
		}

	tri.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();
}

void ShellOfRevolution::generateInnerSurface(const IMAGE_2D& image, StaticTriangularSurface& surface, ParametricCurve& curve, const float inner_radius, const float front_print_thickness, const int smoothing_repeat)
{
	const T circumference = 2.0f*PI*inner_radius;
	const T height_scale = circumference / (T)(image.res_x_ + 1) * (T)image.res_y_ / curve.getLength();	// one more x cell

	GridUniform2D grid(0, 0, image.res_x_, image.res_y_, 0, 0, 1, 1);

	Array2D<T> height_field_node;
	grid.InitializeNodeArray(height_field_node);

	height_field_node.assignAllValues(0.0f);

	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			// average 4 image cell colors to determine 1 node height
			height_field_node(i, j) = (image.data_.getIRepeatedJClamped(i - 1, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray()
				+ image.data_.getIRepeatedJClamped(i, j - 1).GetReversedGray() + image.data_.getIRepeatedJClamped(i, j).GetReversedGray()) * 0.25f;
		}

	SMOOTHING_UNIFORM_2D smoothing;
	smoothing.SmoothIRepeatJClamp(grid, height_field_node, smoothing_repeat);		//TODO: implement

	surface.vertex_positions_.initialize(grid.GetNumAllNodes() * 2);

	// define vertices
	int ix = 0;

	// vertices of outer cylinder
	for (int j = 0; j < height_field_node.j_res_; ++j)
		for (int i = 0; i < height_field_node.i_res_; ++i)
		{
			const T outer_radius = inner_radius;
			const TV2 uv = TV2((T)i / (T)height_field_node.i_res_, (T)j / (height_field_node.j_res_ - 1));			// v = 1 when j = j_res_ - 1

			surface.vertex_positions_[ix++] = getShellX(uv, curve, 0);
		}

	LinkedArray<TV_INT> tri;

	// triangles of outer cylinder
	for (int j = 0; j < height_field_node.j_res_ - 1; j++)
		for (int i = 0; i < height_field_node.i_res_; i++)
		{
			const int i_plus1 = (i == height_field_node.i_res_ - 1) ? 0 : i + 1;		// repeated i + 1

			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j), height_field_node.get1DIndex(i_plus1, j), height_field_node.get1DIndex(i, j + 1)).getReversedCW();
			tri.PushBack() = TV_INT(height_field_node.get1DIndex(i, j + 1), height_field_node.get1DIndex(i_plus1, j), height_field_node.get1DIndex(i_plus1, j + 1)).getReversedCW();
		}

	tri.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.determineFaceAveragedVertexNormals();
}

const TV ShellOfRevolution::cast(const TV2 v, const T z) const
{
	return TV(v.x_, v.y_, z);
}

const bool ShellOfRevolution::isInterfacial(const T& phi0, const T& phi1, const T& th) const
{
	if (phi0 > th && phi1 > th) return false;
	else if (!(phi0 > th) && !(phi1 > th)) return false;

	return true;
}

const bool ShellOfRevolution::isInterfacial(const Array2D<T>& field, const int i, const int j, const T th) const
{
	if (field(i, j) > th)
	{
		if (field.getClamped(i + 1, j) <= th || field.getClamped(i - 1, j) <= th || field.getClamped(i, j + 1) <= th || field.getClamped(i, j - 1) <= th)
		{
			return true;
		}
	}
	else if (field(i, j) <= th)
	{
		if (field.getClamped(i + 1, j) > th || field.getClamped(i - 1, j) > th || field.getClamped(i, j + 1) > th || field.getClamped(i, j - 1) > th)
		{
			return true;
		}
	}

	return false;
}

void ShellOfRevolution::generateFreefromCurvedCylindericalLithophane(const IMAGE_2D& shape_image, const IMAGE_2D& height_image, StaticTriangularSurface& surface, ParametricCurve& curve, const float base_thickness, const float front_print_thickness, const float back_print_thickness, const float radius, const int smoothing_repeat, const float height_scale_input)
{
//	const T circumference = 2.0f*PI*radius;
//	const T height_scale = circumference / (T)(shape_image.res_x_ + 1) * (T)shape_image.res_y_ / curve.getLength() * height_scale_input;	// one more x cell
//	const T height_scale = circumference / (T)(height_image.res_x_ + 1) * (T)height_image.res_y_ * height_scale_input;

	const float iso_threshold = 0.5f;
	const float min_edge_length = 1e-6;

	Array2D<T> shape_field_node, thickness_field_node;
	Array2D<int> u_face, v_face, vix_node;

	shape_field_node.initialize(0, 0, shape_image.res_x_, shape_image.res_y_, false);
	thickness_field_node.initialize(0, 0, shape_image.res_x_, shape_image.res_y_, false);
	u_face.initialize(0, 0, shape_image.res_x_, shape_image.res_y_ - 1, false);
	v_face.initialize(0, 0, shape_image.res_x_, shape_image.res_y_, false);
	vix_node.initialize(0, 0, shape_image.res_x_, shape_image.res_y_, false);

	shape_field_node.assignAllValues(0.0f);
	thickness_field_node.assignAllValues(0.0f);
	v_face.assignAllValues(-1);
	u_face.assignAllValues(-1);
	vix_node.assignAllValues(-1);

	for (int j = 0; j < shape_field_node.j_res_; ++j)
	for (int i = 0; i < shape_field_node.i_res_; ++i)
	{
		shape_field_node(i, j) = CLAMP(shape_image.data_(i,j).GetReversedGray(), 0.0f, 1.0f);
	}

	GridUniform2D grid(0, 0, shape_image.res_x_, shape_image.res_y_, 0, 0, 1, 1);

 	SMOOTHING_UNIFORM_2D smoothing;
//	smoothing.SmoothIRepeatJClampCell(grid, shape_field_node, smoothing_repeat);

	for (int j = 0; j < shape_field_node.j_res_; ++j)
		for (int i = 0; i < shape_field_node.i_res_; ++i)
	{
		thickness_field_node(i, j) = CLAMP(height_image.data_(i, j).GetReversedGray(), 0.0f, 1.0f);
	}

	smoothing.SmoothIRepeatJClampCell(grid, thickness_field_node, smoothing_repeat);

	LinkedArray<TV> vertices;
	LinkedArray<TV_INT> triangles;

	// make bottom node vertices
	int offset = 0;
	for (int j = vix_node.j_start_; j <= vix_node.j_end_; ++j)
	for (int i = vix_node.i_start_; i <= vix_node.i_end_; ++i)
	{
		if (shape_field_node(i, j) > iso_threshold)
		{
			const TV2 uv = grid.GetNodeUV(i, j);
			const TV plane_pos = TV(uv.u_, uv.v_, -thickness_field_node(i, j));// back node vertices are initialized as -thickness field and then modified when generating front vertices

			vertices.PushBack() = plane_pos;

			vix_node(i, j) = offset++;
		}
	}

	// make u-edge vertices
	for (int j = u_face.j_start_; j <= u_face.j_end_; ++j)
	for (int i = u_face.i_start_; i <= u_face.i_end_; ++i)
	{
		const T phi0 = shape_field_node(i, j), phi1 = shape_field_node(i, j + 1);

		if (isInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			const T  litho_depth = thickness_field_node(i, j)*(1.0f - alpha) + thickness_field_node(i, j + 1)*alpha;	//TODO: check alpha
			const TV2 uv = grid.GetNodeUV(i, j)*alpha + grid.GetNodeUV(i, j + 1)*(1.0f - alpha);
			const TV plane_pos = cast(uv, -litho_depth);

			vertices.PushBack() = plane_pos;

			u_face(i, j) = offset++;
		}
	}

	// make v-edge vertices
	for (int j = v_face.j_start_; j <= v_face.j_end_; ++j)
	for (int i = v_face.i_start_; i <= v_face.i_end_; ++i)
	{
		const T phi0 = shape_field_node(i, j), phi1 = shape_field_node.getIRepeatedJClamped(i + 1, j);

		if (isInterfacial(phi0, phi1, iso_threshold) == true)
		{
			const T alpha = CLAMP(ABS(phi1 - iso_threshold) / (ABS(phi0 - iso_threshold) + ABS(phi1 - iso_threshold)), min_edge_length, 1.0f - min_edge_length);

			if (alpha < 0.0f || alpha > 1.0f)
			{
				const float temp = 1.0f;

				std::cout << "alpha = " << alpha << " phi0 " << phi0 << " phi1 " << phi1 << std::endl;
			}

			// thickness_field_node(i,j) = 1.0f for inter-facial vertices
			const T  litho_depth = thickness_field_node(i, j)*(1.0f - alpha) + thickness_field_node.getIRepeatedJClamped(i + 1, j)*alpha;
			const TV2 uv = grid.GetNodeUV(i, j)*alpha + grid.GetNodeUV(i + 1, j)*(1.0f - alpha);
			const TV plane_pos = cast(uv, -litho_depth);

			vertices.PushBack() = plane_pos;

			v_face(i, j) = offset++;
		}
	}

	vertices.CopyToArray(surface.vertex_positions_);	// to check triangle edge length.

	// make bottom faces
	for (int j = shape_field_node.j_start_; j <= shape_field_node.j_end_-1; ++j)
	for (int i = shape_field_node.i_start_; i <= shape_field_node.i_end_; ++i)
	{
		int flag(0);

		if (shape_field_node.getIRepeatedJClamped(i, j) > iso_threshold) flag += 1;
		if (shape_field_node.getIRepeatedJClamped(i + 1, j) > iso_threshold) flag += 2;
		if (shape_field_node.getIRepeatedJClamped(i + 1, j + 1) > iso_threshold) flag += 4;
		if (shape_field_node.getIRepeatedJClamped(i, j + 1) > iso_threshold) flag += 8;

		if (flag == 15)		// full
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node.getIRepeatedJClamped(i + 1, j), vix_node(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), vix_node.getIRepeatedJClamped(i + 1, j + 1), vix_node(i, j + 1)).getReversedCW();
		}
		else if (flag == 5) // n0 n2 (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j + 1), v_face(i, j + 1), u_face.getIRepeatedJClamped(i + 1, j)).getReversedCW();
		}
		else if (flag == 10)// n1 n3  (TODO:check ambiguity)
		{
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), u_face.getIRepeatedJClamped(i + 1, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 14)	// except n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), vix_node.getIRepeatedJClamped(i + 1, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j + 1), u_face(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j + 1), v_face(i, j), vix_node.getIRepeatedJClamped(i + 1, j)).getReversedCW();
		}
		else if (flag == 13)	// except n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), u_face.getIRepeatedJClamped(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face.getIRepeatedJClamped(i + 1, j), vix_node.getIRepeatedJClamped(i + 1, j + 1)).getReversedCW();
		}
		else if (flag == 11)	// except n2
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), u_face.getIRepeatedJClamped(i + 1, j), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node.getIRepeatedJClamped(i + 1, j), u_face.getIRepeatedJClamped(i + 1, j)).getReversedCW();
		}
		else if (flag == 7)		// except n3
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node.getIRepeatedJClamped(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), v_face(i, j + 1), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), vix_node.getIRepeatedJClamped(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 1)		// n0 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), v_face(i, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 2)		// n1 only
		{
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), u_face.getIRepeatedJClamped(i + 1, j), v_face(i, j)).getReversedCW();
		}
		else if (flag == 4)		// n2 only
		{
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j + 1), v_face(i, j + 1), u_face.getIRepeatedJClamped(i + 1, j)).getReversedCW();
		}
		else if (flag == 8)		// n3 only
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
		else if (flag == 3) // n0 n1
		{
			triangles.PushBack() = TV_INT(vix_node(i, j), vix_node.getIRepeatedJClamped(i + 1, j), u_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), u_face.getIRepeatedJClamped(i + 1, j), u_face(i, j)).getReversedCW();
		}
		else if (flag == 6) // n1 n2
		{
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), vix_node.getIRepeatedJClamped(i + 1, j + 1), v_face(i, j + 1)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j), v_face(i, j + 1), v_face(i, j)).getReversedCW();
		}
		else if (flag == 12) // n2 n3
		{
			triangles.PushBack() = TV_INT(vix_node.getIRepeatedJClamped(i + 1, j + 1), vix_node(i, j + 1), u_face.getIRepeatedJClamped(i + 1, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), u_face(i, j), u_face.getIRepeatedJClamped(i + 1, j)).getReversedCW();
		}
		else if (flag == 9) // n3 n0
		{
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), vix_node(i, j), v_face(i, j)).getReversedCW();
			triangles.PushBack() = TV_INT(vix_node(i, j + 1), v_face(i, j), v_face(i, j + 1)).getReversedCW();
		}
	}

	// find boundary edges to make side later
	vertices.CopyToArray(surface.vertex_positions_);
	triangles.CopyToArray(surface.triangles_);

	surface.findAdjacentTrianglesOfVertices();
	surface.findEdgeTrianglesOfTriangles();

	LinkedArray<TV2_INT> boundary_edges_linked_array;
	surface.findBoundaryEdges(boundary_edges_linked_array);

	Array1D<TV2_INT> boundary_edges;
	boundary_edges_linked_array.CopyToArray(boundary_edges);

	// make upper node vertices (copy bottom vertices)
	for (int vix = 0; vix < surface.vertex_positions_.num_elements_; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		vertices.PushBack() = getShellX(TV2(bottom_pos.x_, bottom_pos.y_), curve,  - bottom_pos.z_*front_print_thickness + base_thickness);	// bottom_pos.z_ is -thickness_field_node(i, j)
	}

	vertices.CopyToArray(surface.vertex_positions_);

	// reset depth of bottom plane
	for (int vix = 0; vix < offset; ++vix)
	{
		const TV bottom_pos = surface.vertex_positions_[vix];

		surface.vertex_positions_[vix] = getShellX(TV2(bottom_pos.x_, bottom_pos.y_), curve, (T)0); // bottom_pos.z_ is -thickness_field_node(i, j)
	}

	// make upper faces
	for (int tri = 0; tri < surface.triangles_.num_elements_; ++tri)
	{
		triangles.PushBack() = TV_INT(surface.triangles_[tri].i_ + offset, surface.triangles_[tri].k_ + offset, surface.triangles_[tri].j_ + offset);
	}

	// make side from boundary edges
	for (int be = 0; be < boundary_edges.num_elements_; ++be)
	{
		triangles.PushBack() = TV_INT(boundary_edges[be].i_, boundary_edges[be].i_ + offset, boundary_edges[be].j_);
		triangles.PushBack() = TV_INT(boundary_edges[be].j_, boundary_edges[be].i_ + offset, boundary_edges[be].j_ + offset);
	}

	triangles.CopyToArray(surface.triangles_);

	// aspect ratio sampling
	T aspect_ratio_sum = (T)0;
	for (int j = 0; j < thickness_field_node.j_res_ - 1; ++j)
	{
		int i = 0;
		const TV2 uv_00 = TV2((T)i / (T)thickness_field_node.i_res_, (T)j / (thickness_field_node.j_res_ - 1));
		const TV2 uv_10 = TV2((T)(i + 1) / (T)thickness_field_node.i_res_, (T)j / (thickness_field_node.j_res_ - 1));
		const TV2 uv_01 = TV2((T)i / (T)thickness_field_node.i_res_, (T)(j + 1) / (thickness_field_node.j_res_ - 1));
		const TV2 uv_11 = TV2((T)(i + 1) / (T)thickness_field_node.i_res_, (T)(j + 1) / (thickness_field_node.j_res_ - 1));

		const TV2 uv_0h = (uv_00 + uv_01) * (T)0.5;
		const TV2 uv_1h = (uv_10 + uv_11) * (T)0.5;

		const T dx = (getShellX(uv_0h, curve, (T)0) - getShellX(uv_1h, curve, (T)0)).getMagnitude();

		const T dy = (getShellX(uv_00, curve, (T)0) - getShellX(uv_01, curve, (T)0)).getMagnitude();

		aspect_ratio_sum += dx / dy;
	}

	T average_aspect_ratio = aspect_ratio_sum / (T)(thickness_field_node.j_res_ - 1);

	std::cout << "Average aspect ratio (dx / dy) " << average_aspect_ratio << std::endl;

	const T circumference = 2.0f*PI*radius;
	std::cout << "Recommended resolution for aspect ratio = 1 is " << thickness_field_node.i_res_ << " x " << (int)((T)thickness_field_node.j_res_ / average_aspect_ratio) << std::endl;
}

