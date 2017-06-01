#include "SLICER.h"

void SLICER::SliceOutlines(const StaticTriangularSurface& surface)
{
	//num_slices_ = 900;
	const T z_min = surface.OOBB_.z_min_, z_max = surface.OOBB_.z_max_, dz = (z_max - z_min) / (T)num_slices_;

	z_min_ = z_min;	// for rendering
	dz_ = dz;		// for rendering

	Array1D<LinkedArray<int>> z_sorted_triangles(num_slices_);

	for (int tri_ix = 0; tri_ix < surface.triangles_.num_elements_; ++tri_ix)
	{
		StaticTriangle triangle;
		surface.getTriangle(tri_ix, triangle);
		
		T tri_z_min, tri_z_max;
		triangle.getZMinMax(tri_z_min, tri_z_max);

		for (int z_ix = floorf((tri_z_min - z_min) / dz); z_ix < ceilf((tri_z_max - z_min) / dz); ++z_ix)
			z_sorted_triangles[z_ix].PushBack() = tri_ix;
	}

	Array1D<Array1D<int>> z_sorted_triangles_array(num_slices_);

	for (int z_ix = 0; z_ix < num_slices_; ++z_ix)
	{
		z_sorted_triangles[z_ix].CopyToArray(z_sorted_triangles_array[z_ix]);
	}

	Array1D<LinkedArray<LINE_SEGMENT>> &z_sorted_contours(z_sorted_contours_);	// TODO: initialize
	z_sorted_contours.initialize(num_slices_);

	for (int z_ix = 0; z_ix < num_slices_; ++z_ix)
	{
		PLANE z_plane; z_plane.n_ = TV(0, 0, 1); z_plane.p_ = TV(0, 0, z_min + (T)z_ix*dz);
		for (int tri_ix = 0; tri_ix < z_sorted_triangles_array[z_ix].num_elements_; ++tri_ix)
		{
			StaticTriangle triangle;
			surface.getTriangle(z_sorted_triangles_array[z_ix][tri_ix], triangle);

			SliceTriangle(triangle, z_plane, z_sorted_contours[z_ix]);
		}
	}

	for (int z_ix = 0; z_ix < num_slices_; ++z_ix)
	{
		SortContours(z_sorted_contours[z_ix]);
	}

// 	Array1D<LINE_SEGMENT> temp_contour;
// 	y_sorted_contours_[0].CopyToArray(temp_contour);
// 	for (int i = 0; i < temp_contour.num_elements_; ++i)
// 	{
// 		std::cout << "G1 " << temp_contour[i].p1_ << std::endl;
// 	}
}

void SLICER::SliceDLP(const StaticTriangularSurface& surface)
{
	grid_.Initialize(0, 0, 1000, 1000, surface.OOBB_.x_min_, surface.OOBB_.z_min_, surface.OOBB_.x_max_, surface.OOBB_.z_max_);	//TODO z to y
	grid_.Enlarge(2);

	signed_distance_layers_.initialize(z_sorted_contours_.num_elements_);

	for (int l = 0; l < signed_distance_layers_.num_elements_; ++l)
	{
		grid_.InitializeCellArray(signed_distance_layers_[l]);

		signed_distance_layers_[l].assignAllValues(0);
	}

	for (int l = 0; l < signed_distance_layers_.num_elements_; ++l)
	{
		Array1D<LINE_SEGMENT> contour;
		z_sorted_contours_[l].CopyToArray(contour);

		for (int s = 0; s < contour.num_elements_; ++s)
		{
			const TV2_INT start_ix = grid_.Cell(TV2(contour[s].p0_.x_, contour[s].p0_.z_));		//TODO: check y - z
			const TV2_INT end_ix = grid_.Cell(TV2(contour[s].p1_.x_, contour[s].p1_.z_));

			const T y = contour[s].p0_.y_;

			const int i_min = MIN2(start_ix.i_, end_ix.i_);
			const int i_max = MAX2(start_ix.i_, end_ix.i_);
			const int j_min = MIN2(start_ix.j_, end_ix.j_);
			const int j_max = MAX2(start_ix.j_, end_ix.j_);

			for (int i = i_min; i <= i_max; ++i)
			for (int j = j_min; j <= j_max; ++j)
			{
				const TV2 cell_center = grid_.GetCellCenter(i, j);
				const TV v = TV(cell_center.x_, y, cell_center.y_);

				const T distance = contour[s].getDistance(v);

				if(distance <= grid_.dx_) signed_distance_layers_[l](i, j) = -1;
			}
		}

		Array2D<T> &phi = signed_distance_layers_[l];

		phi.FloodFill(0, 0, -1, 1);

		for(int j = phi.j_start_; j <= phi.j_end_; ++j)
		for (int i = phi.i_start_; i <= phi.i_end_; ++i)
		{
			if (phi(i, j) == 0.0f) phi(i, j) = -1;
		}			
	}
}

const T e_coeff = 0.062336f;

void SLICER::GenerateRaft(std::ofstream& file, T& extrusion) const
{
	const T x_min = 40, x_max = 90, y_min = 40, y_max = 90;
	const T dx = 0.8f;

	T current_z = 0.3f;

	TV pos;
	pos = TV(x_min, y_min, current_z);
	file << "G0 F2400 X" << pos.x_ << " Y" << pos.y_ << " Z" << current_z << std::endl;

	for (T x = x_min; x < x_max; x+=2.0*dx)
	{
		TV new_pos(x, y_min, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;

		new_pos = TV(x, y_max, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;

		new_pos = TV(x + dx, y_max, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;

		new_pos = TV(x + dx, y_min, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;
	}

	current_z += 0.2f;
	pos = TV(x_min, y_min, current_z);
	file << "G0 F2400 X" << pos.x_ << " Y" << pos.y_ << " Z" << current_z << std::endl;

	for (T y = y_min; y < y_max; y += 2.0*dx)
	{
		TV new_pos(x_min, y, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;

		new_pos = TV(x_max, y, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;

		new_pos = TV(x_max, y + dx, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;

		new_pos = TV(x_min, y + dx, current_z);
		extrusion += (new_pos - pos).getMagnitude()*e_coeff;
		pos = new_pos;
		file << "G1 F3600 X" << new_pos.x_ << " Y" << new_pos.y_ << " E" << extrusion << std::endl;
	}
}

void SLICER::GenerateGCode() const
{
	std::cout << "# Writing " << "Test.gcode" << std::endl;
	std::ofstream file("Test.gcode");

	file << "; Generated with Cura_SteamEngine 13.12" << std::endl;
	file << "M140 S95.000000" << std::endl;
	file << "M109 T0 S240.000000" << std::endl;
	file << "T0" << std::endl;
	file << "M190 S95.000000" << std::endl;
	file << "; Sliced at : Sat 14 - 03 - 2015 17 : 07 : 20" << std::endl;
	file << "; Basic settings : Layer height : 0.2 Walls : 1 Fill : 0" << std::endl;
	file << "; Print time : #P_TIME#" << std::endl;
	file << "; Filament used : #F_AMNT#m #F_WGHT#g" << std::endl;
	file << "; Filament cost : #F_COST#" << std::endl;
	file << "G29" << std::endl;
	file << "G1 Y169 F8000" << std::endl;
	file << "G92 E0" << std::endl;
	file << "G1 E20 F300" << std::endl;
	file << "G92 E0" << std::endl;
	file << "; Put printing message on LCD screen" << std::endl;
	file << "M117 Printing..." << std::endl;
	file << std::endl;

	file << ";Layer Count: " << z_sorted_contours_.num_elements_ << std::endl;
	file << "M107" << std::endl;

	LinkedArray<LINE_SEGMENT> all_contours_temp;

	T bed_z_base = 0.3f + 0.2*2.0f;	// two raft layers
	T bed_dz = 0.2f;

	TV prev(0.2, 0.2, bed_z_base);
	T extrusion = 0.0f;

	// move to the first position of the layer (translate in z direction)
	file << "G0 F2400 X" << prev.x_ << " Y" << prev.y_ << " Z" << bed_z_base << std::endl;
	file << "G1 F2400 E0.00000" << std::endl;

	GenerateRaft(file, extrusion);

	//	exit(1);

	for (int y_ix = 0; y_ix < z_sorted_contours_.num_elements_; ++y_ix)
		//	for (int y_ix = 0; y_ix < 3; ++y_ix)
	{
		file << "; LAYER:" << y_ix << std::endl;

		Array1D<LINE_SEGMENT> one_contours_temp;
		z_sorted_contours_[y_ix].CopyToArray(one_contours_temp);

		//transform all coordinates
		for (int i = 0; i < one_contours_temp.num_elements_; ++i)
		{
			one_contours_temp[i].p0_.swapYZ();
			one_contours_temp[i].p1_.swapYZ();

			//			one_contours_temp[i].p0_.x_ = -one_contours_temp[i].p0_.x_;
			one_contours_temp[i].p0_.x_ += 120;
			one_contours_temp[i].p0_.y_ += 75;

			//			one_contours_temp[i].p1_.x_ = -one_contours_temp[i].p1_.x_;
			one_contours_temp[i].p1_.x_ += 120;
			one_contours_temp[i].p1_.y_ += 75;
		}

		if(one_contours_temp.num_elements_ > 0)
		{
			file << "G0 X" << one_contours_temp[0].p0_.x_ << " Y" << one_contours_temp[0].p0_.y_ << " Z" << bed_z_base << std::endl;
		}

		for (int i = 0; i < one_contours_temp.num_elements_; ++i)
		{
			all_contours_temp.PushBack() = one_contours_temp[i];
			extrusion += (one_contours_temp[i].p1_ - one_contours_temp[i].p0_).getMagnitude()*e_coeff;

			file << "G1 F1200 X" << one_contours_temp[i].p1_.x_ << " Y" << one_contours_temp[i].p1_.y_ << " E" << extrusion << std::endl;
		}

		bed_z_base += bed_dz;
	}

	file << "; End GCode" << std::endl;
	file << "M104 S0; extruder heater off" << std::endl;
	file << "M140 S0; heated bed heater off(if you have it)" << std::endl;
	file << "G92 E0" << std::endl;
	file << "G1 E - 3 F500" << std::endl;
	file << "G92 E0" << std::endl;
	file << "G28 Y" << std::endl;
	file << "M84; steppers off" << std::endl;
	file << "G90; absolute positioning" << std::endl;

	file.close();
}

void SLICER::CopyAllContoursToArray(Array1D<LINE_SEGMENT>& all_contours) const
{
	LinkedArray<LINE_SEGMENT> all_contours_temp;

	for (int y_ix = 0; y_ix < z_sorted_contours_.num_elements_; ++y_ix)
	{
		Array1D<LINE_SEGMENT> one_contours_temp;
		z_sorted_contours_[y_ix].CopyToArray(one_contours_temp);

		for (int i = 0; i < one_contours_temp.num_elements_; ++i)
		{
			all_contours_temp.PushBack() = one_contours_temp[i];
		}
	}

	all_contours_temp.CopyToArray(all_contours);
}

const TV SLICER::GetSignedDistanceWeightedAverage(const T& phi0, const T& phi1, const TV& v0, const TV& v1)
{
	const T alpha = ABS(phi0) / (ABS(phi0) + ABS(phi1));

	return v0 * ((T)1 - alpha) + v1 * alpha;
}

void SLICER::SliceTriangle(const StaticTriangle triangle, const PLANE plane, LinkedArray<LINE_SEGMENT>& contour)
{
	// signed distances of triangle vertices to the plane
	const T sd_v0 = plane.GetSignedDistance(triangle.v0_);
	const T sd_v1 = plane.GetSignedDistance(triangle.v1_);
	const T sd_v2 = plane.GetSignedDistance(triangle.v2_);

	// TODO: use epsilon
	// TODO: use triangle normal to determine contour segment direction
	if (sd_v0 == 0.0f && sd_v1 == 0.0f && sd_v2 == 0.0f)		// if this triangle is on the plane 
	{
		contour.PushBack().Initialize(triangle.v0_, triangle.v1_);
		contour.PushBack().Initialize(triangle.v1_, triangle.v2_);
		contour.PushBack().Initialize(triangle.v2_, triangle.v0_);
	}
	else if (sd_v0 == 0.0f && sd_v1 == 0.0f)		// if edge 2 is on the plane
	{
		contour.PushBack().Initialize(triangle.v0_, triangle.v1_);
	}
	else if (sd_v1 == 0.0f && sd_v2 == 0.0f)		// if edge 0 is on the plane
	{
		contour.PushBack().Initialize(triangle.v1_, triangle.v2_);
	}
	else if (sd_v2 == 0.0f && sd_v0 == 0.0f)		// if edge 1 is on the plane
	{
		contour.PushBack().Initialize(triangle.v2_, triangle.v0_);
	}	
	else if ((sd_v0 < 0.0f && (sd_v1 > 0.0f && sd_v2 > 0.0f)) || (sd_v0 > 0.0f && (sd_v1 < 0.0f && sd_v2 < 0.0f)))		// edge 1 and 2 are intersecting
	{
		const TV e1 = GetSignedDistanceWeightedAverage(sd_v0, sd_v1, triangle.v0_, triangle.v1_);
		const TV e2 = GetSignedDistanceWeightedAverage(sd_v0, sd_v2, triangle.v0_, triangle.v2_);

		contour.PushBack().Initialize(e1, e2);
	}
	else if ((sd_v1 < 0.0f && (sd_v0 > 0.0f && sd_v2 > 0.0f)) || (sd_v1 > 0.0f && (sd_v0 < 0.0f && sd_v2 < 0.0f)))		// edge 0 and 2 are intersecting
	{
		const TV e2 = GetSignedDistanceWeightedAverage(sd_v1, sd_v0, triangle.v1_, triangle.v0_);
		const TV e0 = GetSignedDistanceWeightedAverage(sd_v1, sd_v2, triangle.v1_, triangle.v2_);

		contour.PushBack().Initialize(e2, e0);
	}
	else if ((sd_v2 < 0.0f && (sd_v0 > 0.0f && sd_v1 > 0.0f)) || (sd_v2 > 0.0f && (sd_v0 < 0.0f && sd_v1 < 0.0f)))		// edge 0 and 1 are intersecting
	{
		const TV e1 = GetSignedDistanceWeightedAverage(sd_v2, sd_v0, triangle.v2_, triangle.v0_);
		const TV e0 = GetSignedDistanceWeightedAverage(sd_v2, sd_v1, triangle.v2_, triangle.v1_);

		contour.PushBack().Initialize(e1, e0);
	}
	else if ((sd_v0 == 0.0f && (sd_v1 > 0.0f && sd_v2 < 0.0f)) || (sd_v0 == 0.0f && (sd_v1 < 0.0f && sd_v2 > 0.0f)))		// v0 and edge 0 are intersecting
	{
		const TV e0 = GetSignedDistanceWeightedAverage(sd_v1, sd_v2, triangle.v1_, triangle.v2_);

		contour.PushBack().Initialize(e0, triangle.v0_);
	}
	else if ((sd_v1 == 0.0f && (sd_v0 > 0.0f && sd_v2 < 0.0f)) || (sd_v1 == 0.0f && (sd_v0 < 0.0f && sd_v2 > 0.0f)))		// v1 and edge 1 are intersecting
	{
		const TV e1 = GetSignedDistanceWeightedAverage(sd_v0, sd_v2, triangle.v0_, triangle.v2_);

		contour.PushBack().Initialize(e1, triangle.v1_);
	}
	else if ((sd_v2 == 0.0f && (sd_v0 > 0.0f && sd_v1 < 0.0f)) || (sd_v2 == 0.0f && (sd_v0 < 0.0f && sd_v1 > 0.0f)))		// v2 and edge 2 are intersecting
	{
		const TV e2 = GetSignedDistanceWeightedAverage(sd_v0, sd_v1, triangle.v0_, triangle.v1_);

		contour.PushBack().Initialize(e2, triangle.v2_);
	}
	else if (sd_v0 == 0.0f)
	{
		// ignore
	}
	else if (sd_v1 == 0.0f)
	{
		// ignore
	}
	else if (sd_v2 == 0.0f)
	{
		// ignore
	}
}

const int SLICER::FindConnectedContour(const int i_start, Array1D<LINE_SEGMENT>& contour_temp, const Array1D<bool>& contour_check, LinkedArray<LINE_SEGMENT>& contour_, const double ep_sqr)
{
	LINE_SEGMENT current_contour = contour_temp[i_start];

	for (int j = 0; j < contour_temp.num_elements_; ++j)
	{
		if (contour_check[j] == true) continue;

		LINE_SEGMENT candidate_contour = contour_temp[j];

		if ((current_contour.p1_ - candidate_contour.p0_).getSqrMagnitudeDouble() < ep_sqr)
		{
// 			contour_check[j] = true;
// 			contour_.PushBack() = LINE_SEGMENT(candidate_contour.p0_, candidate_contour.p1_);
			
			return j;
		}

		if ((current_contour.p1_ - candidate_contour.p1_).getSqrMagnitudeDouble() < ep_sqr)
		{

			contour_temp[j] = LINE_SEGMENT(candidate_contour.p1_, candidate_contour.p0_);		// swap order

// 			contour_check[j] = true;
// 			contour_.PushBack() = ;

			return j;
		}
	}

	return -1;
}

void SLICER::SortContours(LinkedArray<LINE_SEGMENT>& contour)
{
	if (contour.num_elements_ == 0) return;

	const double ep_sqr = POW2(1e-4);

	Array1D<LINE_SEGMENT> contour_temp;
	contour.CopyToArray(contour_temp);

	Array1D<bool> contour_check(contour_temp.num_elements_, false);

	contour.Reset();	// Empty contour so that it can receive sorted contours

	int i = 0, i_flag = 0;
	while (true)
	{
		if (contour_check[i] == true)
		{
			++i;
			if (i == contour_temp.num_elements_) break;
			else continue;
		}

		LINE_SEGMENT current_contour = contour_temp[i];

		contour_check[i] = true;
		contour.PushBack() = current_contour;		// start of a connected contour

		const int j = FindConnectedContour(i, contour_temp, contour_check, contour, ep_sqr);

		if (j == -1)	// lost connection
		{
			i = i_flag;
		}
		else
		{
			i_flag = i;
			i = j;
		}

		if (i == contour_temp.num_elements_) break;
	}
}