#include "GL_OBJECT.h"
#include "GL_WORLD.h"
#include "GL_TOOLS.h"

// GL_OBJECT::GL_OBJECT()
// 	: polygon_mode_face_(GL_FRONT_AND_BACK), polygon_mode_mode_(GL_FILL), draw_surface_(true), draw_edges_(false), draw_short_edges_(false), draw_high_curvature_vertices_(false), draw_slicer_(false)
// {
// 	initializeGLFunctions();
// }

GL_OBJECT::GL_OBJECT(Object* object_input)
	: polygon_mode_face_(GL_FRONT_AND_BACK), polygon_mode_mode_(GL_FILL), draw_surface_(true), draw_edges_(false), draw_short_edges_(false), draw_high_curvature_vertices_(false), draw_slicer_(false), is_locked_(false)
{
	object_ = object_input;

	initializeGLFunctions();

	UpdateSurface();
}

GL_OBJECT::~GL_OBJECT()
{
// 	for (std::list<GL_ELEMENT*>::iterator itr = gl_element_list_.begin(); itr != gl_element_list_.end(); itr++)
// 	{
// 		SAFE_DELETE(*itr);
// 	}
}

void GL_OBJECT::UpdateHighCurvatureVertices(const T& kappa)
{
	object_->surface_.determineVertexMeanCurvatures();

	std::cout << "UpdateHighCurvatureVertices " << kappa << std::endl;

	LinkedArray<TV> positions, colors;

	object_->surface_.copyRenderingDataHighCurvatureVertices(positions, colors, kappa, 1e-2);

	gl_high_curvature_vertices_.InitializeVBO(positions, colors);
}

void GL_OBJECT::UpdateShortEdges(const T& min_length)
{
	LinkedArray<TV> positions, normals;
	object_->surface_.copyRenderingDataShortEdgeTriangles(positions, normals, min_length, (T)1e-4);		//TODO: option normal offset

	gl_short_edge_triangles_.begin(GL_TRIANGLES, USE_POSITION | USE_NORMAL);
	gl_short_edge_triangles_.bindArray(USE_POSITION, positions);
	gl_short_edge_triangles_.bindArray(USE_NORMAL, normals);
	gl_short_edge_triangles_.end();
}

void GL_OBJECT::DrawShortEdgeTriangles(const GL_VIEW& camera)
{
	glPolygonMode(polygon_mode_face_, GL_FILL);

//	gl_short_edge_triangles_.BindShader(&GL_WORLD::GetInstance().shader_programs_.program_phong_);
	gl_short_edge_triangles_.bindShader(&GL_WORLD::GetInstance().shader_programs_.program_p_);
	gl_short_edge_triangles_.applyMaterial(surface_material_);

	gl_short_edge_triangles_.setUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
	gl_short_edge_triangles_.setUniformColor(glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
	gl_short_edge_triangles_.drawArrays();

	gl_short_edge_triangles_.releaseShader();

//	glLineWidth(1.0f);

// 	gl_short_edges_.BindShader(&GL_WORLD::GetInstance().shader_programs_.program_p_);
// 
// 	gl_short_edges_.SetUniformMVPMatrix(camera.GetTransformMatrix()*object_->m_matrix_);
// 	gl_short_edges_.SetUniformColor(glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
// 	gl_short_edges_.DrawArrays();
// 
// 	gl_short_edges_.ReleaseShader();
}

void GL_OBJECT::DrawHighCurvatureVertices(const GL_VIEW& camera)
{
	glPointSize(5.0f);

	gl_high_curvature_vertices_.bindShader(&GL_WORLD::GetInstance().shader_programs_.program_pc_);

	gl_high_curvature_vertices_.setUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
//	gl_high_curvature_vertices_.SetUniformColor(glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
	gl_high_curvature_vertices_.drawArrays();

	gl_high_curvature_vertices_.releaseShader();
}

void GL_OBJECT::update()
{
	if (object_->data_depot_.updated_ == false) return;

	object_->data_depot_.lock();

	is_locked_ = true;

	resetElementList();

	SinglyLinkedList<ColoredParticlesData*> &colored_particles_list(object_->data_depot_.colored_particles_list_);
	for (colored_particles_list.begin(); colored_particles_list.valid(); colored_particles_list.next())
	{
		ColoredParticlesData *colored_particles = colored_particles_list.getItr();
		GL_ELEMENT *new_gl_element = new GL_ELEMENT;
		
		new_gl_element->name_ = std::string("colored_particles");
//		std::cout << new_gl_element->name_ << std::endl;
		new_gl_element->point_size_ = colored_particles->point_size_;
		new_gl_element->begin(GL_POINTS, USE_POSITION | USE_COLOR);
		new_gl_element->bindArray(USE_POSITION, colored_particles->position_);
		new_gl_element->bindArray(USE_COLOR, colored_particles->color_);
		new_gl_element->end();

		gl_element_list_.pushBack(new_gl_element);

//		std::cout << "adding colored particles data gl_element" << std::endl;
	}

	SinglyLinkedList<LinesData*> &lines_data_list(object_->data_depot_.lines_list_);
	for (lines_data_list.begin(); lines_data_list.valid(); lines_data_list.next())
	{
		LinesData *lines = lines_data_list.getItr();
		GL_ELEMENT *new_gl_element = new GL_ELEMENT;

		new_gl_element->name_ = std::string("lines");
//		std::cout << new_gl_element->name_ << std::endl;
		new_gl_element->line_width_ = lines->line_width_;
		new_gl_element->color_ = lines->color_;
		new_gl_element->begin(GL_LINES, USE_POSITION);
		new_gl_element->bindArray(USE_POSITION, lines->vertices_);
		new_gl_element->end();

		gl_element_list_.pushBack(new_gl_element);

//		std::cout << "adding line data gl_element" << std::endl;
	}

	SinglyLinkedList<PhongTrianglesData*> &phong_triangles_list(object_->data_depot_.phong_triangles_list_);
	for (phong_triangles_list.begin(); phong_triangles_list.valid(); phong_triangles_list.next())
	{
		PhongTrianglesData *phong_triangles_ = phong_triangles_list.getItr();

		GL_ELEMENT *new_gl_element = new GL_ELEMENT;

		new_gl_element->name_ = std::string("phong_triangles");
		new_gl_element->rasterize_mode_ = (GLenum)phong_triangles_->rasterize_mode_;

//		std::cout << new_gl_element->name_ << std::endl;
		new_gl_element->begin(GL_TRIANGLES, USE_POSITION | USE_NORMAL);
		new_gl_element->bindArray(USE_POSITION, phong_triangles_->positions_);
		new_gl_element->bindArray(USE_NORMAL, phong_triangles_->normals_);
		new_gl_element->end();

		gl_element_list_.pushBack(new_gl_element);
	}

	object_->data_depot_.unlock();

	object_->data_depot_.updated_ = false;

	is_locked_ = false;
}

void GL_OBJECT::DrawSurface(const GL_VIEW& camera)
{
	// Phong-shade triangular surface
	glPolygonMode(polygon_mode_face_, GL_FILL);
	gl_triangles_.bindShader(&GL_WORLD::GetInstance().shader_programs_.program_phong_);
	gl_triangles_.applyMaterial(surface_material_);
	gl_triangles_.setUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
	gl_triangles_.drawArrays();
	gl_triangles_.releaseShader();

	// draw OOBB
// 	gl_surface_oobb_.BindShader(&GL_WORLD::GetInstance().shader_programs_.program_p_);
// 	gl_surface_oobb_.SetUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
// 	gl_surface_oobb_.DrawArrays();
// 	gl_surface_oobb_.ReleaseShader();
}

void GL_OBJECT::DrawSlicer(const GL_VIEW& camera)
{
	// outline
	gl_slicer_.program_->bind();
	gl_slicer_.setUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
	
	glLineWidth(3);

	static int num_vertices = gl_slicer_.num_vertices_;

	gl_slicer_.drawArrays(num_vertices);
//	gl_slicer_.DrawArrays();

	if (num_vertices >= gl_slicer_.num_vertices_) num_vertices = 0;
	else num_vertices += 10;

	gl_slicer_.program_->release();

	// phi layers
// 	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
// 
// 	gl_phi_layers_.program_->bind();
// 	gl_phi_layers_.SetUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
// 
// 	gl_phi_layers_.DrawArrays();
// 
// 	gl_phi_layers_.program_->release();
}

void GL_OBJECT::DrawEdges(const GL_VIEW& camera)
{
	glLineWidth(1.0f);

	glPolygonMode(polygon_mode_face_, GL_LINE);

	gl_triangles_.bindShader(&GL_WORLD::GetInstance().shader_programs_.program_pn_);

	gl_triangles_.setUniformMVPMatrix(camera.GetWorldViewMatrix()*object_->m_matrix_);
	gl_triangles_.setUniformColor(glm::vec4(0.0f, 0.0f, 1.0f, 1.0f));
	gl_triangles_.setUniformNormalOffset(0.001f);
	
	gl_triangles_.drawArrays();

	gl_triangles_.releaseShader();
}

void GL_OBJECT::Draw(const GL_VIEW& camera)
{
	//glShadeModel(GL_SMOOTH);	// meaningless with fragment phong shader

	gl_triangles_.applyLighting(lighting_parameters_);	//TODO: to world

	if (draw_surface_) DrawSurface(camera);

	if (draw_edges_) DrawEdges(camera);

	if (draw_short_edges_) DrawShortEdgeTriangles(camera);

	if (draw_high_curvature_vertices_) DrawHighCurvatureVertices(camera);

	if (!draw_surface_ && draw_slicer_) DrawSlicer(camera);

	if(is_locked_ == false)
		for (gl_element_list_.begin(); gl_element_list_.valid(); gl_element_list_.next())
		{
//			std::cout << "Drawing GL_Element " << gl_element_list_.getItr()->name_ << std::endl;

			gl_element_list_.getItr()->draw(camera.GetWorldViewMatrix()*object_->m_matrix_);
		}
}

void GL_OBJECT::UpdateSurface()
{
	LinkedArray<TV> v_pos_buffer, v_normal_buffer;
	object_->surface_.copyRenderingData(v_pos_buffer, v_normal_buffer);

	gl_triangles_.begin(GL_TRIANGLES, USE_POSITION | USE_NORMAL);
	gl_triangles_.bindArray(USE_POSITION, v_pos_buffer);
	gl_triangles_.bindArray(USE_NORMAL, v_normal_buffer);
	gl_triangles_.end();
}

void GL_OBJECT::UpdateSurface(StaticTriangularSurface& surface_)
{
	LinkedArray<TV> v_pos_buffer, v_normal_buffer;
	object_->surface_.copyRenderingData(v_pos_buffer, v_normal_buffer);

	gl_triangles_.begin(GL_TRIANGLES, USE_POSITION | USE_NORMAL);
	gl_triangles_.bindArray(USE_POSITION, surface_.vertex_positions_);
	gl_triangles_.bindArray(USE_NORMAL, surface_.vertex_normals_);
	gl_triangles_.end();
}

void GL_OBJECT::UpdateSlicer()
{
//	SAFE_DELETE(slicer_element_);

	// outlines
	{ gl_slicer_.bindShader(&GL_WORLD::GetInstance().shader_programs_.program_pc_);

	Array1D<LINE_SEGMENT> all_contours;

	object_->slicer_.CopyAllContoursToArray(all_contours);
//	object_->slicer_.GenerateGCode();

	Array1D<TV> positions(all_contours.num_elements_ * 2);

	for (int i = 0; i < all_contours.num_elements_; ++i)
	{
		positions[2 * i + 0] = all_contours[i].p0_;
		positions[2 * i + 1] = all_contours[i].p1_;

		positions[2 * i + 0].y_ += 1e-4;
		positions[2 * i + 1].y_ += 1e-4;
	}

	Array1D<glm::vec4> colors(all_contours.num_elements_ * 2);

	for (int i = 0; i < all_contours.num_elements_; ++i)
	{
 		colors[2 * i + 0] = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
 		colors[2 * i + 1] = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);

// 		colors[2 * i + 0] = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
// 		colors[2 * i + 1] = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);
	}

	gl_slicer_.begin(GL_LINES, USE_POSITION | USE_COLOR);
	gl_slicer_.bindArray(USE_POSITION, positions);
	gl_slicer_.bindArray(USE_COLOR, colors);
	gl_slicer_.end(); }

	// phi
// 	LINKED_ARRAY<TV> positions;
// 	LINKED_ARRAY<glm::vec4> colors;
// 
// 	const GRID_UNIFORM_2D &grid(object_->slicer_.grid_);
// 
// 	T z = object_->slicer_.z_min_;
// 	for (int l = 0; l < object_->slicer_.signed_distance_layers_.num_elements_; ++l)
// 	{
// 		for (int i = grid.i_start_; i <= grid.i_end_; ++i)
// 		for (int j = grid.j_start_; j <= grid.j_end_; ++j)
// 		{
// 			if (object_->slicer_.signed_distance_layers_[l](i, j) <= (T)0)
// 			{
// 				const TV2 center = grid.GetCellCenter(i, j);
// 				const T dx = grid.dx_;
// 
// 				positions.PushBack() = TV(center.x_ - dx, center.y_ - dx, z).GetSwapedYZ();		//TODO: don't swap
// 				positions.PushBack() = TV(center.x_ - dx, center.y_ + dx, z).GetSwapedYZ();
// 				positions.PushBack() = TV(center.x_ + dx, center.y_ - dx, z).GetSwapedYZ();
// 				positions.PushBack() = TV(center.x_ + dx, center.y_ - dx, z).GetSwapedYZ();
// 				positions.PushBack() = TV(center.x_ - dx, center.y_ + dx, z).GetSwapedYZ();
// 				positions.PushBack() = TV(center.x_ + dx, center.y_ + dx, z).GetSwapedYZ();
// 				
// 				colors.PushBack() = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
// 				colors.PushBack() = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
// 				colors.PushBack() = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
// 				colors.PushBack() = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
// 				colors.PushBack() = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
// 				colors.PushBack() = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
// 			}
// 		}
// 
// 		z += object_->slicer_.dz_;
// 	}
// 
// 	gl_phi_layers_.BindShader(&GL_WORLD::GetInstance().shader_programs_.program_pc_);
// 	gl_phi_layers_.Begin(GL_TRIANGLES, USE_POSTION | USE_COLOR);
// 	gl_phi_layers_.BindArray(USE_POSTION, positions);
// 	gl_phi_layers_.BindArray(USE_COLOR, colors);
// 	gl_phi_layers_.End();
}

void GL_OBJECT::UpdateAMR()
{
	LinkedArray<TV> v_pos_buffer;

	const GridAdaptive3D &amr_grid(object_->amr_grid_);
	const GridUniform3D &grid(object_->amr_grid_.grid_);
	const Array3D<bool> &interfacial(object_->amr_grid_.interfacial_);
	const Array3D<int>  &level(object_->amr_grid_.level_);

	// draw intersecting cells
	for (int k = grid.k_start_; k <= grid.k_end_; ++k)
	for (int j = grid.j_start_; j <= grid.j_end_; ++j)
	for (int i = grid.i_start_; i <= grid.i_end_; ++i)
	{
		const int l = amr_grid.level_(i, j, k);

		if (amr_grid.interfacial_(i, j, k) == true)
		{
			GL_TOOLS::AddCubeEdges(BOX_3D<T>(grid.getCellCenter(i, j, k), grid.dx_), v_pos_buffer);
		}

		if (l > 0)
		{
			GridUniform3D subgrid;
			subgrid.initialize(grid.getCellMinMax(i, j, k), POW_OF_TWO(l));

			for (int n = subgrid.k_start_; n <= subgrid.k_end_; ++n)
			for (int m = subgrid.j_start_; m <= subgrid.j_end_; ++m)
			for (int l = subgrid.i_start_; l <= subgrid.i_end_; ++l)
			{
				if (amr_grid.inter_.values_[amr_grid.first_ix_(i, j, k) + subgrid.get1DIndex(l, m, n)] == true)
				{
					GL_TOOLS::AddCubeEdges(subgrid.getCellMinMax(l, m, n), v_pos_buffer);
				}
			}
		}


		//		if (interfacial(i, j, k) == true)
		//			GL_TOOLS::AddCubeEdges(BOX_3D<T>(grid.getCellCenter(i,j,k), grid.dx_), v_pos_buffer);
		// 		if (level(i,j,k) > 0)
		// 		{
		// 			int res;
		// 			if (level(i, j, k) == 1) res = 2;
		// 			else if (level(i, j, k) == 2) res = 4;
		// 
		// 			const GRID_UNIFORM_3D subgrid(grid.getCellMinMax(i, j, k), TV_INT(res, res, res));
		// 
		// 			GL_TOOLS::AddUniformGridEdges(subgrid, v_pos_buffer);
		// 		}
	}

	GL_TOOLS::AddCubeEdges(grid.getMimMax(), v_pos_buffer);

	// draw all cells
	//	GL_TOOLS::AddUniformGridEdges(object_->amr_grid_.grid_, v_pos_buffer);

	gl_surface_oobb_.begin(GL_LINES, USE_POSITION);
	gl_surface_oobb_.bindArray(USE_POSITION, v_pos_buffer);
	gl_surface_oobb_.end();
}