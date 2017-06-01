#include "GL_TOOLS.h"

namespace GL_TOOLS
{
	void addCircleLines(const TV& center, const T& radius, const int num_segments, LinkedArray<Vector3D<T> >& v_pos_buffer)
	{
		const T dth = (T)2 * (T)PI / (T)num_segments;

		T theta = (T)0;
		for (int i = 0; i < num_segments; i++)
		{
			v_pos_buffer.PushBack() = TV(radius * cos(theta) + center.x_, -radius * sin(theta) + center.y_, center.z_);

			theta += dth;
		}

		v_pos_buffer.PushBack() = TV(radius * cos(0), -radius * sin(0), (T)0);
	}

	void AddCubeEdges(const BOX_3D<T>& cube, LinkedArray<Vector3D<T> >& v_pos_buffer)
	{
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_min_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_min_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_max_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_max_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_min_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_min_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_max_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_max_, cube.z_max_);

		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_min_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_max_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_min_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_max_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_min_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_max_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_min_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_max_, cube.z_max_);

		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_min_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_min_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_min_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_min_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_max_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_min_, cube.y_max_, cube.z_max_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_max_, cube.z_min_);
		v_pos_buffer.PushBack() = TV(cube.x_max_, cube.y_max_, cube.z_max_);
	}

	void AddUniformGridEdges(const GridUniform3D& grid, LinkedArray<Vector3D<T> >& v_pos_buffer)
	{
		const TV cell_origin = grid.getMin();

		for (int k = 0; k <= grid.k_res_; ++k)
		for (int j = 0; j <= grid.j_res_; ++j)
		{
			v_pos_buffer.PushBack() = TV(cell_origin + TV(0.0f, (T)j*grid.dy_, (T)k*grid.dz_));
			v_pos_buffer.PushBack() = TV(cell_origin + TV((T)grid.i_res_*grid.dx_, (T)j*grid.dy_, (T)k*grid.dz_));
		}

		for (int k = 0; k <= grid.k_res_; ++k)
		for (int i = 0; i <= grid.i_res_; ++i)
		{
			v_pos_buffer.PushBack() = TV(cell_origin + TV((T)i*grid.dx_, (T)0.0f, (T)k*grid.dz_));
			v_pos_buffer.PushBack() = TV(cell_origin + TV((T)i*grid.dx_, (T)grid.j_res_*grid.dy_, (T)k*grid.dz_));
		}

		for (int j = 0; j <= grid.j_res_; ++j)
		for (int i = 0; i <= grid.i_res_; ++i)
		{
			v_pos_buffer.PushBack() = TV(cell_origin + TV((T)i*grid.dx_, (T)j*grid.dx_, 0.0f));
			v_pos_buffer.PushBack() = TV(cell_origin + TV((T)i*grid.dx_, (T)j*grid.dx_, (T)grid.k_res_*grid.dz_));
		}
	}
}