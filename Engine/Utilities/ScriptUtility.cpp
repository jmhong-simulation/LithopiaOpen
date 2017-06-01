#include "ScriptUtility.h"

namespace ScriptUtility
{
	void initialize(ScriptBlock& sb, GridUniform3D& grid)
	{
		const Vector3D<float>	min = sb.getBlock("GridUniform3D").getValue("min", Vector3D<float>());
		const Vector3D<float>	max = sb.getBlock("GridUniform3D").getValue("max", Vector3D<float>());
		const Vector3D<int>		ix_start = sb.getBlock("GridUniform3D").getValue("ix_start", Vector3D<int>());
		const Vector3D<int>		res = sb.getBlock("GridUniform3D").getValue("res", Vector3D<int>());

		grid.initialize(ix_start, res, min, max);
	}

	BOX_3D<float> initializeBox3D(ScriptBlock& sb)
	{
		const Vector3D<float>	min = sb.getBlock("Box3D").getValue("min", Vector3D<float>());
		const Vector3D<float>	max = sb.getBlock("Box3D").getValue("max", Vector3D<float>());

		return BOX_3D<float>(min, max);
	}
}