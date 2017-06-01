// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "PARAMETRIC_CURVE_SEGMENT.h"
#include <vector>
#include "DataStructure/Array1D.h"
#include "Utilities/ScriptParameterList.h"
#include "DataStructure/NonlinearMapping1D.h"

class ParametricCurve
{
public:
	Array1D<T> t_domain_;
	Array1D<ParametricCurveSegment> segments_;
	Array1D<NonlinearMapping1D> u2t_;

	int mapping_mode_;

public:
	void initializeFromScript(ScriptParameterList& par, const T x_scale, const T y_scale, const T z_scale, const int _mapping_mode)
	{
//		ScriptParameterList& script = par.script_reader_;

		// check the number of curve segments
		int s = 0;
		while (true)
		{
			std::string curveN_x = std::string("curve") + std::to_string(s) + std::string("_x");

			if (par.isValid(curveN_x) == false) break;

			s++;
		}

		segments_.initialize(s);

		for (int s = 0; s < segments_.num_elements_; ++s)
		{
			std::string curveN_x = std::string("curve") + std::to_string(s) + std::string("_x");
			std::string curveN_y = std::string("curve") + std::to_string(s) + std::string("_y");
			std::string curveN_z = std::string("curve") + std::to_string(s) + std::string("_z");

			Vector4D<T> curve_x, curve_y, curve_z;
			curve_x = par.getValue(curveN_x, Vector4D<T>());
			curve_y = par.getValue(curveN_y, Vector4D<T>());
			curve_z = par.getValue(curveN_z, Vector4D<T>());

			curve_x *= x_scale;
			curve_y *= y_scale;
			curve_z *= z_scale;

			segments_[s].x_.InitializeZeroOneDomain(curve_x);// radius direction
			segments_[s].y_.InitializeZeroOneDomain(curve_y);
			segments_[s].z_.InitializeZeroOneDomain(curve_z);// height
		}

		normalizeTRange();

		mapping_mode_ = _mapping_mode;

		if (mapping_mode_ > 0)
			prepareForOptimizedMapping();
	}

	void normalizeTRange()
	{
		t_domain_.initialize(segments_.num_elements_, (T)0);

		const int num_length_segments = 100000;	// accuracy parameter

		T length = 0;

		for (int s = 0; s < segments_.num_elements_; s++)
			length += segments_[s].getLength(0, 1, num_length_segments);

		for (int s = 0; s < segments_.num_elements_; ++s)
		{
			t_domain_[s] = segments_[s].getLength(0, 1, num_length_segments) / length;
		}

		for (int s = 1; s < t_domain_.num_elements_; ++s)
		{
			t_domain_[s] = t_domain_[s - 1] + t_domain_[s];
		}
	}

	void prepareForOptimizedMapping()
	{
		const int num_samples = 100;	// TODO: make an option
		
		u2t_.initialize(segments_.num_elements_);

		for (int s = 0; s < segments_.num_elements_; s++)
		{
			u2t_[s].initialize(num_samples, 0, 1);

			if(mapping_mode_ == 1)
				u2t_[s].normalizeByLength(segments_[s]);
			else if(mapping_mode_ == 2)
				u2t_[s].normalizeByLengthOverX(segments_[s]);
			else
			{
				std::cout << "Unknown mapping mode " << mapping_mode_ << std::endl;
				exit(1);
			}
		}
	}

	int findSegment(const T t) const
	{
		if (t <= t_domain_[0]) return 0;
		else
		{
			for (int s = 1; s < t_domain_.num_elements_; ++s)
			{
				if (t_domain_[s-1] <= t && t <= t_domain_[s]) return s;		// curve segments need to continuous
			}
		}

		return -1;
	}

	T getLength(const T t0, const T t1) const
	{
		const TV x0 = getPosition(t0), x1 = getPosition(t1);

		return (T)(x1 - x0).getMagnitudeDouble();
	}

	T getLength() const
	{
		const T t0 = 0, t1 = 1;
		const int num_length_segments = 100000;	// accuracy parameter

		const T dt = (t1 - t0) / (T)num_length_segments;

		T length = (T)0;
		for (T t = t0; t < t1; t += dt)
		{
			length += getLength(t0, t0 + dt);
		}

		return length;
	}

	Vector3D<T> getPosition(const T t) const
	{
		const int s = findSegment(t);
		
		if (s == -1)
		{
			std::cout << "ParametricCurve::getPosition Cannot find t " << std::endl;
			exit(1);

			return Vector3D<T>();
		}
		else if (s == 0)
		{
			const T seg_length = t_domain_[s];
			const T t_seg = t / t_domain_[s];

			if(mapping_mode_ == 0)
				return segments_[s].getPosition(t_seg);
			else
				return segments_[s].getPosition(u2t_[s].getT(t_seg));

		}
		else // s > 0
		{
			const T seg_length = t_domain_[s] - t_domain_[s - 1];
			const T t_seg = (t - t_domain_[s - 1]) / seg_length;

			if (mapping_mode_ == 0)
				return segments_[s].getPosition(t_seg);
			else
				return segments_[s].getPosition(u2t_[s].getT(t_seg));
		}
	}	
};
