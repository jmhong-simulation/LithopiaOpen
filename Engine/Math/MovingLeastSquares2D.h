// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "MATRIXMXN.h"
#include "MATRIXNXN.h"

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/VectorND.h"
#include "DataStructure/Vector3D.h"
#include "DataStructure/Array1D.h"
#include "MLSPointConstraint2D.h"
#include "MLSLineConstraint2D.h"
#include <vector>

#define MLS_VALUE						0
#define MLS_VALUE_AND_NORMAL			1
#define MLS_DIFFUSIVE_NORMAL_COMPONENT	2
#define MLS_DIFFUSIVE_NORMAL			3

class MovingLeastSquares2D
{
public:
	const static int dimension_ = 2;
	bool delete_constraints_when_destroyed_;

private:
	T epsilon, epsilon_squared, one_over_epsilon_squared;
	int size_of_basis, degree_of_basis;

public:
	std::vector<MLSPointConstraint2D*> point_constraints_;
	std::vector<MLSLineConstraint2D*> line_constraints_;

public:
	MovingLeastSquares2D(const int& degree_of_basis_input = 0, const T& epsilon_input = (T)0.001, const int& number_of_constraints_input = 0)
	{
		delete_constraints_when_destroyed_ = true;

		SetDegreeOfBasis(degree_of_basis_input);
		SetEpsilon(epsilon_input);

		point_constraints_.reserve(number_of_constraints_input);
		for(int i=0;i<number_of_constraints_input;i++) point_constraints_.push_back(new MLSPointConstraint2D(TV2(), (T)0));
	}

	~MovingLeastSquares2D()
	{
		if(delete_constraints_when_destroyed_ == true)
		{
			const int number_of_constraints = (int)point_constraints_.size();
			for(int i = 0; i < number_of_constraints; i++) 
				delete point_constraints_[i];
		}
	}

	inline void Reset()
	{
		int number_of_constraints = (int)point_constraints_.size();
		for(int i = 0; i < number_of_constraints; i++) 
			delete point_constraints_[i];

		point_constraints_.clear();
	}

public:
	void AddConstraint(const TV2& position,const T& value)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position,value));
	}

	void AddConstraint(const TV2& position,const T& value,const int& type)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position, value, type));
	}

	void AddConstraint(const TV2& position, const T& value, const TV2& normal)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position, value, normal));
	}

	void AddConstraint(const TV2& position, const T& value, const TV2& normal_input, const int& type)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position, value, normal_input, type));
	}

	void AddConstraint(const TV2& position, const T& value, const TV2& normal_input, const int& type,const T& weight)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position, value, normal_input, type, weight));
	}

	void AddConstraint(const TV2& position, const T& value,const TV2& normal_input, const int& type,const T& weight,const T& divergence_input)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position, value,normal_input, type, weight, divergence_input));
	}

	void AddConstraint(const TV2& position, const T& value, const TV2& normal_input, const int& type, const T& weight, const T& divergence_input, const Array1D<TV2>& vector_gradient_input)
	{
		point_constraints_.push_back(new MLSPointConstraint2D(position, value,normal_input, type, weight, divergence_input, vector_gradient_input));
	}

	void AddConstraint(MLSPointConstraint2D *mls_constraint_)
	{
		point_constraints_.push_back(mls_constraint_);
	}

	void addLineConstraint(const TV2& p0, const TV2& p1, const T& v0, const T& v1, const T& w)
	{
		line_constraints_.push_back(new MLSLineConstraint2D(p0, p1, v0, v1, w));
	}

	void SetDegreeOfBasis(const int& degree_of_basis_input);
	void SetEpsilon(const T& epsilon_input)
	{
		epsilon = epsilon_input;
		epsilon_squared = epsilon * epsilon;
		one_over_epsilon_squared = (T)1 / epsilon_squared;
	}

	VectorND<T> GetBasis(const TV2& position);
	VectorND<D> GetBasis(const DV2& position);
	VectorND<T> GetBasisDerivative(const TV2& position, const int& d);// dimension should start from 0
	VectorND<D> GetBasisDerivative(const DV2& position, const int& d);// dimension should start from 0
	VectorND<T> GetBasisSecondDerivative(const TV2& position, const int& d);// dimension should start from 0

	T getScalar(const TV2& x);// scalar interpolation
	//T GetScalar_Old(const TV2& x);// scalar interpolation
	//TV2 Get_Gradient(const TV2& x);// gradient of the interpolated scalar	
	//TV2 Get_Diffusive_Second_Derivatives(const TV2& x);
	//TV2 Get_Polygonal_Second_Derivatives(const TV2& x);
	//TV2 Get_Numerical_Second_Derivatives(const TV2& x,const T& dx);
	//TV2 Get_Diffuse_Gradient(const TV2& x);
	//TV2 Get_Polynomial_Gradient(const TV2& x);
	//TV2 Get_Numerical_Gradient(const TV2& x,const T& dx);
	//T Get_Laplacian(const TV2& x);
	//T Get_Diffuse_Laplacian(const TV2& x);// laplacian of the interpolated scalar
	//T Get_Polynomial_Laplacian(const TV2& x);
	//T Get_Numerical_Laplacian(const TV2& x,const T& dx);
	//void Get_Poisson_Coefficients(const TV2& x,VECTOR_ND<T>& coefficients);

	TV2 GetVector(const TV2& x);// simple vector interpolation
	//TV2 Get_Vector_From_Normal_Components(const TV2& x);// normal component constraints
	//TV2 Get_Vector_PDF1(const TV2& x);// div-free condition at all ctr points
	//TV2 Get_Vector_PDF2(const TV2& x);// div-free condition at all ctr points and x
	//TV2 Get_Vector_DIV(const TV2& x,const T& divergence=(T)0);// only one div-free constraint at x (moving div-free constraint)
	//TV2 Get_Vector_DIV_W3(const TV2& x,const T& divergence=(T)0);// only one div-free constraint at x with higher-order w
	//TV2 Get_Vector_DIV_W4(const TV2& x,const T& divergence=(T)0);// only one div-free constraint at x with higher-order w
	//TV2 Get_Vector_DIV_And_Normal_Components(const TV2& x);// normal component constraints and divergence free condition at x
	//TV2 Get_Vector_Numerical_DIV(const TV2& x,const T& dx);// div-free at x by numerical differentiation
	//T Get_Divergence(const TV2& x);

	//void Get_CMIP1(const TV2& x,T& value,TV2& derivatives);//Cubic scalar interpolation
	//void Get_CMIP2(const TV2& x,T& value,TV2& derivatives);//Cubic scalar interpolation
	//void Get_CMIP(const TV2& x,TV2& value,ARRAY<TV2>& derivatives,
	//			  const bool use_value_constraints,const bool use_pdf_constraint,const bool use_indirect_derivatives,
	//			  const T value_constraint_weight,const T derivatives_weight,const T pdf_constraint_weight);//Cubic vector interpolation with PDF constraint

	//T Get_Numerical_Divergence(const TV2& x,const T& dx);
	//TV2 Get_Numerical_Jacobian(const TV2& x,const TV2& velocity,const T& dx=1e-6);

	//void Print_Constraints();
	//void Report_Scalar_Error();//Errors at constraint positions
	//void Report_Vector_Error();//Errors at constraint positions

	////SPH implementation. This is temporary. Define SPH class later
	//T Get_SPH_Scalar(const TV2& x);

};//End of class MLS

