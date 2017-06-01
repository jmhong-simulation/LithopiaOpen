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

#include <vector>

#define MLS_VALUE						0
#define MLS_VALUE_AND_NORMAL			1
#define MLS_DIFFUSIVE_NORMAL_COMPONENT	2
#define MLS_DIFFUSIVE_NORMAL			3

class MOVING_LEAST_SQUARES
{
public:
	const static int dimension_ = 3;
	bool delete_constraints_when_destroyed_;

public:
	class MLS_CONSTRAINT
	{
	public:
		int type_;//type of given constraint (0:value,1:value and normal,2:normal component)
		T value, weight_,divergence_;
		TV position_,normal_;
		Array1D<TV> vector_gradient_;
		T epsilon;

	public:
		MLS_CONSTRAINT(const TV& position_input, const T& value_input)
		{
			position_ = position_input;
			value = value_input;
			type_ = 0;
			weight_ = (T)1;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const TV& position_input, const T& value_input, const int& type_input)
		{
			position_ = position_input;
			value = value_input;
			type_=type_input;
			weight_ = (T)1;
			epsilon = -1;
		}
		
		MLS_CONSTRAINT(const TV& position_input, const T& value_input, const TV& normal_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = 1;
			weight_ = (T)1;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = (T)1;
			epsilon= -1;
		}

		MLS_CONSTRAINT(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input, const T& weight_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = weight_input;
			epsilon = -1;
		}

		MLS_CONSTRAINT(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input, const T& weight_input, const T& divergence_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = weight_input;
			divergence_ = divergence_input;
			epsilon =- 1;
		}

		MLS_CONSTRAINT(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input, const T& weight_input, const T& divergence_input, const Array1D<TV>& vector_gradient_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = weight_input;
			divergence_ = divergence_input;
			vector_gradient_ = vector_gradient_input;
			epsilon = -1;
		}

		void Set(const TV& position_input, const T& value_input)
		{
			position_ = position_input;
			value = value_input;
			type_ = 0;
			weight_ = (T)1;
		}
		
		void Set(const TV& position_input, const T& value_input, const int& type_input)
		{
			position_ = position_input;
			value = value_input;
			type_ = type_input;
			weight_ = (T)1;
		}
		
		void Set(const TV& position_input, const T& value_input, const TV& normal_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = 1;
			weight_ = (T)1;
		}

		void Set(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = (T)1;
		}

		void Set(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input, const T& weight_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = weight_input;
		}

		void Set(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input, const T& weight_input, const T& divergence_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = weight_input;
			divergence_ = divergence_input;
		}

		void Set(const TV& position_input, const T& value_input, const TV& normal_input, const int& type_input, const T& weight_input, const T& divergence_input, const Array1D<TV>& vector_gradient_input)
		{
			position_ = position_input;
			value = value_input;
			normal_ = normal_input;
			type_ = type_input;
			weight_ = weight_input;
			divergence_ = divergence_input;
			vector_gradient_ = vector_gradient_input;
		}
	
	};//End of class MLS_CONSTRAINT

private:
	T epsilon, epsilon_squared, one_over_epsilon_squared;
	int size_of_basis, degree_of_basis;

public:
	std::vector<MLS_CONSTRAINT*> constraints;

public:
	MOVING_LEAST_SQUARES(const int& degree_of_basis_input = 0, const T& epsilon_input = (T)0.001, const int& number_of_constraints_input = 0)
	{
		delete_constraints_when_destroyed_ = true;

		SetDegreeOfBasis(degree_of_basis_input);
		SetEpsilon(epsilon_input);

		constraints.reserve(number_of_constraints_input);
		for(int i=0;i<number_of_constraints_input;i++) constraints.push_back(new MLS_CONSTRAINT(TV(), (T)0));
	}

	~MOVING_LEAST_SQUARES()
	{
		if(delete_constraints_when_destroyed_ == true)
		{
			const int number_of_constraints = (int)constraints.size();
			for(int i = 0; i < number_of_constraints; i++) 
				delete constraints[i];
		}
	}

	inline void Reset()
	{
		int number_of_constraints = (int)constraints.size();
		for(int i = 0; i < number_of_constraints; i++) 
			delete constraints[i];

		constraints.clear();
	}

public:
	inline void AddConstraint(const TV& position,const T& value)
	{
		constraints.push_back(new MLS_CONSTRAINT(position,value));
	}

	inline void AddConstraint(const TV& position,const T& value,const int& type)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, type));
	}

	inline void AddConstraint(const TV& position, const T& value, const TV& normal)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal));
	}

	inline void AddConstraint(const TV& position, const T& value, const TV& normal_input, const int& type)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal_input, type));
	}

	inline void AddConstraint(const TV& position, const T& value, const TV& normal_input, const int& type,const T& weight)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value, normal_input, type, weight));
	}

	inline void AddConstraint(const TV& position, const T& value,const TV& normal_input, const int& type,const T& weight,const T& divergence_input)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value,normal_input, type, weight, divergence_input));
	}

	inline void AddConstraint(const TV& position, const T& value, const TV& normal_input, const int& type, const T& weight, const T& divergence_input, const Array1D<TV>& vector_gradient_input)
	{
		constraints.push_back(new MLS_CONSTRAINT(position, value,normal_input, type, weight, divergence_input, vector_gradient_input));
	}

	inline void AddConstraint(MLS_CONSTRAINT *mls_constraint_)
	{
		constraints.push_back(mls_constraint_);
	}

	void SetDegreeOfBasis(const int& degree_of_basis_input);
	void SetEpsilon(const T& epsilon_input)
	{
		epsilon = epsilon_input;
		epsilon_squared = epsilon * epsilon;
		one_over_epsilon_squared = (T)1 / epsilon_squared;
	}

	VectorND<T> GetBasis(const TV& position);
	VectorND<D> GetBasis(const DV& position);
	VectorND<T> GetBasisDerivative(const TV& position, const int& d);// dimension should start from 0
	VectorND<D> GetBasisDerivative(const DV& position, const int& d);// dimension should start from 0
	VectorND<T> GetBasisSecondDerivative(const TV& position, const int& d);// dimension should start from 0

	T GetScalar(const TV& x);// scalar interpolation
	//T GetScalar_Old(const TV& x);// scalar interpolation
	//TV Get_Gradient(const TV& x);// gradient of the interpolated scalar	
	//TV Get_Diffusive_Second_Derivatives(const TV& x);
	//TV Get_Polygonal_Second_Derivatives(const TV& x);
	//TV Get_Numerical_Second_Derivatives(const TV& x,const T& dx);
	//TV Get_Diffuse_Gradient(const TV& x);
	//TV Get_Polynomial_Gradient(const TV& x);
	//TV Get_Numerical_Gradient(const TV& x,const T& dx);
	//T Get_Laplacian(const TV& x);
	//T Get_Diffuse_Laplacian(const TV& x);// laplacian of the interpolated scalar
	//T Get_Polynomial_Laplacian(const TV& x);
	//T Get_Numerical_Laplacian(const TV& x,const T& dx);
	//void Get_Poisson_Coefficients(const TV& x,VECTOR_ND<T>& coefficients);

	TV GetVector(const TV& x);// simple vector interpolation
	//TV Get_Vector_From_Normal_Components(const TV& x);// normal component constraints
	//TV Get_Vector_PDF1(const TV& x);// div-free condition at all ctr points
	//TV Get_Vector_PDF2(const TV& x);// div-free condition at all ctr points and x
	//TV Get_Vector_DIV(const TV& x,const T& divergence=(T)0);// only one div-free constraint at x (moving div-free constraint)
	//TV Get_Vector_DIV_W3(const TV& x,const T& divergence=(T)0);// only one div-free constraint at x with higher-order w
	//TV Get_Vector_DIV_W4(const TV& x,const T& divergence=(T)0);// only one div-free constraint at x with higher-order w
	//TV Get_Vector_DIV_And_Normal_Components(const TV& x);// normal component constraints and divergence free condition at x
	//TV Get_Vector_Numerical_DIV(const TV& x,const T& dx);// div-free at x by numerical differentiation
	//T Get_Divergence(const TV& x);

	//void Get_CMIP1(const TV& x,T& value,TV& derivatives);//Cubic scalar interpolation
	//void Get_CMIP2(const TV& x,T& value,TV& derivatives);//Cubic scalar interpolation
	//void Get_CMIP(const TV& x,TV& value,ARRAY<TV>& derivatives,
	//			  const bool use_value_constraints,const bool use_pdf_constraint,const bool use_indirect_derivatives,
	//			  const T value_constraint_weight,const T derivatives_weight,const T pdf_constraint_weight);//Cubic vector interpolation with PDF constraint

	//T Get_Numerical_Divergence(const TV& x,const T& dx);
	//TV Get_Numerical_Jacobian(const TV& x,const TV& velocity,const T& dx=1e-6);

	//void Print_Constraints();
	//void Report_Scalar_Error();//Errors at constraint positions
	//void Report_Vector_Error();//Errors at constraint positions

	////SPH implementation. This is temporary. Define SPH class later
	//T Get_SPH_Scalar(const TV& x);

};//End of class MLS

