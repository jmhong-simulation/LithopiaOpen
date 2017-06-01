// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "MovingLeastSquares.h"

void MOVING_LEAST_SQUARES::SetDegreeOfBasis(const int& degree_of_basis_input)
{
	degree_of_basis = degree_of_basis_input;
//	size_of_basis = Factorial(dimension+degree)/(Factorial(dimension)*Factorial(degree));// analytic calculation

	if(degree_of_basis == 0) size_of_basis = 1;
	else if(degree_of_basis == 1) 
		size_of_basis = 1 + dimension_;
	else if(degree_of_basis == 2)
	{
		if(dimension_ == 2) size_of_basis = 6;
		else if(dimension_==3) size_of_basis = 10;
	}
	else if(degree_of_basis == 3)
	{
		if(dimension_ == 2) size_of_basis = 10;
		else if(dimension_ == 3) size_of_basis = 20;
	}
	else
	{
		std::cout<<"MLS degree of basis is not defined with "<< degree_of_basis<<" degree"<<std::endl;
		exit(1); 
	}
}

VectorND<T> MOVING_LEAST_SQUARES::GetBasis(const TV& position)
{
	VectorND<T> basis;
	basis.Initialize(size_of_basis, true);//size_of_basis=1,3 in 2D and 1,4 in 3D
	if(degree_of_basis == 0) basis[0] = (T)1;
	else if(degree_of_basis == 1)
	{
		basis[0] = (T)1;
		for(int i = 1; i < size_of_basis; i++) basis[i] = position.values_[i-1];
	}
	else if(degree_of_basis == 2)
	{
		if(dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = x*x;
			basis[4] = y*y;
			basis[5] = x*y;
		}
		else if(dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = z;
			basis[4] = x*x;
			basis[5] = y*y;
			basis[6] = z*z;
			basis[7] = x*y;
			basis[8] = y*z;
			basis[9] = x*z;
		}
	}
	else if(degree_of_basis == 3)
	{
		if(dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = x*x;
			basis[4] = y*y;
			basis[5] = x*y;
			basis[6] = x*y*y;
			basis[7] = x*x*y;
			basis[8] = x*x*x;
			basis[9] = y*y*y;
		}
		else if(dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = z;
			basis[4] = x*x;
			basis[5] = y*y;
			basis[6] = z*z;
			basis[7] = x*y;
			basis[8] = y*z;
			basis[9] = x*z;
			basis[10] = x*x*y;
			basis[11] = x*x*z;
			basis[12] = x*y*y;
			basis[13] = y*y*z;
			basis[14] = x*z*z;
			basis[15] = y*z*z;
			basis[16] = x*y*z;
			basis[17] = x*x*x;
			basis[18] = y*y*y;
			basis[19] = z*z*z;
		}
	}

	return basis;
}

VectorND<D> MOVING_LEAST_SQUARES::GetBasis(const DV& position)
{
	VectorND<D> basis;
	basis.Initialize(size_of_basis, true);//size_of_basis=1,3 in 2D and 1,4 in 3D
	if (degree_of_basis == 0) basis[0] = (T)1;
	else if (degree_of_basis == 1)
	{
		basis[0] = (T)1;
		for (int i = 1; i < size_of_basis; i++) basis[i] = position.values_[i - 1];
	}
	else if (degree_of_basis == 2)
	{
		if (dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = x*x;
			basis[4] = y*y;
			basis[5] = x*y;
		}
		else if (dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = z;
			basis[4] = x*x;
			basis[5] = y*y;
			basis[6] = z*z;
			basis[7] = x*y;
			basis[8] = y*z;
			basis[9] = x*z;
		}
	}
	else if (degree_of_basis == 3)
	{
		if (dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = x*x;
			basis[4] = y*y;
			basis[5] = x*y;
			basis[6] = x*y*y;
			basis[7] = x*x*y;
			basis[8] = x*x*x;
			basis[9] = y*y*y;
		}
		else if (dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			basis[0] = 1;
			basis[1] = x;
			basis[2] = y;
			basis[3] = z;
			basis[4] = x*x;
			basis[5] = y*y;
			basis[6] = z*z;
			basis[7] = x*y;
			basis[8] = y*z;
			basis[9] = x*z;
			basis[10] = x*x*y;
			basis[11] = x*x*z;
			basis[12] = x*y*y;
			basis[13] = y*y*z;
			basis[14] = x*z*z;
			basis[15] = y*z*z;
			basis[16] = x*y*z;
			basis[17] = x*x*x;
			basis[18] = y*y*y;
			basis[19] = z*z*z;
		}
	}

	return basis;
}

T MOVING_LEAST_SQUARES::GetScalar(const TV& x)
{
	const int number_of_constraints = (int)constraints.size();

	int number_of_rows = 0;
	for(int i = 0; i < number_of_constraints; i ++)
	{
		switch(constraints[i]->type_)
		{
		case MLS_VALUE:
			number_of_rows ++;
			break;
		case MLS_VALUE_AND_NORMAL:
			number_of_rows ++;
			break;
		case MLS_DIFFUSIVE_NORMAL:
			number_of_rows += dimension_;
			break;
		}
	}

	MATRIX_MXN<D> B(number_of_rows, size_of_basis);
	VectorND<D> W(number_of_rows), WWphi(number_of_rows);

	int row = 0;
	for(int i = 0; i < number_of_constraints; i++)
	{
		MLS_CONSTRAINT *ctr = constraints[i];
//		assert(ctr->epsilon > 0);
		ctr->epsilon = epsilon;

		Vector3D<T> dev_t = ctr->position_ - x;
		Vector3D<D> deviation((D)dev_t.x_, (D)dev_t.y_, (D)dev_t.z_);
		D w = (D)1/(deviation.getSqrMagnitude() + ctr->epsilon*ctr->epsilon)*ctr->weight_;
		
		if(ctr->type_ == MLS_VALUE)
		{
			W.values_[row] = w;
			B.Set_Row(row, GetBasis(deviation));

			WWphi.values_[row] = w*w*ctr->value;

			row ++;
		}
		else if(ctr->type_ == MLS_VALUE_AND_NORMAL)
		{
			W.values_[row] = w;
			B.Set_Row(row, GetBasis(deviation));

//			std::cout<< GetBasis(deviation).values_[0] << " " << GetBasis(deviation).values_[1] << " " << GetBasis(deviation).values_[2] << " " << std::endl;

			WWphi.values_[row] = w*w*(ctr->value - dotProduct(dev_t, ctr->normal_));
			
			row ++;
		}
		else if(ctr->type_ == MLS_DIFFUSIVE_NORMAL)
		{
			for(int d = 0; d < dimension_; d++)
			{
				W.values_[row] = w;
				VectorND<D> value = GetBasisDerivative(deviation, d);

				B.Set_Row(row, GetBasisDerivative(deviation, d));

				WWphi.values_[row] = w*w*ctr->normal_.values_[d];
				row ++;
			}
		}
	}

// 	std::cout << B << std::endl;
// 	std::cout << " -----------"<<std::endl;
// 	std::cout << (B.Weighted_Normal_Equations_Matrix(W)).Cholesky_Solve(B.Transposed_Multiply(WWphi)) << std::endl;
// 
// 	exit(1);

	VectorND<D> c = (B.Weighted_Normal_Equations_Matrix(W)).Cholesky_Solve(B.Transposed_Multiply(WWphi));

//	std::cout<< c.values_[0] <<std::endl;

	return c.values_[0];
}

TV MOVING_LEAST_SQUARES::GetVector(const TV& x)// simple vector interpolation
{
	const int number_of_constraints = (int)constraints.size();

	if (number_of_constraints == 0)return TV();

	MATRIX_MXN<T> B(number_of_constraints*dimension_, size_of_basis*dimension_), WB(number_of_constraints*dimension_, size_of_basis*dimension_);
	VectorND<T> WWphi(number_of_constraints*(dimension_ + 1));

	int i = 0;
	for (int k = 0; k < number_of_constraints; k++)
	{
		MLS_CONSTRAINT *ctr = constraints[k];

		TV deviation = ctr->position_ - x;

		T w = (T)1 / (deviation.getSqrMagnitude() + epsilon_squared);
		VectorND<T> basis = GetBasis(deviation);
		
		for (int d = 0; d < dimension_; d++)
		{
			i++;

			for (int k = 0; k < size_of_basis; k++)
			{
				B(i, k*dimension_ + d) = basis[k];
				WB(i, k*dimension_ + d) = w*B(i, k*dimension_ + d);
			}

			WWphi[i] = w*w*ctr->value*ctr->normal_[d];
		}
	}

	VectorND<T> c = ((WB).Normal_Equations_Matrix()).Cholesky_Solve(B.Transpose()*WWphi);

	TV answer; for (int d = 0; d < dimension_; d++) answer.values_[d] = c[d];
	return answer;
}

VectorND<T> MOVING_LEAST_SQUARES::GetBasisDerivative(const TV& position, const int& d)// dimension should start from 0
{
	VectorND<T> basis_derivative;
	basis_derivative.Initialize(size_of_basis, true);

	if(degree_of_basis == 0)
		basis_derivative[0] = (T)0;
	else if(degree_of_basis == 1)
	{
		basis_derivative[0] = (T)0;
		if(size_of_basis > 1)
			basis_derivative[d+1] = (T)1;
	}
	else if(degree_of_basis == 2)
	{
		if(dimension_ == 2)
		{
			T x = position.values_[0],y = position.values_[1];

			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 2 * x;
				basis_derivative[4] = 0;
				basis_derivative[5] = y;
			}
			else if(d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * y;
				basis_derivative[5] = x;
			}
		}
		else if(dimension_==3)
		{
			T x=position.values_[0], y = position.values_[1], z = position.values_[2];
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * x;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = y;
				basis_derivative[8] = 0;
				basis_derivative[9] = z;
			}
			else if(d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 2 * y;
				basis_derivative[6] = 0;
				basis_derivative[7] = x;
				basis_derivative[8] = z;
				basis_derivative[9] = 0;
			}
			else if(d == 2)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 1;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2 * z;
				basis_derivative[7] = 0;
				basis_derivative[8] = y;
				basis_derivative[9] = x;
			}
		}
	}
	else if(degree_of_basis == 3)
	{
		if(dimension_ == 2)
		{
			T x=position.values_[0], y = position.values_[1];
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 2 * x;
				basis_derivative[4] = 0;
				basis_derivative[5] = y;
				basis_derivative[6] = y*y;
				basis_derivative[7] = 2*x*y;
				basis_derivative[8] = 3*x*x;
				basis_derivative[9] = 0;
			}
			else if(d==1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2*y;
				basis_derivative[5] = x;
				basis_derivative[6] = 2*x*y;
				basis_derivative[7] = x*x;
				basis_derivative[8] = 0;
				basis_derivative[9] = 3 * y * y;
			}
		}
		else if(dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * x;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = y;
				basis_derivative[8] = 0;
				basis_derivative[9] = x;
				basis_derivative[10] = 2 * x * y;
				basis_derivative[11] = 2 * x * z;
				basis_derivative[12] = y * y;
				basis_derivative[13] = 0;
				basis_derivative[14] = z * z;
				basis_derivative[15] = 0;
				basis_derivative[16] = y * z;
				basis_derivative[17]= 3 * x * x;
				basis_derivative[18]= 0;
				basis_derivative[19]= 0;
			}
			else if(d==1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 2 * y;
				basis_derivative[6] = 0;
				basis_derivative[7] = x;
				basis_derivative[8] = y;
				basis_derivative[9] = 0;
				basis_derivative[10] = x * x;
				basis_derivative[11] = 0;
				basis_derivative[12] = 2 * x * y;
				basis_derivative[13] = 2 * y * z;
				basis_derivative[14] = 0;
				basis_derivative[15] = z * z;
				basis_derivative[16] = x * z;
				basis_derivative[17] = 0;
				basis_derivative[18] = 3 * y * y;
				basis_derivative[19] = 0;
			}
			else if(d==2)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 1;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2 * z;
				basis_derivative[7] = 0;
				basis_derivative[8] = y;
				basis_derivative[9] = x;
				basis_derivative[10] = 0;
				basis_derivative[11] = x * x;
				basis_derivative[12] = 0;
				basis_derivative[13] = y * y;
				basis_derivative[14] = 2 * x * z;
				basis_derivative[15] = 2 * y * z;
				basis_derivative[16] = x * y;
				basis_derivative[17] = 0;
				basis_derivative[18] = 0;
				basis_derivative[19] = 3 * z * z;
			}
		}
	}
	return basis_derivative;
}

VectorND<D> MOVING_LEAST_SQUARES::GetBasisDerivative(const DV& position, const int& d)// dimension should start from 0
{
	VectorND<D> basis_derivative;
	basis_derivative.Initialize(size_of_basis, true);

	if (degree_of_basis == 0)
		basis_derivative[0] = (T)0;
	else if (degree_of_basis == 1)
	{
		basis_derivative[0] = (T)0;
		if (size_of_basis > 1)
			basis_derivative[d + 1] = (T)1;
	}
	else if (degree_of_basis == 2)
	{
		if (dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];

			if (d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 2 * x;
				basis_derivative[4] = 0;
				basis_derivative[5] = y;
			}
			else if (d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * y;
				basis_derivative[5] = x;
			}
		}
		else if (dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			if (d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * x;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = y;
				basis_derivative[8] = 0;
				basis_derivative[9] = z;
			}
			else if (d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 2 * y;
				basis_derivative[6] = 0;
				basis_derivative[7] = x;
				basis_derivative[8] = z;
				basis_derivative[9] = 0;
			}
			else if (d == 2)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 1;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2 * z;
				basis_derivative[7] = 0;
				basis_derivative[8] = y;
				basis_derivative[9] = x;
			}
		}
	}
	else if (degree_of_basis == 3)
	{
		if (dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];
			if (d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 2 * x;
				basis_derivative[4] = 0;
				basis_derivative[5] = y;
				basis_derivative[6] = y*y;
				basis_derivative[7] = 2 * x*y;
				basis_derivative[8] = 3 * x*x;
				basis_derivative[9] = 0;
			}
			else if (d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * y;
				basis_derivative[5] = x;
				basis_derivative[6] = 2 * x*y;
				basis_derivative[7] = x*x;
				basis_derivative[8] = 0;
				basis_derivative[9] = 3 * y * y;
			}
		}
		else if (dimension_ == 3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			if (d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 1;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2 * x;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = y;
				basis_derivative[8] = 0;
				basis_derivative[9] = x;
				basis_derivative[10] = 2 * x * y;
				basis_derivative[11] = 2 * x * z;
				basis_derivative[12] = y * y;
				basis_derivative[13] = 0;
				basis_derivative[14] = z * z;
				basis_derivative[15] = 0;
				basis_derivative[16] = y * z;
				basis_derivative[17] = 3 * x * x;
				basis_derivative[18] = 0;
				basis_derivative[19] = 0;
			}
			else if (d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 1;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 2 * y;
				basis_derivative[6] = 0;
				basis_derivative[7] = x;
				basis_derivative[8] = y;
				basis_derivative[9] = 0;
				basis_derivative[10] = x * x;
				basis_derivative[11] = 0;
				basis_derivative[12] = 2 * x * y;
				basis_derivative[13] = 2 * y * z;
				basis_derivative[14] = 0;
				basis_derivative[15] = z * z;
				basis_derivative[16] = x * z;
				basis_derivative[17] = 0;
				basis_derivative[18] = 3 * y * y;
				basis_derivative[19] = 0;
			}
			else if (d == 2)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 1;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2 * z;
				basis_derivative[7] = 0;
				basis_derivative[8] = y;
				basis_derivative[9] = x;
				basis_derivative[10] = 0;
				basis_derivative[11] = x * x;
				basis_derivative[12] = 0;
				basis_derivative[13] = y * y;
				basis_derivative[14] = 2 * x * z;
				basis_derivative[15] = 2 * y * z;
				basis_derivative[16] = x * y;
				basis_derivative[17] = 0;
				basis_derivative[18] = 0;
				basis_derivative[19] = 3 * z * z;
			}
		}
	}
	return basis_derivative;
}

VectorND<T> MOVING_LEAST_SQUARES::GetBasisSecondDerivative(const TV& position, const int& d)// dimension should start from 0
{
	VectorND<T> basis_derivative(size_of_basis);
	if(degree_of_basis == 0)
		basis_derivative[0] = (T)0;
	else if(degree_of_basis == 1)
	{
		basis_derivative[0] = (T)0;
		if(size_of_basis>1)
			basis_derivative[d+1] = (T)0;
	}
	else if(degree_of_basis == 2)
	{
		if(dimension_ == 2)
		{
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 2;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
			}
			else if(d == 1)
			{					
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2;
				basis_derivative[5] = 0;
			}
		}
		else if(dimension_ == 3)
		{
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = 0;
				basis_derivative[8] = 0;
				basis_derivative[9] = 0;
			}
			else if(d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 2;
				basis_derivative[6] = 0;
				basis_derivative[7] = 0;
				basis_derivative[8] = 0;
				basis_derivative[9] = 0;
			}
			else if(d==2)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2;
				basis_derivative[7] = 0;
				basis_derivative[8] = 0;
				basis_derivative[9] = 0;
			}
		}
	}
	else if(degree_of_basis == 3)
	{
		if(dimension_ == 2)
		{
			T x = position.values_[0], y = position.values_[1];
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 2;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = 2 * y;
				basis_derivative[8] = 6 * x;
				basis_derivative[9] = 0;
			}
			else if(d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2 * x;
				basis_derivative[7] = 0;
				basis_derivative[8] = 0;
				basis_derivative[9] = 6 * y;
			}
		}
		else if(dimension_==3)
		{
			T x = position.values_[0], y = position.values_[1], z = position.values_[2];
			if(d == 0)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 2;
				basis_derivative[5] = 0;
				basis_derivative[6] = 0;
				basis_derivative[7] = 0;
				basis_derivative[8] = 0;
				basis_derivative[9] = 1;
				basis_derivative[10] = 2 * y;
				basis_derivative[11] = 2 * z;
				basis_derivative[12] = 0;
				basis_derivative[13] = 0;
				basis_derivative[14] =0;
				basis_derivative[15] = 0;
				basis_derivative[16] = 0;
				basis_derivative[17] = 6 * x;
				basis_derivative[18] = 0;
				basis_derivative[19] = 0;
			}
			else if(d == 1)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0; 
				basis_derivative[4] = 0;
				basis_derivative[5] = 2;
				basis_derivative[6] = 0;
				basis_derivative[7] = 0;
				basis_derivative[8] = 1;
				basis_derivative[9] = 0;
				basis_derivative[10] = 0;
				basis_derivative[11] = 0;
				basis_derivative[12] = 2 * x;
				basis_derivative[13] = 2 * z;
				basis_derivative[14] = 0;
				basis_derivative[15] = 0;
				basis_derivative[16] = 0;
				basis_derivative[17] = 0;
				basis_derivative[18] = 6 * y;
				basis_derivative[19] = 0;
			}
			else if(d == 2)
			{
				basis_derivative[0] = 0;
				basis_derivative[1] = 0;
				basis_derivative[2] = 0;
				basis_derivative[3] = 0;
				basis_derivative[4] = 0;
				basis_derivative[5] = 0;
				basis_derivative[6] = 2;
				basis_derivative[7] = 0;
				basis_derivative[8] = 0;
				basis_derivative[9] = 0;
				basis_derivative[10] = 0;
				basis_derivative[11] = 0;
				basis_derivative[12] = 0;
				basis_derivative[13] = 0;
				basis_derivative[14] = 2 * x;
				basis_derivative[15] = 2 * y;
				basis_derivative[16] = 0;
				basis_derivative[17] = 0; 
				basis_derivative[18] = 0;
				basis_derivative[19] = 6 * z;
			}
		}
	}
	return basis_derivative;
}
