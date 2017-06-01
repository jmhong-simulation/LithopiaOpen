#pragma once

#include "GENERIC_DEFINITIONS.h"

#include <list>
#include <vector>
#include <iostream>

class DynamicTriangle;
class DynamicEdge;

class DynamicVertex
{
public:
	TV3	x_;
	TV3	n_;
	TV3	s_;// tension stress

public:
	TV3	velocity_;
	TV3 force_;
	T mass_;
	int grain_boundary_number_;
	T phi_;

public:
	T	reaction_speed_, reaction_speed_dot_;
	T	stress_;// det of tension stress
	TV3	deviation_;
	T	normal_deviation_;
	TV3	curvature_normal_;
	T	normalizer_;
	
	int		sampled_direction_;
	int		index_;

	bool	feature_;
	bool	is_mc_vertex_;
	bool	is_boundary_;

	std::list<DynamicEdge*> edges_;
	std::list<DynamicTriangle*> triangles_;
	
public:
	DynamicVertex();
	DynamicVertex(const T* x);
	DynamicVertex(TV3* x); 
	DynamicVertex(T* x,T* n);
	DynamicVertex(TV3* x,TV3* n);
	DynamicVertex(const T x, const T y, const T z, const T nx, const T ny, const T nz);

	~DynamicVertex();

public:
	void	SetPosition(T* x);
	void	SetNormal(T* n);
	T*		GetPosition();
	T*		GetNormal();
	void	Reset();
	void	AddTriangle(DynamicTriangle*  triangle);
	void	DelTriangle(DynamicTriangle*  triangle);
	void	AddEdge(DynamicEdge*  edge);
	void	DelEdge(DynamicEdge*  edge);
	void	DetermineNormal();
	void	DetermineAverageNormal();
	void	DrawNormal();
	void	DrawCurvatureNormal(const T& scalar = 1);
	void	DrawNeighborConnectivity();
	void	DrawDeviation();
	void	DrawVelocity();
	void	Replace(DynamicVertex* v_after);
	void	DetCurvatureNormal();		

	TV3 &GetTVPosition() { return x_;}
	TV3 &GetTVNormal()   { return n_;}

public:
	void PrintPos()
	{
		printf("%f %f %f \n", x_.x_, x_.y_, x_.z_);
	}
};