#pragma once

#include <list>
#include <vector>
#include <iostream>
#include "DynamicEdge.h"
#include "GENERIC_DEFINITIONS.h"

class DynamicVertex;

////////////////////////////////////////////////////////////////////
// TRIANGLE
////////////////////////////////////////////////////////////////////
class DynamicTriangle
{
public:
	TV3 n_;
	TV3 uv_[3];// texture coordinates
	T area_;
	DynamicVertex* vertices_[3];
	DynamicVertex* edge_vertex_[3];
	DynamicTriangle* triangles_[3];
	DynamicEdge *edges_[3];
	bool is_old_;
	bool wrong_;
	bool is_surface_;
	TV3_INT vertex_indices_;

public:
	DynamicTriangle();
	DynamicTriangle(DynamicVertex* v0,DynamicVertex* v1,DynamicVertex* v2);
	~DynamicTriangle();

public:
	void		DelTriangle(DynamicTriangle* triangle);
	void		DetermineNormal();
	void		DetermineAverageNormal();
	T*			GetNormal();
	void		SetNormal(T* n);
	void		DrawNormal();
	void		Draw();
	void		DrawBack();
	void		DrawEdges();
	void		DrawCenter();
	void		DrawNeighborConnectivity();
	void		CorrectCCW();
	int			CountFeatureVertex();
	void		FindIntersection(int direction,T* p,T* p_intersection,T* uv);		// ?? ???⿡?? ???�ٺ? p ??ǥ (??Į?? 2??), ?׸??? uv??ǥ (??Į?? 2??)	
	bool		IsInside(DynamicVertex* v);
	void		Flip();//?ڱ? ?ڽ??? ?޸??? �?Ŵ? Triangular Surface???? ?Ѵ?.
	void		Flip(DynamicVertex* v0,DynamicVertex* v1);
	DynamicVertex*	FindAnotherVertex(DynamicVertex* v0,DynamicVertex* v1);
	int			GetNeighborIndex(DynamicTriangle* triangle);
	int			GetVertexIndex(DynamicVertex* v);
	void		ChkNeighborConnectivity();
	T			GetOppositeEdgeLength(DynamicVertex* v);
	void		AddLocalCurvatureNormal();
	void		SetEdge(DynamicEdge *e0, DynamicEdge *e1, DynamicEdge *e2);

	TV3 &GetTVNormal();
	TV3 GetTVCenter();

	void GetCenter(T* center);
	T SignedDistance(const TV3& location) const;

	TV3 ClosestPoint(const TV3& location) const;

	TV3 ClosestPoint(const TV3& location, bool& on_triangle) const;

	static TV3 BarycentricCoordinates(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3);

	static TV3 ClosestPointFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3);

	static TV3 ClosestPointFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3, bool& on_triangle);

	static T SignedDistanceFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3);

	static T SignedDistanceFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3, bool& on_triangle);

	static TV3 ClosestPointFromLine(const TV3& location, const TV3& x1, const TV3& x2);

	static bool IntersectTriangle(const TV3& i0, const TV3& i1, const DynamicTriangle& triangle, TV3& intersection_position);

	static bool RayThruTriangle(const DynamicTriangle& triangle, const TV3& R1, const TV3& R2,  TV3& intersection_position_out, const T& offset);
};