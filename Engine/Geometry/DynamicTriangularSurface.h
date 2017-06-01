#pragma once

#include <iostream>
#include <list>
#include <vector>
#include <glm/glm.hpp>

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"
#include "Geometry/BOX_3D.h"
#include "Math/Matrix3X3.h"
#include "Math/Quaternion.h"

// Triangular geometry components
class DynamicVertex;
class DynamicEdge;
class DynamicTriangle;
class DynamicHole;

class DynamicTriangularSurface  
{
public:
	std::vector<DynamicVertex*>	vertices_;
	std::list<DynamicEdge*>		edges_;
	std::list<DynamicTriangle*>	triangles_;
	std::list<DynamicHole*>		holes_;

	BOX_3D<T> bounding_box_;
//	BOX bounding_box_;

	Array1D<DynamicTriangle*> triangle_array_;

public:
	DynamicTriangularSurface();	
	virtual ~DynamicTriangularSurface();

public:
	DynamicEdge* 		AddEdge(DynamicVertex* v0,DynamicVertex* v1);
	DynamicTriangle*	AddTriangle(int v0,int v1,int v2);
	DynamicTriangle*	AddTriangle(DynamicVertex* v0,DynamicVertex* v1,DynamicVertex* v2);
	DynamicTriangle*	AddTriangleWithConnectivity(DynamicVertex* v0,DynamicVertex* v1,DynamicVertex* v2);
	DynamicVertex* 	AddVertex(T *x);
	DynamicVertex* 	AddVertex(T *x, T *n);
	DynamicVertex* 	AddVertex(DynamicVertex* vertex);
	void		DelAllTriangles();
	void		DelTriangle(DynamicTriangle* );
	void		DetCurvatureNormals();
	void		DetVertexNormals();
	void		DetFaceNormals();
	void		AverageDuplexPositionNormal(DynamicTriangularSurface* neighbor, float dx);
	void		DetermineNormalDeviation();
	void		DetTextureCoordinates(T* xy,T* uv,T s);
	void		DetTextureCoordinates(T x,T y,T u,T v,T s);
	void		ChkBoundary();
	void		ChkBoundaryVertices(T* size);
	void		ChkTrianglesNeighborConnectivity();
	DynamicVertex* 	ChkNearestTextureCoordinate(T u, T v);
// 	void		DrawEdges();
// 	void		DrawFaceNormals();	
// 	void		DrawHoles();
// 	void		DrawHoles(int index);
// 	void		DrawTriangles(const bool& draw_front = true);
// 	void		DrawVertices();
// 	void		DrawVertexNormals();
// 	void		DrawVertexDeviation();
// 	void		DrawTrianglesNeighborConnectivity();
// 	void		DrawTrianglesCenter();
// 	void		DrawVerticesNeighborConnectivity();
// 	void		DrawCurvatureNormal(const T& scalar = 1);
// 	void		DrawVertexVelocity();
// 	void		DrawTrianglesColor(const bool& draw_front);
	void		Filtering(T lamda);
	void		Filtering2(T lamda);
	void		FilteringTangent(T lamda);
	void		FilteringMinimumVariation(T lamda);
	void		FilteringDeviation(T lamda);
	void		Collapse(DynamicTriangle* triangle);
	void		Collapse(DynamicVertex* v0,DynamicVertex* v1);
	void		FlipEdges();
	void		CorrectCCW();
	void		RemoveSmallTriangles(T threshold);
	void		RemoveSmallEdges(T threshold);
	void		RemoveLongTriangles(T ratio);
	void		RemoveVertexDeviation();
	void		UpdateBoundingBox();

	void		ButterflySubdivision();
	void		SimpleSubdivision();
	void		SplitLongEdges(const T& max_l);
	void		CollapseShortEdges(const T& min_l);
	void		Scale(const T& scale);
	void		Scale(const TV3& scale);
	void		Rotate(const QUATERNION& q);
	void		Rotate(const glm::mat4& r);
	void		Translate(T x,T y,T z);
	void		Translate(const TV3& trans);
	void		FillHoles();
	void		RemoveNonManifold();
	void		Reset();

/*
	void		TransS2W(const int thread_id, const DOMAIN_UNIFORM_3D* domain);
	void		TransW2S(const int thread_id, const DOMAIN_UNIFORM_3D* domain);

	void		TransS2WThreaded(const DOMAIN_UNIFORM_3D* domain);
	void		TransW2SThreaded(const DOMAIN_UNIFORM_3D* domain);
*/

	void		MinMaxPosition(TV3& min, TV3& max);

	// File I/O functions
	void		Read(const char *filename);
	void		ReadOBJ(const char *filename, const T angle = (T)0, const TV3 axis = TV3(1,1,1), const TV3 scale = TV3(1,1,1), const TV3 translation = TV3());
	void		ReadSMF(const char *filename);
	void		Write(const char *filename);
	void		WriteOBJ(const char *filename);
	void		WriteOBJ(const char *filename, TV3& position, QUATERNION& quat);
	void		WritePOLY(const char* filename);
	static void	WriteOBJ(Array1D<DynamicTriangularSurface*> &arr, const char *filename);
	// realflow bin mesh file exporter
	static void	WriteBinMeshFile(Array1D<DynamicTriangularSurface*> &arr, int i_frame, const char *filename);// realflow BIN mesh file을 export
	
	void		GetGeometryInfo(std::vector<TV3> &vertices,std::vector<TV3_INT> &faces);
	void		ArrangeTriangles();
	void		DertermineTriangleNormal();

	TV3			Center();
	void        InertiaTensor(MATRIX_3X3& inertia_mat, const T& mass);

	void		SamplingPointsOnSurface(std::vector<TV3>& arr, const T dw);
	void		SamplingPointsOnSurface(std::vector<TV3>& vert_arr, std::vector<TV3>& nor_arr, const T dw);
};
