#include <algorithm>
#include <direct.h>
#include <fstream>
#include <string.h>
#include <boost/format.hpp>

#include "io.h"
#include "math.h"
#include "DynamicEdge.h"
#include "DynamicHole.h"
#include "DynamicTriangle.h"
#include "DynamicTriangularSurface.h"
#include "DynamicVertex.h"
#include "DataStructure/Vector3D.h"

/*
#include "GL/glut.h"
#include "GL/gl.h"
*/

//////////////////////////////////////////////////////////////////////
// Triangular Surface Construction/Destruction
//////////////////////////////////////////////////////////////////////
DynamicTriangularSurface::DynamicTriangularSurface()
{
}

DynamicTriangularSurface::~DynamicTriangularSurface()
{
	Reset();
}

//////////////////////////////////////////////////////////////////////
// Triangular Surface Public methods --> Interfaces
//////////////////////////////////////////////////////////////////////
void DynamicTriangularSurface::Reset()
{
	TRAVERSE_VERTICES
	{
		SAFE_DELETE(*itr_vertex);
	}	

	TRAVERSE_EDGES
	{
		SAFE_DELETE(*itr_edge);
	}

	TRAVERSE_TRIANGLES
	{
		SAFE_DELETE(*itr_triangle);
	}

	vertices_.clear();
	edges_.clear();
	triangles_.clear();
}

DynamicVertex* DynamicTriangularSurface::AddVertex(T *x)
{
	DynamicVertex*  vertex = new DynamicVertex(x);
	DynamicTriangularSurface::vertices_.push_back(vertex);
	return vertex;
}

DynamicVertex* DynamicTriangularSurface::AddVertex(T *x, T *n)
{
	DynamicVertex*  vertex = new DynamicVertex(x, n);
	DynamicTriangularSurface::vertices_.push_back(vertex);
	return vertex;
}

DynamicVertex* DynamicTriangularSurface::AddVertex(DynamicVertex*  vertex)
{
	DynamicTriangularSurface::vertices_.push_back(vertex);
	return vertex;
}

DynamicEdge* DynamicTriangularSurface::AddEdge(DynamicVertex*  v0, DynamicVertex*  v1)
{
	DynamicEdge*  edge = new DynamicEdge(v0, v1);
	DynamicTriangularSurface::edges_.push_back(edge);
	return edge;
}

DynamicTriangle* DynamicTriangularSurface::AddTriangle(int v0,int v1,int v2)
{
	return DynamicTriangularSurface::AddTriangle(DynamicTriangularSurface::vertices_[v0], DynamicTriangularSurface::vertices_[v1], DynamicTriangularSurface::vertices_[v2]);
}

void DynamicTriangularSurface::Collapse(DynamicVertex*  v0, DynamicVertex*  v1)
{
	TV3 center;
	center = v0->x_; 
	center += v1->x_;
	center /= 2.0f;

	TV3 center_normal;
	center_normal = v0->n_;
	center_normal += v1->n_;
	center_normal /= 2.0f;
	
	TV3 deviation;
	deviation = v0->deviation_;
	deviation += v1->deviation_;
	deviation /= 2.0f;

	DynamicVertex*  v_after = DynamicTriangularSurface::AddVertex(center.values_, center_normal.values_);
	v_after->reaction_speed_ = (v0->reaction_speed_+v1->reaction_speed_)*(T)0.5;

	v_after->deviation_ = deviation;

	if(v0->is_mc_vertex_ == true || v1->is_mc_vertex_ == true)
	{
		v_after->is_mc_vertex_ = true;
	}
	// v0?? v1?? ?ﰢ?? ???? �???? v_after???? ???????ش?.
	std::list <DynamicTriangle* >::iterator itr_triangle;
	for(itr_triangle = v0->triangles_.begin(); itr_triangle != v0->triangles_.end(); itr_triangle ++)
	{
		v_after->AddTriangle(*itr_triangle);
	}
	v0->triangles_.erase(v0->triangles_.begin(), v0->triangles_.end());

	for(itr_triangle = v1->triangles_.begin(); itr_triangle != v1->triangles_.end(); itr_triangle ++)
	{
		v_after->AddTriangle(*itr_triangle);
	}
	v1->triangles_.erase(v1->triangles_.begin(), v1->triangles_.end());

	for(itr_triangle = v_after->triangles_.begin(); itr_triangle != v_after->triangles_.end(); itr_triangle ++)
	{	
		DynamicTriangle*  triangle = *itr_triangle;

		for(int d = 0; d < 3; d ++)
		{
			if((triangle->vertices_[d] == v0) || (triangle->vertices_[d] == v1))
			{
				triangle->vertices_[d] = v_after;
			}		
		}
	}

	for(itr_triangle = v_after->triangles_.begin(); itr_triangle != v_after->triangles_.end(); itr_triangle ++)
	{
		DynamicTriangle*  triangle = *itr_triangle;
		int num = 0;
		for(int d = 0; d < 3; d ++)
		{
			if(triangle->vertices_[d] == v_after)
				num ++;
		}

		if(num >= 2)
		{
			v_after->triangles_.remove(triangle);
			DynamicTriangularSurface::DelTriangle(triangle);
			itr_triangle = v_after->triangles_.begin();
		}
	}

	// triangle neighbor �?? ó??
	DynamicTriangularSurface::ChkTrianglesNeighborConnectivity();
}

void DynamicTriangularSurface::Collapse(DynamicTriangle*  triangle)
{
	DynamicVertex*  v[3];
	v[0] = triangle->vertices_[0];
	v[1] = triangle->vertices_[1];
	v[2] = triangle->vertices_[2];

	TV3 center; 
	center = v[0]->x_;
	center += v[1]->x_;
	center += v[2]->x_;
	center /= 3.0f;
	
	TV3 normal;
	normal = v[0]->n_;
	normal += v[1]->n_;
	normal += v[2]->n_;
	normal /= 3.0f;

	TV3 deviation;
	deviation = v[0]->deviation_;
	deviation += v[1]->deviation_;
	deviation += v[2]->deviation_;
	deviation /= 3.0f;

	DynamicVertex*  v_after = DynamicTriangularSurface::AddVertex(center.values_, normal.values_);
	v_after->deviation_ = deviation;
	
	if(v[0]->is_mc_vertex_ == true || v[1]->is_mc_vertex_ == true || v[2]->is_mc_vertex_ == true)
	{
		v_after->is_mc_vertex_ = true;
	}

	// ??ĥ vertex???? ?̿? triangle �???? ???????ش?.
	std::list <DynamicTriangle* >::iterator itr_triangle;
	for(int d = 0; d < 3; d ++)
	{
		for(itr_triangle = v[d]->triangles_.begin(); itr_triangle != v[d]->triangles_.end(); itr_triangle ++)
		{
			v_after->AddTriangle(*itr_triangle);			
		}

		v[d]->triangles_.erase(v[d]->triangles_.begin(), v[d]->triangles_.end());
	}

//	cout << v_after->triangles_.size() << endl;

	// ?????? vertex?鿡 ???ؼ? ?????? ?ִ? ?ﰢ?????? vertex�???? ?????Ͽ??ش?.
	for(itr_triangle = v_after->triangles_.begin(); itr_triangle != v_after->triangles_.end(); itr_triangle ++)
	{	
		DynamicTriangle*  triangle = *itr_triangle;		

		for(int d = 0; d < 3; d ++)
		{			
			if((triangle->vertices_[d] == v[0]) || (triangle->vertices_[d] == v[1]) || (triangle->vertices_[d] == v[2]))
			{
				triangle->vertices_[d] = v_after;
			}
		}
	}

	// v_after?? ?? ?? ?̻? ?????? ?ִ? ?ﰢ????� ???? ???δ?.
	for(itr_triangle = v_after->triangles_.begin(); itr_triangle != v_after->triangles_.end(); itr_triangle ++)
	{
		DynamicTriangle*  triangle = *itr_triangle;
		int num = 0;
		for(int d = 0; d <= 3; d ++)
		{
			if(triangle->vertices_[d] == v_after)
				num ++;
		}

		if(num >= 2)
		{
			v_after->triangles_.remove(triangle);
			DynamicTriangularSurface::DelTriangle(triangle);
			itr_triangle = v_after->triangles_.begin();
		}
	}

	// edge???? ???? ã???ֱ?
}

void DynamicTriangularSurface::DelTriangle(DynamicTriangle*  triangle)
{
	DynamicTriangularSurface::triangles_.remove(triangle);
	
	for(int i = 0; i < 3; i ++)
	{
		triangle->vertices_[i]->DelTriangle(triangle);
		if(triangle->triangles_[i] != NULL)
		{
			triangle->triangles_[i]->DelTriangle(triangle);
//			triangle->triangles_[i] = NULL;
		}
	}

	delete triangle;
}

void DynamicTriangularSurface::DelAllTriangles()
{
	std::list <DynamicTriangle* >::iterator itr_triangle = DynamicTriangularSurface::triangles_.begin();
	while(itr_triangle != DynamicTriangularSurface::triangles_.end())
	{
		DynamicTriangle*  triangle = *itr_triangle;
		for(int i = 0; i < 3; i ++)
		{
			triangle->vertices_[i]->DelTriangle(triangle);
			if(triangle->triangles_[i] != NULL)
			{
				triangle->triangles_[i]->DelTriangle(triangle);
				triangle->triangles_[i] = NULL;
			}
		}
		
		delete triangle;

		DynamicTriangularSurface::triangles_.erase(itr_triangle);

		itr_triangle = DynamicTriangularSurface::triangles_.begin();
	}
}

DynamicTriangle* DynamicTriangularSurface::AddTriangle(DynamicVertex* v0,DynamicVertex* v1,DynamicVertex* v2)
{
	DynamicVertex *v[3]={v0,v1,v2};
	DynamicTriangle *triangle=new DynamicTriangle(v[0],v[1],v[2]);
	triangles_.push_back(triangle);

	// update vertex-triangle connectivity
	v[0]->AddTriangle(triangle);
	v[1]->AddTriangle(triangle);
	v[2]->AddTriangle(triangle);

	return triangle;	
}

DynamicTriangle* DynamicTriangularSurface::AddTriangleWithConnectivity(DynamicVertex* v0,DynamicVertex* v1,DynamicVertex* v2)
{
	DynamicVertex *v[3]={v0,v1,v2};
	DynamicTriangle *triangle=new DynamicTriangle(v[0],v[1],v[2]);
	triangles_.push_back(triangle);

	// update triangle-triangle connectivity
	std::list<DynamicTriangle*>::iterator itr_triangle, itr_triangle2;
	DynamicTriangle *triangle2;
	for(int i = 0; i < 3; ++i)
	{
		triangle2 = NULL;
		int index = 0;
		for(itr_triangle = v[i]->triangles_.begin(); itr_triangle != v[i]->triangles_.end(); itr_triangle++)
		{
			for(itr_triangle2 = v[(i+1)%3]->triangles_.begin(); itr_triangle2 != v[(i+1)%3]->triangles_.end(); ++itr_triangle2)
			{
				if((*itr_triangle) == (*itr_triangle2))
				{
					triangle2 = (*itr_triangle);
					break;
				}
			}
			if(triangle2 != NULL)
			{
				if(triangle2->vertices_[0] == v[i])			index+=0;
				else if(triangle2->vertices_[1] == v[i])	index+=1;
				else if(triangle2->vertices_[2] == v[i])	index+=2;

				if(triangle2->vertices_[0] == v[(i+1)%3])		index+=0;
				else if(triangle2->vertices_[1] == v[(i+1)%3])	index+=1;
				else if(triangle2->vertices_[2] == v[(i+1)%3])	index+=2;

				if(index == 1)		triangle2->triangles_[0] = triangle;
				else if(index == 2) triangle2->triangles_[2] = triangle;
				else if(index == 3)	triangle2->triangles_[1] = triangle;

				triangle->triangles_[i]=triangle2;
				break;
			}
		}

		triangle->triangles_[i] = triangle2;
	}

	// update vertex-triangle connectivity
	v[0]->AddTriangle(triangle);
	v[1]->AddTriangle(triangle);
	v[2]->AddTriangle(triangle);

	return triangle;	
}

/*
void DynamicTriangularSurface::DrawVertices()
{	
	glBegin(GL_POINTS);
	for(int i=0;i<(int)vertices_.size();i++)
	{		
		glNormal3fv((float*)vertices_[i]->n_.values_);
		glVertex3fv((float*)vertices_[i]->x_.values_);
	}
	glEnd();
}
*/

void DynamicTriangularSurface::Scale(const T& scale)
{
	for(int i = 0; i < (int)vertices_.size(); i++)
	{
		for(int d = 0; d < 3; d++)
		{
			vertices_[i]->x_.values_[d] *= scale;
		}
	}
	bounding_box_.Scale(scale);
}

void DynamicTriangularSurface::Scale(const TV3& scale)
{
	for(int i = 0; i < (int)vertices_.size(); i++)
	{
		vertices_[i]->x_.x_ *= scale.x_;
		vertices_[i]->x_.y_ *= scale.y_;
		vertices_[i]->x_.z_ *= scale.z_;
	}
	bounding_box_.Scale(scale);
}

void DynamicTriangularSurface::Rotate(const glm::mat4& r)
{
	for(int i = 0; i < (int)vertices_.size(); i++)
	{
		glm::vec4 vt(vertices_[i]->x_.x_, vertices_[i]->x_.y_, vertices_[i]->x_.z_, 1.0f);
		vt = r*vt;

		vertices_[i]->x_.x_ = vt.x;
		vertices_[i]->x_.y_ = vt.y;
		vertices_[i]->x_.z_ = vt.z;
	}

	UpdateBoundingBox();
}

void DynamicTriangularSurface::Rotate(const QUATERNION& q)
{
	for(int i = 0; i < (int)vertices_.size(); i++)
	{
		TV3 vt(vertices_[i]->x_.x_, vertices_[i]->x_.y_, vertices_[i]->x_.z_);
		vt = q.Rotate(vt);

		vertices_[i]->x_.x_ = vt.x_;
		vertices_[i]->x_.y_ = vt.y_;
		vertices_[i]->x_.z_ = vt.z_;
	}

	UpdateBoundingBox();
}

void DynamicTriangularSurface::Translate(T x,T y,T z)
{
	for(int i = 0; i < (int)vertices_.size(); i++)
	{
		vertices_[i]->x_.values_[0] += x;
		vertices_[i]->x_.values_[1] += y;
		vertices_[i]->x_.values_[2] += z;
	}

	bounding_box_.Translate(TV3(x, y, z));
}

void DynamicTriangularSurface::Translate(const TV3& trans)
{
	for(int i = 0; i < (int)vertices_.size(); i++)
	{
		vertices_[i]->x_.values_[0] += trans.x_;
		vertices_[i]->x_.values_[1] += trans.y_;
		vertices_[i]->x_.values_[2] += trans.z_;
	}

	bounding_box_.Translate(trans);
}

/*
void TRIANGULAR_SURFACE::TransS2W(const int thread_id, const DOMAIN_UNIFORM_3D* domain)
{
	PREPARE_FOR_1D_ITERATION((int)vertices_.size());

	BEGIN_1D_ITERATION
	{
		vertices_[p]->x_ = domain->TransPosS2W(vertices_[p]->x_);
	}
	END_1D_ITERATION;

	BEGIN_HEAD_THREAD_WORK
	{
		TV3 box_min = domain->TransPosS2W(TV3(bounding_box_.x_min_,bounding_box_.y_min_,bounding_box_.z_min_));
		TV3 box_max = domain->TransPosS2W(TV3(bounding_box_.x_max_,bounding_box_.y_max_,bounding_box_.z_max_));

		bounding_box_.x_min_ = box_min.x_;
		bounding_box_.y_min_ = box_min.y_;
		bounding_box_.z_min_ = box_min.z_;

		bounding_box_.x_max_ = box_max.x_;
		bounding_box_.y_max_ = box_max.y_;
		bounding_box_.z_max_ = box_max.z_;	
	}
	END_HEAD_THREAD_WORK;
}
*/

/*
void TRIANGULAR_SURFACE::TransW2S(const int thread_id, const DOMAIN_UNIFORM_3D* domain)
{
	PREPARE_FOR_1D_ITERATION((int)vertices_.size());

	BEGIN_1D_ITERATION
	{
		vertices_[p]->x_ = domain->TransPosW2S(vertices_[p]->x_);
	}
	END_1D_ITERATION;

	BEGIN_HEAD_THREAD_WORK
	{
		TV3 box_min = domain->TransPosW2S(TV3(bounding_box_.x_min_,bounding_box_.y_min_,bounding_box_.z_min_));
		TV3 box_max = domain->TransPosW2S(TV3(bounding_box_.x_max_,bounding_box_.y_max_,bounding_box_.z_max_));

		bounding_box_.x_min_ = box_min.x_;
		bounding_box_.y_min_ = box_min.y_;
		bounding_box_.z_min_ = box_min.z_;

		bounding_box_.x_max_ = box_max.x_;
		bounding_box_.y_max_ = box_max.y_;
		bounding_box_.z_max_ = box_max.z_;		
	}
	END_HEAD_THREAD_WORK;
}
*/

/*
void TRIANGULAR_SURFACE::TransS2WThreaded(const DOMAIN_UNIFORM_3D* domain)
{
	GET_MULTITHREADING->RunThreads(&TRIANGULAR_SURFACE::TransS2W, this, domain);
}
*/

/*
void TRIANGULAR_SURFACE::TransW2SThreaded(const DOMAIN_UNIFORM_3D* domain)
{
	GET_MULTITHREADING->RunThreads(&TRIANGULAR_SURFACE::TransW2S, this, domain);
}
*/

void DynamicTriangularSurface::UpdateBoundingBox()
{
	if((int)vertices_.size() < 1)
	{
		return;
	}

	bounding_box_ = BOX_3D<T>(TV3(vertices_[0]->x_), TV3(vertices_[0]->x_));

	for(int i = 1; i < (int)vertices_.size(); i++)
	{
		bounding_box_.Include(TV3(vertices_[i]->x_));
	}
}

void DynamicTriangularSurface::ChkBoundaryVertices(T *size)
{
	TRAVERSE_VERTICES
	{
		for(int d=0;d<3;d++)
		{
			if((*itr_vertex)->GetPosition()[d]<=(T)0 || (*itr_vertex)->GetPosition()[d]>=(T)size[d])
			{
				(*itr_vertex)->is_boundary_=true;
			}
		}
	}
}

/*
void DynamicTriangularSurface::DrawVertexNormals()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawNormal();
	}
}

void DynamicTriangularSurface::DrawVertexDeviation()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawDeviation();
	}
}

void DynamicTriangularSurface::DrawFaceNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawNormal();
	}
}

void DynamicTriangularSurface::DrawCurvatureNormal(const T& scalar)
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawCurvatureNormal(scalar);
	}
}

void DynamicTriangularSurface::DrawVertexVelocity()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawVelocity();
	}
}

void DynamicTriangularSurface::DrawEdges()
{
	for_each(triangles_.begin(), triangles_.end(), std::mem_fun(&DynamicTriangle::DrawEdges));
}// Remove () from function name to use mem_fun

void DynamicTriangularSurface::DrawTriangles(const bool& draw_front)
{
	if(draw_front)
	{
		TRAVERSE_TRIANGLES(*itr_triangle)->Draw();
	}
	else
	{
		TRAVERSE_TRIANGLES(*itr_triangle)->DrawBack();
	}
}
*/

/*
void DynamicTriangularSurface::DrawTrianglesColor(const bool& draw_front)
{
	if(draw_front == true)
	{
		glBegin(GL_TRIANGLES);
		{
			TRAVERSE_TRIANGLES
			{
				DynamicTriangle& triangle = *(*itr_triangle);

				glColor3fv(triangle.vertices_[0]->curvature_normal_.values_);
				glNormal3fv(triangle.vertices_[0]->GetNormal());
				glVertex3fv(triangle.vertices_[0]->GetPosition());

				glColor3fv(triangle.vertices_[1]->curvature_normal_.values_);
				glNormal3fv(triangle.vertices_[1]->GetNormal());
				glVertex3fv(triangle.vertices_[1]->GetPosition());

				glColor3fv(triangle.vertices_[2]->curvature_normal_.values_);
				glNormal3fv(triangle.vertices_[2]->GetNormal());
				glVertex3fv(triangle.vertices_[2]->GetPosition());
			}
		}
		glEnd();
	}
	else // draw_front == false
	{
		glBegin(GL_TRIANGLES);
		{
			TRAVERSE_TRIANGLES
			{
				DynamicTriangle& triangle = *(*itr_triangle);

				glColor3fv(triangle.vertices_[0]->curvature_normal_.values_);
				glNormal3fv(triangle.vertices_[0]->GetNormal());
				glVertex3fv(triangle.vertices_[0]->GetPosition());

				glColor3fv(triangle.vertices_[2]->curvature_normal_.values_);
				glNormal3fv(triangle.vertices_[2]->GetNormal());
				glVertex3fv(triangle.vertices_[2]->GetPosition());

				glColor3fv(triangle.vertices_[1]->curvature_normal_.values_);
				glNormal3fv(triangle.vertices_[1]->GetNormal());
				glVertex3fv(triangle.vertices_[1]->GetPosition());
			}
		}
		glEnd();
	}
}
*/

/*
void DynamicTriangularSurface::DrawTrianglesNeighborConnectivity()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawNeighborConnectivity();
	}
}
*/

/*
void DynamicTriangularSurface::DrawTrianglesCenter()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DrawCenter();
	}
}
*/

/*
void DynamicTriangularSurface::DrawVerticesNeighborConnectivity()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DrawNeighborConnectivity();
	}
}
*/

/*
void DynamicTriangularSurface::DrawHoles()
{
	std::list<DynamicHole*>::iterator itr_hole;
	for(itr_hole=holes_.begin();itr_hole!=holes_.end();itr_hole++)
	{
		(*itr_hole)->Draw();
	}
}
*/

/*
void DynamicTriangularSurface::DrawHoles(int index)
{
	std::list<DynamicHole*>::iterator itr_hole = holes_.begin();
	for(int i=0;i<index;i++)
	{
		itr_hole++;
	}
	(*itr_hole)->Draw();	
}
*/

void DynamicTriangularSurface::RemoveVertexDeviation()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->x_ += (*itr_vertex)->deviation_;
		for(int d=0;d<3;d++)
		{
			(*itr_vertex)->deviation_.values_[d]=(T)0;
		}
	}
}

void DynamicTriangularSurface::RemoveSmallTriangles(T th)
{
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle=*itr_triangle;
	
		if(triangle->triangles_[0]==NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles_[0]->triangles_[0]==NULL)continue;
			if(triangle->triangles_[0]->triangles_[1]==NULL)continue;
			if(triangle->triangles_[0]->triangles_[2]==NULL)continue;
		}
		if(triangle->triangles_[1]==NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles_[1]->triangles_[0]==NULL)continue;
			if(triangle->triangles_[1]->triangles_[1]==NULL)continue;
			if(triangle->triangles_[1]->triangles_[2]==NULL)continue;
		}		
		if(triangle->triangles_[2]==NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles_[2]->triangles_[0]==NULL)continue;
			if(triangle->triangles_[2]->triangles_[1]==NULL)continue;
			if(triangle->triangles_[2]->triangles_[2]==NULL)continue;
		}		

		TV3 l0,l1,l2;
		l0 = triangle->vertices_[0]->x_ - triangle->vertices_[1]->x_;
		l1 = triangle->vertices_[1]->x_ - triangle->vertices_[2]->x_;
		l2 = triangle->vertices_[2]->x_ - triangle->vertices_[0]->x_;

		T det0,det1,det2;
		det0 = l0.getMagnitude();
		det1 = l1.getMagnitude();
		det2 = l2.getMagnitude();
//		ARRAY_VECTOR3::det<T>(l0.values_,&det0);
//		ARRAY_VECTOR3::det<T>(l1.values_,&det1);
//		ARRAY_VECTOR3::det<T>(l2.values_,&det2);

		if(det0<th && det1<th && det2<th)
		{
			DynamicTriangularSurface::Collapse(triangle);			
			itr_triangle=triangles_.begin();
		}
	}
}

void DynamicTriangularSurface::RemoveLongTriangles(T r)
{
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle=*itr_triangle;
		
		// ?????? ?׵θ??? ?پ????? ??� ?????Ѵ?.
		if(triangle->triangles_[0]==NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles_[0]->triangles_[0]==NULL)continue;
			if(triangle->triangles_[0]->triangles_[1]==NULL)continue;
			if(triangle->triangles_[0]->triangles_[2]==NULL)continue;
		}
		if(triangle->triangles_[1]==NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles_[1]->triangles_[0]==NULL)continue;
			if(triangle->triangles_[1]->triangles_[1]==NULL)continue;
			if(triangle->triangles_[1]->triangles_[2]==NULL)continue;
		}		
		if(triangle->triangles_[2] == NULL)
		{
			continue;
		}
		else
		{
			if(triangle->triangles_[2]->triangles_[0]==NULL)continue;
			if(triangle->triangles_[2]->triangles_[1]==NULL)continue;
			if(triangle->triangles_[2]->triangles_[2]==NULL)continue;
		}		

		TV3 l0,l1,l2;
		l0 = triangle->vertices_[0]->x_ - triangle->vertices_[1]->x_;
		l1 = triangle->vertices_[1]->x_ - triangle->vertices_[2]->x_;
		l2 = triangle->vertices_[2]->x_ - triangle->vertices_[0]->x_;

		T det0,det1,det2;
		det0 = l0.getMagnitude();
		det1 = l1.getMagnitude();
		det2 = l2.getMagnitude();
//		ARRAY_VECTOR3::det<T>(l0.values_,&det0);
//		ARRAY_VECTOR3::det<T>(l1.values_,&det1);
//		ARRAY_VECTOR3::det<T>(l2.values_,&det2);

		T max=det0;
		if(det1>max)
		{
			max=det1;
		}
		if(det2>max)
		{
			max=det2;
		}

		T temp=det0+det1+det2-max;

		T ratio=max/temp;

		if(ratio>(T)0.99)
		{			
			Collapse(triangle);
			itr_triangle=triangles_.begin();
		}
	}
}

void DynamicTriangularSurface::RemoveSmallEdges(T th)
{
	TRAVERSE_TRIANGLES
	{
//		TRIANGLE *begin=(*triangles.begin());
//		TRIANGLE *end=(*triangles.end());
		DynamicTriangle *triangle=*itr_triangle;

		// boundary?? triangle??� �?????? ?ʴ´?.
		bool skip=false;
		for(int i=0;i<3;i++)
		{
			if(triangle->triangles_[i]==NULL)
			{
				skip=true;
				break;
			}
			else
			{
				for(int j=0;j<3;j++)
				{
					if(triangle->triangles_[i]->triangles_[j]==NULL)
					{
						skip=true;
						break;
					}
				}
			}
			if(skip==true)break;
		}
		if(skip==true)continue;

		for(int i=0;i<3;i++)
		{
			TV3 l;// edge 0 vector
			l = triangle->vertices_[i]->x_ - triangle->vertices_[(i+1)%3]->x_; 
			T det;// length of edge 0
//			ARRAY_VECTOR3::det<T>(l.values_,&det);
			det = l.getMagnitude();
			if(det<th)
			{	
				if(itr_triangle!=triangles_.begin())
				{
					itr_triangle--;
					Collapse(triangle->vertices_[i],triangle->vertices_[(i+1)%3]);
				}
				else
				{
					Collapse(triangle->vertices_[i],triangle->vertices_[(i+1)%3]);
					itr_triangle=triangles_.begin();
				}
			}
			continue;
		}
	}
}


void DynamicTriangularSurface::AverageDuplexPositionNormal(DynamicTriangularSurface* neighbor, float dx)
{
	T critical_value = dx/(T)4.0;
	for(unsigned int i=0; i < vertices_.size(); i++)
	{	
		DynamicVertex* i_vertex = vertices_[i];
	
		for(unsigned int j=0; j < neighbor->vertices_.size(); j++)
		{
			DynamicVertex* j_vertex = neighbor->vertices_[j];
			TV3 dis(i_vertex->x_ - j_vertex->x_); 
			if(abs(dis.values_[0]) < critical_value && abs(dis.values_[1]) < critical_value && abs(dis.values_[2]) < critical_value) 
			{
				TV3 avrNormal((i_vertex->n_ + j_vertex->n_) * 0.5f);
				i_vertex->n_ = avrNormal;
				j_vertex->n_ = avrNormal;
			}
		}
	}
}

void DynamicTriangularSurface::Read(const char *filename)
{
	std::ifstream ist(filename, std::ios::binary);

	if(ist.bad())
	{
		std::cout<<"Failed to open file: "<<filename<<std::endl;
	}
	int number_of_vertices=0;
	int number_of_triangles=0;

	// read vertex positions and normals
	ist.read((char*)&number_of_vertices,sizeof(int));
	std::cout<<"Reading vertices "<<number_of_vertices<<std::endl;
	for(int i=0;i<number_of_vertices;i++)
	{
		T position[3],normal[3],velocity[3];
		ist.read((char*)&position[0],sizeof(T));
		ist.read((char*)&position[1],sizeof(T));
		ist.read((char*)&position[2],sizeof(T));
		ist.read((char*)&normal[0],sizeof(T));
		ist.read((char*)&normal[1],sizeof(T));
		ist.read((char*)&normal[2],sizeof(T));
		ist.read((char*)&velocity[0],sizeof(T));
		ist.read((char*)&velocity[1],sizeof(T));
		ist.read((char*)&velocity[2],sizeof(T));		
		AddVertex(position,normal);
	}

	// read triangles
	ist.read((char*)&number_of_triangles,sizeof(int));
	std::cout<<"Reading triangles "<<number_of_triangles<<std::endl;
	for(int i=0;i<number_of_triangles;i++)
	{
		int index[3];
		ist.read((char*)&index[0],sizeof(int));
		ist.read((char*)&index[1],sizeof(int));
		ist.read((char*)&index[2],sizeof(int));
		AddTriangle(index[0],index[1],index[2]);
	}

	ist.close();
}

void DynamicTriangularSurface::ReadSMF(const char *filename)
{
	std::cout<<"# Reading "<<filename<<std::endl;
	std::ifstream file(filename);
	char c[255];
	int frame=0;
	int flag=0;
	while(1)
	{
		frame++;
		if(frame>20000)
		{
			std::cout << ".";
			frame=0;
		}
		file>>c;
		if(file.eof()!=0)
		{
			break;
		}
		if(strcmp(c,"#")==0)
		{
			file.getline(c,255);// ?ּ? ó??, ?ּ?� 255?? ?̳?
		}
		else if(strcmp(c,"v")==0)
		{
			if(flag==0)
			{
				std::cout << "v"; flag++;
			}
			T x[3];
			file>>x[0]>>x[1]>>x[2];
			DynamicTriangularSurface::AddVertex(x);
		}
		else if(strcmp(c,"f")==0)
		{	
			if(flag==1)
			{
				std::cout << "f"; flag++;
			}
			int v[3];
			file>>v[0]>>v[1]>>v[2];
			v[0]--;
			v[1]--;
			v[2]--;
			DynamicTriangle *triangle=AddTriangle(vertices_[v[0]],vertices_[v[1]],vertices_[v[2]]);
			triangle->DetermineNormal();
		}
	}

	DetVertexNormals();
	file.close();

	std::cout << "# Finished" << std::endl;
}

void DynamicTriangularSurface::GetGeometryInfo(std::vector<TV3> &polygon_vertices,std::vector<TV3_INT> &faces)
{
	int index;
	for(index=1;index<=(int)DynamicTriangularSurface::vertices_.size();index++)
	{
		DynamicTriangularSurface::vertices_[index-1]->index_=index;
		polygon_vertices.push_back(TV3(vertices_[index-1]->x_.values_[0],vertices_[index-1]->x_.values_[1],vertices_[index-1]->x_.values_[2]));
	}

	TRAVERSE_TRIANGLES
	{
		faces.push_back(TV3_INT((*itr_triangle)->vertices_[0]->index_,(*itr_triangle)->vertices_[1]->index_,(*itr_triangle)->vertices_[2]->index_));
	}
}

void DynamicTriangularSurface::DertermineTriangleNormal()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DetermineNormal();
	}
}

void DynamicTriangularSurface::ArrangeTriangles()
{
	triangle_array_.initialize((int)triangles_.size());

	int index = 0;
	TRAVERSE_TRIANGLES
	{
		triangle_array_[index++] = *itr_triangle;
	}
}

void DynamicTriangularSurface::Write(const char *filename)
{
	std::ofstream ost(filename, std::ios::binary);
	if(!ost.is_open())
	{
		std::cout<<"Failed to write "<<filename<<std::endl;
	}
	int number_of_vertices=(int)vertices_.size();
	int number_of_triangles=(int)triangles_.size();

	// write vertex positions and normals
	std::cout<<"Writing vertices "<<number_of_vertices<<std::endl;
	ost.write((char*)&number_of_vertices,sizeof(int));
	for(int index=0;index<number_of_vertices;index++)
	{
		DynamicVertex *v=vertices_[index];
		v->index_=index;// set index of each vertex (from 0)
		ost.write((char*)&v->x_.values_[0],sizeof(T));
		ost.write((char*)&v->x_.values_[1],sizeof(T));
		ost.write((char*)&v->x_.values_[2],sizeof(T));
		ost.write((char*)&v->n_.values_[0],sizeof(T));
		ost.write((char*)&v->n_.values_[1],sizeof(T));
		ost.write((char*)&v->n_.values_[2],sizeof(T));
		ost.write((char*)&v->velocity_.values_[0],sizeof(T));
		ost.write((char*)&v->velocity_.values_[1],sizeof(T));
		ost.write((char*)&v->velocity_.values_[2],sizeof(T));
	}

	// write triangles
	std::cout<<"Writing triangles "<<number_of_triangles<<std::endl;
	ost.write((char*)&number_of_triangles,sizeof(int));
	std::list<DynamicTriangle*>::iterator itr_triangle;
	for(itr_triangle=triangles_.begin();itr_triangle!=triangles_.end();itr_triangle++)
	{
		ost.write((char*)&(*itr_triangle)->vertices_[0]->index_,sizeof(int));
		ost.write((char*)&(*itr_triangle)->vertices_[1]->index_,sizeof(int));
		ost.write((char*)&(*itr_triangle)->vertices_[2]->index_,sizeof(int));
	}

	ost.close();
}

void DynamicTriangularSurface::WriteOBJ(const char *filename, TV3& position, QUATERNION& quat)
{
	std::cout<<"# Writing "<<filename<<std::endl;
	std::ofstream file(filename);
//	int frame=0;
	int index;
	
	TV3 vertex_position;
	TV3 vertex_normal;
	TV3 vertex_velocity;

	for(index=1;index<=(int)vertices_.size();index++)
	{
//		VERTEX *v=vertices[index-1];
		DynamicTriangularSurface::vertices_[index-1]->index_=index;

		TV3 x(vertices_[index-1]->x_);
		vertex_position = quat.Rotate(x);
		vertex_position = vertex_position+position;

		//file<<"v "<<vertices[index-1]->x_[0]<<" "<<vertices[index-1]->x_[1]<<" "<<vertices[index-1]->x_[2]<<endl;
		file << "v " << vertex_position.x_ << " " << vertex_position.y_ << " " << vertex_position.z_ << std::endl;
	}
	
	for(index=1;index<=(int)vertices_.size();index++)
	{
		file << "vt " << "1" << " " << "0" << std::endl;
	}

	for(index=1;index<=(int)vertices_.size();index++)
	{
		TV3 x(vertices_[index-1]->n_);
		vertex_normal = quat.Rotate(x);

		file << "vn " << vertex_normal.x_ << " " << vertex_normal.y_ << " " << vertex_normal.z_ << std::endl;
	}
	for(index=1;index<=(int)vertices_.size();index++)
	{
		T scale=(T)0.1;
		TV3 x(vertices_[index-1]->velocity_);
		vertex_velocity = quat.Rotate(x);
		//file<<"vv "<<vertices[index-1]->velocity[0]*scale<<" "<<vertices[index-1]->velocity[1]*scale<<" "<<vertices[index-1]->velocity[2]*scale<<endl;
		file << "vv " << vertex_velocity.x_*scale << " " << vertex_velocity.y_*scale << " " << vertex_velocity.z_*scale << std::endl;
	}

	TRAVERSE_TRIANGLES
	{
		int index[3]={(*itr_triangle)->vertices_[0]->index_,(*itr_triangle)->vertices_[1]->index_,(*itr_triangle)->vertices_[2]->index_};
		file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << std::endl;
	}

	file.close();
}

void DynamicTriangularSurface::WriteOBJ(Array1D<DynamicTriangularSurface*> &arr, const char *filename)
{
	std::cout<<"# Writing "<<filename<<std::endl;
	std::ofstream file(filename);
	
	Array1D<int> vertex_count(arr.num_elements_);

	int vertex_counter = 0;

	for(int i=0;i<arr.num_elements_;i++)
	{
		DynamicTriangularSurface *triangular_surface = arr.values_[i];
		int vertices_size = (int)triangular_surface->vertices_.size();
		int index;

		for(index=1;index<=vertices_size;index++)
		{
	//		VERTEX *v=vertices[index-1];
			file << "v " << triangular_surface->vertices_[index - 1]->x_.values_[0] << " " << triangular_surface->vertices_[index - 1]->x_.values_[1] << " " << triangular_surface->vertices_[index - 1]->x_.values_[2] << std::endl;
		}		
	}

	for(int i=0;i<arr.num_elements_;i++)
	{
		DynamicTriangularSurface *triangular_surface = arr.values_[i];
		int vertices_size = (int)triangular_surface->vertices_.size();
		int index;

		for(index=1;index<=vertices_size;index++)
		{
			file << "vt " << "1" << " " << "0" << std::endl;
		}
	}

	for(int i=0;i<arr.num_elements_;i++)
	{
		DynamicTriangularSurface *triangular_surface = arr.values_[i];
		int vertices_size = (int)triangular_surface->vertices_.size();
		int index;

		for(index=1;index<=vertices_size;index++)
		{
			file << "vn " << triangular_surface->vertices_[index - 1]->n_.values_[0] << " " << triangular_surface->vertices_[index - 1]->n_.values_[1] << " " << triangular_surface->vertices_[index - 1]->n_.values_[2] << std::endl;
		}
	}

	for(int i=0;i<arr.num_elements_;i++)
	{
		DynamicTriangularSurface *triangular_surface = arr.values_[i];
		int vertices_size = (int)triangular_surface->vertices_.size();
		int index;

		for(index=1;index<=vertices_size;index++)
		{
			T scale=(T)0.1;
			file << "vv " << triangular_surface->vertices_[index - 1]->velocity_.values_[0] * scale << " " << triangular_surface->vertices_[index - 1]->velocity_.values_[1] * scale << " " << triangular_surface->vertices_[index - 1]->velocity_.values_[2] * scale << std::endl;
		}
	}
	
	for(int i=0;i<arr.num_elements_;i++)
	{
		DynamicTriangularSurface *triangular_surface = arr.values_[i];
		std::list <DynamicTriangle*>::iterator itr_triangle;
		for(itr_triangle = triangular_surface->triangles_.begin(); itr_triangle != triangular_surface->triangles_.end(); itr_triangle ++)
		{	
			int index[3]={(*itr_triangle)->vertices_[0]->index_,
							(*itr_triangle)->vertices_[1]->index_,
							(*itr_triangle)->vertices_[2]->index_};
			file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << std::endl;
		}
	}
	file.close();
}
//Realflow Bin (mesh) file exporter
//Eulerian or Sph에서 생성한 fluid surface를 export한다.
//i_frame : 현재 frame number
void DynamicTriangularSurface::WriteBinMeshFile(Array1D<DynamicTriangularSurface*> &surface_arr, int i_frame, const char *dir_path)
{
	//file name, 현재 bin file이 저장되는 폴더와 파일의 이름이 dir_path로 같다.
	std::string file_name;
	std::string folder_name(dir_path);//TODO : folder name scripting
	file_name.append(folder_name.c_str());
	file_name.append("\\"+std::string(dir_path)+"_");//TODO : file name scripting
	file_name.append(boost::str(boost::format("%05d")%i_frame));
	file_name.append(".bin");

	std::cout<<"bin:"<<dir_path<<std::endl;//debug

	if(_access(folder_name.c_str(), 0)==-1)//folder not exist
	{
		_mkdir(folder_name.c_str());
	}
	FILE* fp = fopen(file_name.c_str(),"wb");	
	if(fp==NULL)
	{
		std::cout<<i_frame<<"RF Bin (mesh) file pointer is NULL"<<std::endl; 
		return; 
	}
		

	//HEADER
	unsigned int ID_code = 0xDADADADA;
	fwrite(&ID_code, sizeof(unsigned int), 1, fp);

	unsigned int version = 4;
	fwrite(&version, sizeof(unsigned int), 1, fp);

	unsigned int chunk_code = 0xCCCCCCCC;
	fwrite(&chunk_code, sizeof(unsigned int), 1, fp);


	//1. PARTICLE VTX POSITIONS
	//!!!! x, z좌표가 바뀌어야 정상적으로 돌아감 !!!!	
	int n_tot_vtx = 0;
	for(int i_sur=0; i_sur<surface_arr.num_elements_; i_sur++)//각 surface
	{
		n_tot_vtx += (int)surface_arr[i_sur]->vertices_.size();
	}		
	fwrite(&n_tot_vtx, sizeof(int), 1, fp);


	int n_vtx;//각 surface에서의 vtx갯수
	float x, y, z;
	for(int i_sur=0; i_sur<surface_arr.num_elements_; i_sur++)//각 surface	
	{
		DynamicTriangularSurface *surface = surface_arr[i_sur];
		n_vtx = (int)surface->vertices_.size();
		
		for(int i_vtx=0; i_vtx<n_vtx; i_vtx++)
		{
			surface->vertices_[i_vtx]->index_ = i_vtx;
			x = (float)surface->vertices_[i_vtx]->x_.values_[0];
			y = (float)surface->vertices_[i_vtx]->x_.values_[1];
			z = (float)surface->vertices_[i_vtx]->x_.values_[2];
			fwrite(&x, sizeof(float), 1, fp);
			fwrite(&y, sizeof(float), 1, fp);
			fwrite(&z, sizeof(float), 1, fp);
		}
	}	


	//2. PARTICLE FACES
	int n_tot_face = 0;
	for(int i_sur=0; i_sur<surface_arr.num_elements_; i_sur++)
	{
		n_tot_face += (int)surface_arr[i_sur]->triangles_.size();
	}		               
	fwrite(&n_tot_face, sizeof(int), 1, fp);
	

	int n_accum_vtx=0;
	int i, j, k;
	for(int i_sur=0; i_sur<surface_arr.num_elements_; i_sur++)//각 surface
	{
		DynamicTriangularSurface *surface = surface_arr[i_sur];
		
		std::list <DynamicTriangle*>::iterator itr_tri;
		for(itr_tri = surface->triangles_.begin(); itr_tri != surface->triangles_.end(); itr_tri++)		
		{	
			i = (*itr_tri)->vertices_[0]->index_ + n_accum_vtx;
			j = (*itr_tri)->vertices_[1]->index_ + n_accum_vtx;
			k = (*itr_tri)->vertices_[2]->index_ + n_accum_vtx;
			fwrite(&i, sizeof(int), 1, fp);
			fwrite(&j, sizeof(int), 1, fp);
			fwrite(&k, sizeof(int), 1, fp);
		}
		n_accum_vtx += (int)surface->vertices_.size();
	}


	//3. FLUID TEXTURE
	//   혼합 유체 TEXTURING관련인듯. SKIP..
	/*[unsigned int] ; texture chunk code = 0xCCCCCC00 (**)
		[int] ; number of fluids
		loop for [number of vertices]
	loop for [number of fluids-1] ; version>=3 (***)
		[float] ; texture weight (***)
		endloop
		[float] ; X texture coordinate
		[float] ; Y texture coordinate
		[float] ; Z texture coordinate
		endloop
	//*/	
	
	
	//4. VTX VELOCITY
	unsigned int vel_chunk_code = 0xCCCCCC11;
	fwrite(&vel_chunk_code, sizeof(unsigned int), 1, fp);

	float vel_x, vel_y, vel_z;
	for(int i_sur=0; i_sur<surface_arr.num_elements_; i_sur++)//각 surface
	{
		DynamicTriangularSurface *surface = surface_arr[i_sur];
		n_vtx = (int)surface->vertices_.size();

		for(int i_vtx=0; i_vtx<n_vtx; i_vtx++)
		{			
			vel_x = (float)surface->vertices_[i_vtx]->velocity_.values_[0];
			vel_y = (float)surface->vertices_[i_vtx]->velocity_.values_[1];
			vel_z = (float)surface->vertices_[i_vtx]->velocity_.values_[2];
			fwrite(&vel_x, sizeof(float), 1, fp);
			fwrite(&vel_y, sizeof(float), 1, fp);
			fwrite(&vel_z, sizeof(float), 1, fp);
		}			
			
	}
	
	//end of file mark
	unsigned int eof_mark = 0xDEDEDEDE;
	fwrite(&eof_mark, sizeof(unsigned int), 1, fp);
	
	fflush(fp);
	fclose(fp);
}

void DynamicTriangularSurface::WritePOLY(const char* filename)
{
	std::cout<<"# Writing "<<filename<<std::endl;
	std::ofstream file(filename);

	// #Part 1 - node list
	file << "#Part 1 - node list" << std::endl;
	int node_count = (int)vertices_.size();
	int dimention = 3;
	int attribute = 0;
	int boundary_marker = 0;

	file << node_count << " " << dimention << " " << attribute << " " << boundary_marker << std::endl;

	for(int i=0; i<node_count; i++)
		file << (i+1) << " " << vertices_[i]->x_.values_[0] << " " << vertices_[i]->x_.values_[1] << " " << vertices_[i]->x_.values_[2] << std::endl;	

	// #Part 2 - facet list
	file << "#Part 2 - facet list" << std::endl;
	int facet_count = (int)triangles_.size();

	file << facet_count << " " << boundary_marker << std::endl;

	TRAVERSE_TRIANGLES
	{
		file << 1 << " " << 0 << " " << 1 << std::endl;
		file << 3 << " ";

		int index[3]={(*itr_triangle)->vertex_indices_.i_,(*itr_triangle)->vertex_indices_.j_, (*itr_triangle)->vertex_indices_.k_};
		file << index[0]+1 << " " << index[1]+1 << " " << index[2]+1 << std::endl;
	}


	// # Part 3 - hole list
	file << "# Part 3 - hole list" << std::endl;
	file << 0 << std::endl;


	// # Part 4 - region list
	file << "# Part 4 - region list" << std::endl;
	file << 0 << std::endl;

	file.close();
}

void DynamicTriangularSurface::WriteOBJ(const char *filename)
{
	std::cout<<"# Writing "<<filename<<std::endl;
	std::ofstream file(filename);
//	int frame=0;
	int index;
	for(index=1;index<=(int)vertices_.size();index++)
	{
//		VERTEX *v=vertices[index-1];
		DynamicTriangularSurface::vertices_[index-1]->index_=index;
		file << "v " << vertices_[index - 1]->x_.values_[0] << " " << vertices_[index - 1]->x_.values_[1] << " " << vertices_[index - 1]->x_.values_[2] << std::endl;
	}
	
	for(index=1;index<=(int)vertices_.size();index++)
	{
		file << "vt " << "1" << " " << "0" << std::endl;
	}
	for(index=1;index<=(int)vertices_.size();index++)
	{
		file << "vn " << vertices_[index - 1]->n_.values_[0] << " " << vertices_[index - 1]->n_.values_[1] << " " << vertices_[index - 1]->n_.values_[2] << std::endl;
	}
	for(index=1;index<=(int)vertices_.size();index++)
	{
		T scale=(T)0.1;
		file << "vv " << vertices_[index - 1]->velocity_.values_[0] * scale << " " << vertices_[index - 1]->velocity_.values_[1] * scale << " " << vertices_[index - 1]->velocity_.values_[2] * scale << std::endl;
	}
	TRAVERSE_TRIANGLES
	{
		int index[3]={(*itr_triangle)->vertices_[0]->index_,(*itr_triangle)->vertices_[1]->index_,(*itr_triangle)->vertices_[2]->index_};
		file << "f " << index[0] << "/" << index[0] << "/" << index[0] << " " << index[1] << "/" << index[1] << "/" << index[1] << " " << index[2] << "/" << index[2] << "/" << index[2] << std::endl;
	}

	file.close();
}

void DynamicTriangularSurface::ReadOBJ(const char *filename, const T angle, const TV3 axis, const TV3 scale, const TV3 translation)
{
	std::ifstream file(filename);

	if (file.good() == false)
	{
		std::cout<<filename<<" does not exist. Program terminated."<<std::endl;
		exit(-1);
	}

	std::vector<TV3> normals_temp;
	std::vector<TV3> texture_coordinates_temp;

	char c[255];
	
	bool flag_v(false);
	bool flag_vt(false);
	bool flag_vn(false);
	bool flag_f(false);

	bool is_first_vertex(true);

	QUATERNION quaternion(-angle*(T)0.0174532925,axis);//degree to radian
	
	while(1)
	{
		file >> c;

		if(file.eof() != 0)
		{
			break;
		}

		if(strcmp(c,"#") == 0)
		{
			file.getline(c,255); // comments (less than 255 characters)
		}
		else if(strcmp(c,"v") == 0) // vertices
		{
			if(flag_v == false)
			{
				flag_v = true;
			}
			TV3 vertex_pos;
			file >> vertex_pos.x_ >> vertex_pos.y_ >> vertex_pos.z_;

			if(angle != (T)0)
			{
				vertex_pos = quaternion.Rotate(vertex_pos);
			}
			
			vertex_pos.x_*=scale.x_;
			vertex_pos.y_*=scale.y_;
			vertex_pos.z_*=scale.z_;

			vertex_pos = vertex_pos + translation;

			DynamicVertex *vertex_temp = new DynamicVertex(vertex_pos.values_);
			vertex_temp->is_mc_vertex_ = true;
			vertices_.push_back(vertex_temp);

			if(is_first_vertex == true)
			{
				is_first_vertex = false;
//				bounding_box_.Initialize(vertex_pos, vertex_pos);
				bounding_box_ = BOX_3D<T>(vertex_pos, vertex_pos);
			}
			else
			{
				bounding_box_.Include(vertex_pos);
			}
		}
		else if(strcmp(c,"vt") == 0) // texture coordinate
		{	
			if(flag_vt == false)
			{
				flag_vt = true;
			}

			TV3 texure_vt;
			file >> texure_vt.x_ >> texure_vt.y_;

			texture_coordinates_temp.push_back(texure_vt);			
		}
		else if(strcmp(c,"vn") == 0) // vertex normal
		{	
			if(flag_vn==false)
			{
				flag_vn=true;
			}

			TV3 vertex_n;
			if(flag_vn)
			{
				file >> vertex_n.x_ >> vertex_n.y_ >> vertex_n.z_;
			}

			if(angle!=(T)0)
			{
				vertex_n = quaternion.Rotate(vertex_n);
			}

			normals_temp.push_back(vertex_n);
		}
		else if(strcmp(c,"f") == 0)
		{
			if(flag_f == false)
			{
				flag_f = true;
			}

			int v[3],vt[3],vn[3];
			if(flag_vt == true && flag_vn == true)
			{
				for(int i=0;i<3;i++)
				{
					file>>v[i];file.get(c,2);
					file>>vt[i];file.get(c,2);
					file>>vn[i];

					v[i]--;
					vt[i]--;
					vn[i]--;
				}
			}
			else if(flag_vt==false && flag_vn==true)
			{
				for(int i=0;i<3;i++)
				{
					file>>v[i];file.get(c,2);file.get(c,2);
					file>>vn[i];
					v[i]--;
					vn[i]--;
				}
			}
			else if(flag_vt==false && flag_vn==false)
			{
				for(int i=0;i<3;i++)
				{
					file>>v[i];					
					v[i]--;
				}
			}			

			// add triangle
			DynamicTriangle* triangle=DynamicTriangularSurface::AddTriangle(vertices_[v[0]], vertices_[v[1]], vertices_[v[2]]);
			triangle->vertex_indices_.i_ = v[0];
			triangle->vertex_indices_.j_ = v[1];
			triangle->vertex_indices_.k_ = v[2];

			if(flag_vt==true)
			{
				// set texture coordinate for vertex
				triangle->uv_[0].x_ = texture_coordinates_temp[vt[0]].x_; triangle->uv_[0].y_ = texture_coordinates_temp[vt[0]].y_;
				triangle->uv_[1].x_ = texture_coordinates_temp[vt[1]].x_; triangle->uv_[1].y_ = texture_coordinates_temp[vt[1]].y_;
				triangle->uv_[2].x_ = texture_coordinates_temp[vt[2]].x_; triangle->uv_[2].y_ = texture_coordinates_temp[vt[2]].y_; 
			}

			if(flag_vn==true)
			{
				TV3 n;
//				ARRAY_VECTOR3::set<T>(n.values_,(T)0);
				n = TV3((T)0,(T)0,(T)0);
				n += normals_temp[vn[0]];
				n += normals_temp[vn[1]];
				n += normals_temp[vn[2]];
				n.normalize();
				triangle->SetNormal(n.values_);
			}
			else
			{
				triangle->DetermineNormal();
			}
		}
	}

	DynamicTriangularSurface::DetVertexNormals();

	file.clear();
	file.close();
	
}

void DynamicTriangularSurface::DetVertexNormals()
{
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DetermineNormal();
	}
}

void DynamicTriangularSurface::DetermineNormalDeviation()
{
	TRAVERSE_VERTICES
	{
		T detdeviation=(*itr_vertex)->deviation_.getMagnitude();
		if(dotProduct((*itr_vertex)->n_,(*itr_vertex)->deviation_) < (T)0)
		{
			detdeviation*=-(T)1;		
		}
//		(*itr_vertex)->normal_deviation_=ARRAY_VECTOR3::dot((*itr_vertex)->n_,(*itr_vertex)->deviation_);
		(*itr_vertex)->normal_deviation_=detdeviation;
		(*itr_vertex)->deviation_ = (*itr_vertex)->normal_deviation_*(*itr_vertex)->n_;
	}
}

void DynamicTriangularSurface::DetTextureCoordinates(T *xy, T *uv, T s)
{
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle = *itr_triangle;

		for(int i=0;i<3;i++)
		{
			TV3 &uv = triangle->uv_[i];
			DynamicVertex *vertex = triangle->vertices_[i];
			uv.x_=(vertex->GetPosition()[0]-xy[0])/s+uv.x_;
			uv.y_=(vertex->GetPosition()[1]-xy[1])/s+uv.y_;
			if(uv.x_>(T)1)
			{
				uv.x_=(T)1;
			}
			if(uv.y_>(T)1)
			{
				uv.y_=(T)1;
			}
			if(uv.x_<(T)0)
			{
				uv.x_=(T)0;
			}
			if(uv.y_<(T)0)
			{
				uv.y_=(T)0;
			}
		}
	}
}

void DynamicTriangularSurface::DetTextureCoordinates(T x,T y,T u,T v,T s)
{
	T xy[2]={x,y};
	T uv[2]={u,v};
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle = *itr_triangle;

		for(int i=0;i<3;i++)
		{
			TV3 &uv = triangle->uv_[i];
			DynamicVertex *vertex = triangle->vertices_[i];
			uv.x_=(vertex->GetPosition()[0]-xy[0])/s+uv.x_;
			uv.y_=(vertex->GetPosition()[1]-xy[1])/s+uv.y_;
			if(uv.x_>(T)1)
			{
				uv.x_=(T)1;
			}
			if(uv.y_>(T)1)
			{
				uv.y_=(T)1;
			}
			if(uv.x_<(T)0)
			{
				uv.x_=(T)0;
			}
			if(uv.y_<(T)0)
			{
				uv.y_=(T)0;
			}
		}
	}
}

void DynamicTriangularSurface::DetCurvatureNormals()
{
//	TRAVERSE_VERTICES{ARRAY_VECTOR3::set((*itr_vertex)->curvature_normal_,(T)0);}
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->AddLocalCurvatureNormal();
	}
	TRAVERSE_VERTICES
	{
		(*itr_vertex)->DetCurvatureNormal();
	}
}

void DynamicTriangularSurface::DetFaceNormals()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->DetermineNormal();
	}
}

void DynamicTriangularSurface::Filtering(T lamda)
{
//	T lamda = 0.6307f;
//	T myu = -0.672f;
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles_.size()==0)
		{
			continue;
		}
		if((*itr_vertex)->is_boundary_==true)
		{
			continue;
		}

//		T dx[3]={(T)0,(T)0,(T)0};
//		T w_total=(T)0;
//		T w=(T)0;
//		T x1[3]={(T)0,(T)0,(T)0};

		std::list <DynamicVertex*> v;
		std::list <DynamicVertex*>::iterator itr_v;
		std::list <DynamicTriangle*>::iterator itr_triangle;

		for(itr_triangle=(*itr_vertex)->triangles_.begin();itr_triangle!=(*itr_vertex)->triangles_.end();itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices_[0]);
			v.remove((*itr_triangle)->vertices_[1]);
			v.remove((*itr_triangle)->vertices_[2]);
			v.push_back((*itr_triangle)->vertices_[0]);
			v.push_back((*itr_triangle)->vertices_[1]);
			v.push_back((*itr_triangle)->vertices_[2]);
		}
		v.remove((*itr_vertex));

		TV3 Li((T)0,(T)0,(T)0);
		T E=(T)0;
		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{
			TV3 dist;
			T c=(T)1;
			dist = (*itr_v)->x_ - (*itr_vertex)->x_;
			c = dist.getMagnitude();
//			ARRAY_VECTOR3::sub<T>((*itr_v)->GetPosition(),(*itr_vertex)->GetPosition(),dist);
//			ARRAY_VECTOR3::det<T>(dist,&c);
			E+=c;
			Li += ((*itr_v)->x_ -(*itr_vertex)->x_)/c;
		}

		Li = (T)2 / E * Li;

		(*itr_vertex)->deviation_ -= lamda * Li;

		(*itr_vertex)->x_ += lamda * Li;	

/*		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{
			T dist[3];
			T c=(T)1;
			ARRAY_VECTOR3::sub((*itr_v)->GetPosition(),(*itr_vertex)->GetPosition(),dist);
			ARRAY_VECTOR3::det(dist,&c);
			c=(T)1/c;
			w_total+=c;
			dx[0]+=c*(*itr_v)->GetPosition()[0];
			dx[1]+=c*(*itr_v)->GetPosition()[1];
			dx[2]+=c*(*itr_v)->GetPosition()[2];
		}
		
		dx[0]/=w_total;
		dx[1]/=w_total;
		dx[2]/=w_total;
		dx[0]-=(*itr_vertex)->GetPosition()[0];
		dx[1]-=(*itr_vertex)->GetPosition()[1];
		dx[2]-=(*itr_vertex)->GetPosition()[2];

		(*itr_vertex)->GetPosition()[0]+=lamda*dx[0];
		(*itr_vertex)->GetPosition()[1]+=lamda*dx[1];
		(*itr_vertex)->GetPosition()[2]+=lamda*dx[2];*/
	}
}

void DynamicTriangularSurface::FilteringTangent(T lamda)
{
//	T lamda=0.6307f;
//	T myu=-0.672f;
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles_.size()==0)
		{
			continue;
		}
		if((*itr_vertex)->is_mc_vertex_==true)
		{
			continue;
		}
//		T dx[3]={(T)0,(T)0,(T)0};
//		T w_total=(T)0;
//		T w=(T)0;
//		T x1[3]={(T)0,(T)0,(T)0};

		std::list <DynamicVertex*> v;
		std::list <DynamicVertex*>::iterator itr_v;
		std::list <DynamicTriangle*>::iterator itr_triangle;
		for(itr_triangle=(*itr_vertex)->triangles_.begin();itr_triangle!=(*itr_vertex)->triangles_.end();itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices_[0]);
			v.remove((*itr_triangle)->vertices_[1]);
			v.remove((*itr_triangle)->vertices_[2]);
			v.push_back((*itr_triangle)->vertices_[0]);
			v.push_back((*itr_triangle)->vertices_[1]);
			v.push_back((*itr_triangle)->vertices_[2]);
		}
		v.remove((*itr_vertex));

		TV3 Li((T)0,(T)0,(T)0);
		T E=(T)0;
		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{
			TV3 dist;
			T c=(T)1;
			dist = (*itr_v)->x_ - (*itr_vertex)->x_;
			c = dist.getMagnitude();
//			ARRAY_VECTOR3::sub<T>((*itr_v)->GetPosition(),(*itr_vertex)->GetPosition(),dist);
//			ARRAY_VECTOR3::det<T>(dist,&c);
			E+=c;
			Li += ((*itr_v)->x_-(*itr_vertex)->x_)/c;
		}

		Li = (T)2 / E * Li;

		// Li?? tangent???Ҹ? ?߷�???.
		T temp=dotProduct(Li, (*itr_vertex)->n_);
		TV3 n;
		n = (*itr_vertex)->n_;
		n *= temp;
		Li -= n;
		(*itr_vertex)->x_+=lamda*Li;
	}
}

void DynamicTriangularSurface::FilteringMinimumVariation(T lamda)
{
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles_.size()==0)
		{
			continue;
		}
		if((*itr_vertex)->is_mc_vertex_==true)
		{
			continue;
		}

		std::list <DynamicVertex*> v;
		std::list <DynamicVertex*>::iterator itr_v;
		std::list <DynamicTriangle*>::iterator itr_triangle;
		for(itr_triangle=(*itr_vertex)->triangles_.begin();itr_triangle!=(*itr_vertex)->triangles_.end();itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices_[0]);
			v.remove((*itr_triangle)->vertices_[1]);
			v.remove((*itr_triangle)->vertices_[2]);
			v.push_back((*itr_triangle)->vertices_[0]);
			v.push_back((*itr_triangle)->vertices_[1]);
			v.push_back((*itr_triangle)->vertices_[2]);
		}
		v.remove((*itr_vertex));

		T stress=(T)0;
		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{
			if(dotProduct((*itr_v)->curvature_normal_,(*itr_v)->n_)<(T)0)
			{
				stress -= (*itr_v)->curvature_normal_.getMagnitude();
			}
			else
			{
				stress += (*itr_v)->curvature_normal_.getMagnitude();
			}
		}

		stress/=(T)v.size();

		T thisstress=(T)0;
		if(dotProduct((*itr_vertex)->curvature_normal_,(*itr_vertex)->n_)<(T)0)
		{
			thisstress -= (*itr_vertex)->curvature_normal_.getMagnitude();
		}
		else
		{
			thisstress += (*itr_vertex)->curvature_normal_.getMagnitude();
		}

		T ds=thisstress-stress;

		TV3 dx;
		dx = (*itr_vertex)->curvature_normal_;
		dx.normalize();
		dx *= -ds;
		dx *= lamda;
		
		(*itr_vertex)->x_ += dx;
	}
}

void DynamicTriangularSurface::Filtering2(T lamda)
{
	std::vector<DynamicVertex*>::iterator itr_vertex;
	for(itr_vertex=vertices_.begin();itr_vertex!=vertices_.end();itr_vertex++)
	{
		// tension stress?? ???Ѵ?.
		std::list <DynamicVertex*> v;
		std::list <DynamicVertex*>::iterator itr_v;
		std::list <DynamicTriangle*>::iterator itr_triangle;
		for(itr_triangle=(*itr_vertex)->triangles_.begin();itr_triangle!=(*itr_vertex)->triangles_.end();itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices_[0]);
			v.remove((*itr_triangle)->vertices_[1]);
			v.remove((*itr_triangle)->vertices_[2]);
			v.push_back((*itr_triangle)->vertices_[0]);
			v.push_back((*itr_triangle)->vertices_[1]);
			v.push_back((*itr_triangle)->vertices_[2]);
		}
		v.remove((*itr_vertex));
//		ARRAY_VECTOR3::set<T>((*itr_vertex)->s_.values_,(T)0);
		(*itr_vertex)->s_ = TV3((T)0,(T)0,(T)0);
		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{	
			TV3 dist;
			dist = (*itr_v)->x_ - (*itr_vertex)->x_;
			(*itr_vertex)->s_ = (*itr_vertex)->s_ + dist;
		}

		(*itr_vertex)->s_ /= (T)v.size();
		(*itr_vertex)->stress_ = (*itr_vertex)->s_.getMagnitude();
//		ARRAY_VECTOR3::div<T>((*itr_vertex)->s_.values_,(T)v.size());
//		ARRAY_VECTOR3::det<T>((*itr_vertex)->s_.values_,&(*itr_vertex)->stress_);
		if(((*itr_vertex)->s_.values_[0]*(*itr_vertex)->n_.values_[0]+(*itr_vertex)->s_.values_[1]*(*itr_vertex)->n_.values_[1]+(*itr_vertex)->s_.values_[2]*(*itr_vertex)->n_.values_[2])<(T)0)
		{
			(*itr_vertex)->stress_*=-(T)1;
		}
	}	
	
	for(itr_vertex=vertices_.begin();itr_vertex!=vertices_.end();itr_vertex++)
	{
		if((*itr_vertex)->triangles_.size()==0)
		{
			continue;
		}
//		if((*itr_vertex)->isMCVertex==true)continue;

		// tension stress?? ???Ѵ?.
		std::list <DynamicVertex*> v;
		std::list <DynamicVertex*>::iterator itr_v;
		std::list <DynamicTriangle*>::iterator itr_triangle;
		for(itr_triangle=(*itr_vertex)->triangles_.begin();itr_triangle!=(*itr_vertex)->triangles_.end();itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices_[0]);
			v.remove((*itr_triangle)->vertices_[1]);
			v.remove((*itr_triangle)->vertices_[2]);
			v.push_back((*itr_triangle)->vertices_[0]);
			v.push_back((*itr_triangle)->vertices_[1]);
			v.push_back((*itr_triangle)->vertices_[2]);
		}
		v.remove((*itr_vertex));

		TV3 s_mean((T)0,(T)0,(T)0);
		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{	
//			s_mean+=(*itr_v)->stress;
			s_mean += (*itr_v)->s_;
		}
		s_mean /= (T)v.size();
		TV3 dx;
		dx = (*itr_vertex)->s_ - s_mean;
		(*itr_vertex)->x_ += lamda * dx;
	}
}

void DynamicTriangularSurface::FilteringDeviation(T lamda)
{
	TRAVERSE_VERTICES
	{
		if((*itr_vertex)->triangles_.size()==0)
		{
			continue;
		}
		if((*itr_vertex)->is_boundary_==true)
		{
			continue;
		}

		std::list<DynamicVertex*> v;
		std::list<DynamicVertex*>::iterator itr_v;
		std::list<DynamicTriangle*>::iterator itr_triangle;

		for(itr_triangle=(*itr_vertex)->triangles_.begin();itr_triangle!=(*itr_vertex)->triangles_.end();itr_triangle++)
		{
			v.remove((*itr_triangle)->vertices_[0]);
			v.remove((*itr_triangle)->vertices_[1]);
			v.remove((*itr_triangle)->vertices_[2]);
			v.push_back((*itr_triangle)->vertices_[0]);
			v.push_back((*itr_triangle)->vertices_[1]);
			v.push_back((*itr_triangle)->vertices_[2]);
		}
		v.remove((*itr_vertex));

		T ds=(T)0;
		T dl=(T)0;
		for(itr_v=v.begin();itr_v!=v.end();itr_v++)
		{
			TV3 l;
			l = (*itr_v)->x_ - (*itr_vertex)->x_;
//			ds+=((*itr_v)->normal_deviation_-(*itr_vertex)->normal_deviation_)/ARRAY_VECTOR3::det(l);
			ds+=((*itr_v)->normal_deviation_-(*itr_vertex)->normal_deviation_);
			dl += l.getMagnitude();
		}
		ds/=(T)v.size();
//		ds-=(*itr_vertex)->normal_deviation_;
//		ds/=dl;
		(*itr_vertex)->normal_deviation_+=ds*lamda;
		(*itr_vertex)->deviation_ = (*itr_vertex)->normal_deviation_ * (*itr_vertex)->n_;
	}
}

void DynamicTriangularSurface::FlipEdges()
{
//	int count=0;
	TRAVERSE_TRIANGLES
	{
//		TRIANGLE *ttt=*itr_triangle;
		if((*itr_triangle)->CountFeatureVertex()!=1)
		{
			continue;
		}
		DynamicVertex *v_feature=NULL;
		for(int i=0;i<3;i++)
		{
			if((*itr_triangle)->vertices_[i]->feature_==true)
			{
				v_feature=(*itr_triangle)->vertices_[i];
				break;
			}
		}
		DynamicTriangle *t_neighbor=NULL;
		for(int i=0;i<3;i++)
		{
			if((*itr_triangle)->triangles_[i]!=NULL)
			{
				if((*itr_triangle)->triangles_[i]->CountFeatureVertex()==1)
				{
					if((*itr_triangle)->triangles_[i]->vertices_[0]!=v_feature)
					{
						if((*itr_triangle)->triangles_[i]->vertices_[1]!=v_feature)
						{
							if((*itr_triangle)->triangles_[i]->vertices_[2]!=v_feature)
							{
								t_neighbor = (*itr_triangle)->triangles_[i];break;
							}
						}
					}
				}
			}
		}
		if(t_neighbor==NULL)
		{
			continue;
		}
		DynamicVertex *v[4],*temp;
		v[0]=(*itr_triangle)->vertices_[0];
		v[1]=(*itr_triangle)->vertices_[1];
		v[2]=(*itr_triangle)->vertices_[2];

		// v[0] � feature vertex?? ??????.
		if(v[1]->feature_==true)
		{
			temp=v[1];v[1]=v[0];v[0]=temp;
		}
		if(v[2]->feature_==true)
		{
			temp=v[2];v[2]=v[0];v[0]=temp;
		}
		for(int j=0;j<3;j++)
		{
			if(t_neighbor->vertices_[j]->feature_==true)
			{
				v[3]=t_neighbor->vertices_[j];
				break;
			}
		}
		DynamicTriangle *t=(*itr_triangle);				
		DelTriangle(t);
		DelTriangle(t_neighbor);
		AddTriangle(v[0],v[1],v[3]);
		AddTriangle(v[0],v[2],v[3]);				

		itr_triangle=triangles_.begin();		
//		continue;
	}
/*		if((*itr_triangle)!=NULL)
		if((*itr_triangle)->CountFeatureVertex()==1){
			for(int i=0;i<3;i++)
				if((*itr_triangle)->triangles_[i]!=NULL)
				if((*itr_triangle)->triangles_[i]->CountFeatureVertex()==1){
					TRIANGLE *t0=*itr_triangle;
					TRIANGLE *t1=(*itr_triangle)->triangles_[i];
					VERTEX *v[4],temp;
					v[0]=(*itr_triangle)->vertices_[0];
					v[1]=(*itr_triangle)->vertices_[1];
					v[2]=(*itr_triangle)->vertices_[2];
					// v[0] � feature vertex?? ??????.
					if(v[1]->feature_==true){temp=v[1];v[1]=v[0];v[0]=temp;}
					if(v[2]->feature_==true){temp=v[2];v[2]=v[0];v[0]=temp;}
					for(int j=0;j<3;j++)if(t1->vertices_[j]->feature_==true){v[3]=t1->vertices_[j];break;}
					if(count==4)return;
					if(count==3){int t=0;}
					count++;
					if(v[0]==v[3])continue;
					DelTriangle(t0);
					DelTriangle(t1);
					AddTriangle(v[0],v[1],v[3]);
					AddTriangle(v[0],v[2],v[3]);
					itr_triangle=triangles.begin();
					break;}}}*/
}

void DynamicTriangularSurface::CorrectCCW()
{
	TRAVERSE_TRIANGLES
	{
		(*itr_triangle)->CorrectCCW();
	}
}

void DynamicTriangularSurface::ButterflySubdivision()
{
	std::cout << "Subdivision started" << std::endl;
	std::cout << "triangles are " << (int)triangles_.size() << std::endl;

	std::cout << "making new edge vertices" << std::endl;
	TRAVERSE_TRIANGLES
	{
		int progress_check = 0;
		if(progress_check > 100)
		{
			std::cout << ".";
			progress_check = 0;
		}
		else progress_check ++;

		(*itr_triangle)->is_old_ = true;

		for(int d = 0; d < 3; d++)	// for 3 edges of this triangle
		{
			if((*itr_triangle)->triangles_[0] != NULL && (*itr_triangle)->triangles_[1] != NULL && (*itr_triangle)->triangles_[2])
			{
				// butterfly Subdivision
				DynamicVertex *stencil_vertex[8];
				TV3 x[8];
				TV3 x_weight;
				TV3 n;
//				TRIANGLE *thistriangle=*itr_triangle;
				stencil_vertex[1]=(*itr_triangle)->vertices_[d];
				stencil_vertex[3]=(*itr_triangle)->vertices_[(d+1)%3];
				stencil_vertex[4]=(*itr_triangle)->vertices_[(d+2)%3];
				stencil_vertex[0]=(*itr_triangle)->triangles_[(d+1)%3]->FindAnotherVertex(stencil_vertex[d],stencil_vertex[(d+2)%3]);
				stencil_vertex[2]=(*itr_triangle)->triangles_[(d+2)%3]->FindAnotherVertex(stencil_vertex[d],stencil_vertex[(d+1)%3]);
				stencil_vertex[0]=(*itr_triangle)->triangles_[(d+1)%3]->FindAnotherVertex(stencil_vertex[1],stencil_vertex[4]);
				stencil_vertex[2]=(*itr_triangle)->triangles_[(d+2)%3]->FindAnotherVertex(stencil_vertex[1],stencil_vertex[3]);
				x[1] = stencil_vertex[1]->x_;
				x[3] = stencil_vertex[3]->x_;
				x[4] = stencil_vertex[4]->x_;
				x[0] = stencil_vertex[0]->x_;
				x[2] = stencil_vertex[2]->x_;
				x[1] /= (T)8;
				x[3] /= (T)2*(T)2;
				x[4] /= (T)2*(T)2;
				x[0] /= -(T)16;
				x[2] /= -(T)16;

				x_weight = TV3((T)0,(T)0,(T)0);
//				ARRAY_VECTOR3::set<T>(x_weight.values_,(T)0);

				x_weight += x[1];
				x_weight += x[3];
				x_weight += x[4];
				x_weight += x[0];
				x_weight += x[2];
				n = stencil_vertex[3]->GetNormal();
				n += stencil_vertex[4]->GetNormal();
				n.normalize();

				if((*itr_triangle)->edge_vertex_[d] == NULL)
				{
					(*itr_triangle)->edge_vertex_[d] = AddVertex(x_weight.values_,n.values_);
					(*itr_triangle)->triangles_[d]->edge_vertex_[(*itr_triangle)->triangles_[d]->GetNeighborIndex((*itr_triangle))] = (*itr_triangle)->edge_vertex_[d];
				}
				else
				{
					(*itr_triangle)->edge_vertex_[d]->x_ += x_weight;
				}
			}
			else
			{
				DynamicVertex *v[8];
				TV3 x[8];
				TV3 x_weight;
				TV3 n;
//				TRIANGLE *thistriangle=*itr_triangle;
				v[3]=(*itr_triangle)->vertices_[(d+1)%3];
				v[4]=(*itr_triangle)->vertices_[(d+2)%3];
				x[3] = v[3]->x_;
				x[4] = v[4]->x_;
				x[3] /= (T)2;// ?ٸ? ?ﰢ???? ???ؼ? ?? ?? ?? ?????? ?״ϱ?...
				x[4] /= (T)2;

				x_weight = TV3((T)0,(T)0,(T)0);
//				ARRAY_VECTOR3::set<T>(x_weight.values_,(T)0);

				x_weight += x[3];
				x_weight += x[4];
				n = v[3]->n_;
				n += v[4]->n_;
				n.normalize();
				if((*itr_triangle)->edge_vertex_[d]==NULL)
				{
					(*itr_triangle)->edge_vertex_[d]=AddVertex(x_weight.values_,n.values_);
				}
				else
				{
//					ARRAY_VECTOR3::add<T>((*itr_triangle)->edge_vertex_[d]->GetPosition(), x_weight, (*itr_triangle)->edge_vertex_[d]->GetPosition());
					x_weight = (*itr_triangle)->edge_vertex_[d]->x_;
				}
			}
		}
	}

	std::cout << "splitting triangles" << std::endl;
	std::list<DynamicTriangle*>::iterator itr_triangle2 = triangles_.begin();
	while(itr_triangle2 != triangles_.end())
	{	
		static int progress_check=0;
		if(progress_check > 10000)
		{
			std::cout << ".";
			progress_check=0;
		}
		else
		{
			progress_check++;
		}

		if((*itr_triangle2)->is_old_==false)
		{
			break;
		}

//		if((*itr_triangle2)->edge_vertex_[0]==NULL){T temp=(T)3;}
//		if((*itr_triangle2)->edge_vertex_[1] == NULL){T temp=(T)3;}
//		if((*itr_triangle2)->edge_vertex_[2] == NULL){T temp=(T)3;}

		AddTriangle((*itr_triangle2)->vertices_[0],(*itr_triangle2)->edge_vertex_[2],(*itr_triangle2)->edge_vertex_[1])->is_old_=false;
		AddTriangle((*itr_triangle2)->edge_vertex_[2],(*itr_triangle2)->edge_vertex_[1],(*itr_triangle2)->edge_vertex_[0])->is_old_=false;
		AddTriangle((*itr_triangle2)->vertices_[1],(*itr_triangle2)->edge_vertex_[0],(*itr_triangle2)->edge_vertex_[2])->is_old_=false;
		AddTriangle((*itr_triangle2)->vertices_[2],(*itr_triangle2)->edge_vertex_[0],(*itr_triangle2)->edge_vertex_[1])->is_old_=false;		

		DynamicTriangle *triangle=*itr_triangle2;
		DelTriangle(triangle);
		itr_triangle2=triangles_.begin();
	}

	std::cout << "triangles are " << triangles_.size() << std::endl;
	std::cout << "Subdivision ended" << std::endl;
}

void DynamicTriangularSurface::SimpleSubdivision()
{
	std::cout << "Subdivision started" << std::endl;
	std::cout << "triangles are " << (int)triangles_.size() << std::endl;

	std::cout << "making new edge vertices" << std::endl;
	for (std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(); itr_triangle != triangles_.end(); itr_triangle++)
	{
		int progress_check = 0;
		if(progress_check > 100)
		{
			std::cout << ".";
			progress_check = 0;
		}
		else progress_check ++;

		(*itr_triangle)->is_old_ = true;

		for(int d = 0; d < 3; d++)	// for 3 edges of this triangle
		{
			if((*itr_triangle)->edge_vertex_[d] == NULL)
			{
				TV3 new_pos, new_nor;
				new_pos = ((*itr_triangle)->vertices_[(d+1)%3]->x_ + (*itr_triangle)->vertices_[(d+2)%3]->x_)*(T)0.5;
				new_nor = ((*itr_triangle)->vertices_[(d+1)%3]->n_ + (*itr_triangle)->vertices_[(d+2)%3]->n_)*(T)0.5;

				(*itr_triangle)->edge_vertex_[d] = AddVertex(new_pos.values_,new_nor.values_);
				(*itr_triangle)->triangles_[d]->edge_vertex_[(*itr_triangle)->triangles_[d]->GetNeighborIndex((*itr_triangle))] = (*itr_triangle)->edge_vertex_[d];
			}
		}
	}

	// generate new triangles
	std::list<DynamicTriangle*> new_triangles_;
	for (std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(); itr_triangle != triangles_.end(); itr_triangle++)
	{
		//TODO: check CCW 
		DynamicTriangle* triangle1 = new DynamicTriangle((*itr_triangle)->vertices_[0],(*itr_triangle)->edge_vertex_[2],(*itr_triangle)->edge_vertex_[1]);
		DynamicTriangle* triangle2 = new DynamicTriangle((*itr_triangle)->edge_vertex_[2],(*itr_triangle)->edge_vertex_[1],(*itr_triangle)->edge_vertex_[0]);
		DynamicTriangle* triangle3 = new DynamicTriangle((*itr_triangle)->vertices_[1],(*itr_triangle)->edge_vertex_[0],(*itr_triangle)->edge_vertex_[2]);
		DynamicTriangle* triangle4 = new DynamicTriangle((*itr_triangle)->vertices_[2],(*itr_triangle)->edge_vertex_[0],(*itr_triangle)->edge_vertex_[1]);
		
		new_triangles_.push_back(triangle1);
		new_triangles_.push_back(triangle2);
		new_triangles_.push_back(triangle3);
		new_triangles_.push_back(triangle4);
	}

	// delete old triangles
	std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(), itr_temp;
	while(itr_triangle != triangles_.end())
	{
		(*itr_triangle)->vertices_[0]->DelTriangle(*itr_triangle);
		(*itr_triangle)->vertices_[1]->DelTriangle(*itr_triangle);
		(*itr_triangle)->vertices_[2]->DelTriangle(*itr_triangle);

		delete (*itr_triangle);

		itr_temp = itr_triangle;

		itr_triangle ++;

		triangles_.remove(*itr_temp);
	}

	// delete old triangle list
//	triangles_.erase(triangles_.begin(), triangles_.end());

	// register new triangles to their vertices
	for (std::list <DynamicTriangle*>::iterator itr_triangle = new_triangles_.begin(); itr_triangle != new_triangles_.end(); itr_triangle++)
	{
		(*itr_triangle)->vertices_[0]->AddTriangle(*itr_triangle);
		(*itr_triangle)->vertices_[1]->AddTriangle(*itr_triangle);
		(*itr_triangle)->vertices_[2]->AddTriangle(*itr_triangle);
	}

	// add new triangles to triangle list
	triangles_.insert(triangles_.end(), new_triangles_.begin(), new_triangles_.end());

	std::cout << "triangles are " << triangles_.size() << std::endl;
	std::cout << "Subdivision ended" << std::endl;
}

void DynamicTriangularSurface::CollapseShortEdges(const T& min_edge_length)
{
	const T sqr_min_edge_length = POW2(min_edge_length);

	std::list<DynamicTriangle*> triangles_to_delete;

	std::cout << "check 1" << std::endl;
	for (std::vector<DynamicVertex*>::iterator itr_vertex = vertices_.begin(); itr_vertex != vertices_.end(); itr_vertex++)
	{
		DynamicVertex *v0 = *itr_vertex, *v1 = NULL;
		DynamicTriangle *triangle0 = NULL, *triangle1 = NULL;

		if(v0->triangles_.empty()) continue;	

		for (std::list <DynamicTriangle*>::iterator itr_triangle = v0->triangles_.begin(); itr_triangle != v0->triangles_.end(); itr_triangle++)
		{
			DynamicTriangle *triangle = *itr_triangle;
			
			assert(triangle != NULL);

			// find short edge v0-v1
			for(int d=0;d<3;++d)
			{
				DynamicVertex *v_temp = triangle->vertices_[d];

				if(v_temp == v0) continue;

				// check edge length
				T sqr_edge_length = 0;
				for(int l=0;l<3;++l) sqr_edge_length += POW2(v0->x_.values_[l] - v_temp->x_.values_[l]);

				if(sqr_edge_length < sqr_min_edge_length)
				{
					v1 = v_temp;

					triangle0 = triangle;	// first triangle that contains v0 and v1

					break;
				}
			}

			if(v1 != NULL) break;
		}

		if(v1 == NULL) continue;	// no short edge

		// find one more triangle that contains v0 and v1 simultaneously. (triange0 was found when v1 was found)
		for (std::list<DynamicTriangle*>::iterator itr = v0->triangles_.begin(); itr != v0->triangles_.end(); itr++)
		{
			if(*itr == triangle0) continue;

			if((*itr)->vertices_[0] == v1 || (*itr)->vertices_[1] == v1 || (*itr)->vertices_[2] == v1)
			{
				assert((*itr)->vertices_[0] == v0 || (*itr)->vertices_[1] == v0 || (*itr)->vertices_[2] == v0);

				triangle1 = *itr;
				break;
			}
		}

		assert(triangle1 != NULL);	//TODO: temporary. valid only for closed surface

		// Remove them triangle list of vertex v0 (triangle list of v1 will be empty)
		v0->triangles_.remove(triangle0);
		v0->triangles_.remove(triangle1);
		v1->triangles_.remove(triangle0);
		v1->triangles_.remove(triangle1);

		// Add them to delete list
		triangles_to_delete.push_back(triangle0);
		triangles_to_delete.push_back(triangle1);

		// let v0 replace v1 in all surrounding triangles
		v0->x_ = (v0->x_ + v1->x_)*(T)0.5;
		v0->n_ = (v0->n_ + v1->n_)*(T)0.5;
		v0->reaction_speed_ = (v0->reaction_speed_ + v1->reaction_speed_)*(T)0.5;

		for (std::list<DynamicTriangle*>::iterator itr = v0->triangles_.begin(); itr != v0->triangles_.end(); itr++)
		{
			if((*itr)->vertices_[0] == v1) (*itr)->vertices_[0] = v0;
			if((*itr)->vertices_[1] == v1) (*itr)->vertices_[1] = v0;
			if((*itr)->vertices_[2] == v1) (*itr)->vertices_[2] = v0;
		}

		for (std::list<DynamicTriangle*>::iterator itr = v1->triangles_.begin(); itr != v1->triangles_.end(); itr++)
		{
			if((*itr)->vertices_[0] == v1) (*itr)->vertices_[0] = v0;
			if((*itr)->vertices_[1] == v1) (*itr)->vertices_[1] = v0;
			if((*itr)->vertices_[2] == v1) (*itr)->vertices_[2] = v0;
		}

		// empty triangle list of v1
		v1->triangles_.clear();
		//TODO: delete v1 later
	}

	std::cout << "check 2" << std::endl;

	// delete triangles to be deleted
	int num_deleted = 0;
	std::list <DynamicTriangle*>::iterator itr_triangle = triangles_to_delete.begin();
	while(itr_triangle != triangles_to_delete.end())
	{
		(*itr_triangle)->vertices_[0]->triangles_.remove(*itr_triangle);
		(*itr_triangle)->vertices_[1]->triangles_.remove(*itr_triangle);
		(*itr_triangle)->vertices_[2]->triangles_.remove(*itr_triangle);

		delete (*itr_triangle);
		triangles_.remove(*itr_triangle);
		num_deleted ++;

		itr_triangle ++;
	}

	std::cout << "Edge collapse num of deleted triangles = " << num_deleted << std::endl;

	std::cout << "check 3" << std::endl;

	// rebuild triangle-triangle connectivity
	for (std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(); itr_triangle != triangles_.end(); itr_triangle++)
		(*itr_triangle)->ChkNeighborConnectivity();

	std::cout << "check 5" << std::endl;
}

void DynamicTriangularSurface::SplitLongEdges(const T& max_edge_length)
{
	const T sqr_max_edge_length = POW2(max_edge_length);

	std::cout << "The number of triangles before edge splitting is " << (int)triangles_.size() << std::endl;

	for (std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(); itr_triangle != triangles_.end(); itr_triangle++)
	{
		for(int d = 0; d < 3; d++)	// for 3 edges of this triangle
		{
			if((*itr_triangle)->edge_vertex_[d] == NULL)
			{
				DynamicVertex *v0 = (*itr_triangle)->vertices_[(d+1)%3], *v1 = (*itr_triangle)->vertices_[(d+2)%3];

				// check edge length
				T sqr_edge_length = 0;
				for(int l=0;l<3;++l) sqr_edge_length += POW2(v0->x_.values_[l] - v1->x_.values_[l]);
	
				// generate edge vertex in the middle of it if edge longer than predefined max edge length
				if(sqr_edge_length > sqr_max_edge_length)
				{
					TV3 new_pos, new_nor;
					new_pos = (v0->x_ + v1->x_)*(T)0.5;
					new_nor = (v0->n_ + v1->n_)*(T)0.5;

					DynamicVertex *new_edge_vertex = AddVertex(new_pos.values_,new_nor.values_);

					new_edge_vertex->reaction_speed_ = (v0->reaction_speed_ + v1->reaction_speed_)*(T)0.5;

					(*itr_triangle)->edge_vertex_[d] = new_edge_vertex;

					DynamicTriangle *neighbor_triangle = (*itr_triangle)->triangles_[d];

					if(neighbor_triangle != 0)
						(*itr_triangle)->triangles_[d]->edge_vertex_[neighbor_triangle->GetNeighborIndex((*itr_triangle))] = new_edge_vertex;
				}
			}
		}
	}

	// generate new triangles
	int num_new_triangles = 0, num_delete_triangles = 0;
	std::list<DynamicTriangle*> new_triangles_;
	for (std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(); itr_triangle != triangles_.end(); itr_triangle++)
	{
		DynamicVertex *v0 = (*itr_triangle)->vertices_[0], *v1 = (*itr_triangle)->vertices_[1], *v2 = (*itr_triangle)->vertices_[2];
		DynamicVertex *e0 = (*itr_triangle)->edge_vertex_[0], *e1 = (*itr_triangle)->edge_vertex_[1], *e2 = (*itr_triangle)->edge_vertex_[2];

		num_delete_triangles ++;

		// 0 edge vertex (no split)
		if(e0 == NULL && e1 == NULL && e2 == NULL)
		{
			// do nothing
			num_delete_triangles --;
		}
		// 1 edge vertex (3 cases) split this triangle into 2 new triangles
		else if(e0 != NULL && e1 == NULL && e2 == NULL)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(e0, v0, v1);
			DynamicTriangle* triangle2 = new DynamicTriangle(v0, e0, v2);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);

			num_new_triangles += 2;
		}
		else if(e0 == NULL && e1 != NULL && e2 == NULL)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(v0, v1, e1);
			DynamicTriangle* triangle2 = new DynamicTriangle(e1, v1, v2);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);

			num_new_triangles += 2;
		}
		else if(e0 == NULL && e1 == NULL && e2 != NULL)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(v0, e2, v2);
			DynamicTriangle* triangle2 = new DynamicTriangle(v2, e2, v1);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);

			num_new_triangles += 2;
		}
		// 2 edge vertex (3 cases)
		else if(e0 == NULL && e1 != NULL && e2 != NULL)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(v0, e2, e1);
			DynamicTriangle* triangle2 = new DynamicTriangle(e2, v1, e1);
			DynamicTriangle* triangle3 = new DynamicTriangle(e1, v1, v2);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);
			new_triangles_.push_back(triangle3);

			num_new_triangles += 3;
		}
		else if(e0 != NULL && e1 == NULL && e2 != NULL)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(v0, e2, v2);
			DynamicTriangle* triangle2 = new DynamicTriangle(e2, v1, e0);
			DynamicTriangle* triangle3 = new DynamicTriangle(e2, e0, v2);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);
			new_triangles_.push_back(triangle3);

			num_new_triangles += 3;
		}
		else if(e0 != NULL && e1 != NULL && e2 == NULL)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(v0, v1, e0);
			DynamicTriangle* triangle2 = new DynamicTriangle(v0, e0, e1);
			DynamicTriangle* triangle3 = new DynamicTriangle(e1, e0, v2);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);
			new_triangles_.push_back(triangle3);

			num_new_triangles += 3;
		}
		else // 3 edge vertexes (1 cases)
		{
			DynamicTriangle* triangle1 = new DynamicTriangle(v0,e2,e1);
			DynamicTriangle* triangle2 = new DynamicTriangle(e0,e1,e2);
			DynamicTriangle* triangle3 = new DynamicTriangle(v1,e0,e2);
			DynamicTriangle* triangle4 = new DynamicTriangle(v2,e1,e0);

			new_triangles_.push_back(triangle1);
			new_triangles_.push_back(triangle2);
			new_triangles_.push_back(triangle3);
			new_triangles_.push_back(triangle4);

			num_new_triangles += 4;
		}
	}

	std::cout<<"New triangles = "<<num_new_triangles<<" To be deleted = "<<num_delete_triangles<<std::endl;

	// delete splited old triangles
	int num_deleted = 0;
	std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(), itr_temp;
	while(itr_triangle != triangles_.end())
	{
		if((*itr_triangle)->edge_vertex_[0] == NULL && (*itr_triangle)->edge_vertex_[1] == NULL && (*itr_triangle)->edge_vertex_[2] == NULL)
		{
			itr_triangle ++;
		}
		else
		{
			(*itr_triangle)->vertices_[0]->DelTriangle(*itr_triangle);
			(*itr_triangle)->vertices_[1]->DelTriangle(*itr_triangle);
			(*itr_triangle)->vertices_[2]->DelTriangle(*itr_triangle);

			itr_temp = itr_triangle;

			itr_triangle ++;

			delete (*itr_temp);
			triangles_.remove(*itr_temp);

			num_deleted ++;
		}
	}

	std::cout << "triangles deleted " << num_deleted << std::endl;

	// register new triangles to their vertices
	for (std::list <DynamicTriangle*>::iterator itr_triangle = new_triangles_.begin(); itr_triangle != new_triangles_.end(); itr_triangle++)
	{
		(*itr_triangle)->vertices_[0]->AddTriangle(*itr_triangle);
		(*itr_triangle)->vertices_[1]->AddTriangle(*itr_triangle);
		(*itr_triangle)->vertices_[2]->AddTriangle(*itr_triangle);
	}

	// add new triangles to triangle list
	triangles_.insert(triangles_.end(), new_triangles_.begin(), new_triangles_.end());

	std::cout << "The number of triangles after edge splitting is " << (int)triangles_.size() << std::endl;

	// rebuild triangle-triangle connectivity
	for (std::list <DynamicTriangle*>::iterator itr_triangle = triangles_.begin(); itr_triangle != triangles_.end(); itr_triangle++)
		(*itr_triangle)->ChkNeighborConnectivity();
}

void DynamicTriangularSurface::ChkTrianglesNeighborConnectivity()
{
	{TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle = *itr_triangle;
		triangle->triangles_[0] = NULL;
		triangle->triangles_[1] = NULL;
		triangle->triangles_[2] = NULL;
	}}

	{TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle = *itr_triangle;
		triangle->ChkNeighborConnectivity();
	}}
}

DynamicVertex* DynamicTriangularSurface::ChkNearestTextureCoordinate(T u,T v)
{
	DynamicVertex *nearest=NULL;
	T distance=(T)100000000;
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle = *itr_triangle;
		DynamicVertex **tri_vertices = triangle->vertices_;
		for(int i=0;i<3;i++)
		{
			T uv[2]={triangle->uv_[i].x_,triangle->uv_[i].y_};
			T dis=(uv[0]-u)*(uv[0]-u)+(uv[1]-v)*(uv[1]-v);
			if(dis<distance)
			{
				distance=dis;
				nearest=tri_vertices[i];
			}
		}
	}
	return nearest;
}

void DynamicTriangularSurface::RemoveNonManifold()
{
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle *triangle=*itr_triangle;		
		if(triangle->vertices_[0]->is_boundary_==true)
		{
			continue;
		}
		if(triangle->vertices_[1]->is_boundary_==true)
		{
			continue;
		}
		if(triangle->vertices_[2]->is_boundary_==true)
		{
			continue;
		}
		if(triangle->triangles_[0]==NULL)
		{
			DelTriangle(triangle);
			itr_triangle=triangles_.begin();
			continue;
		}
		if(triangle->triangles_[1]==NULL)
		{
			DelTriangle(triangle);
			itr_triangle=triangles_.begin();
			continue;
		}
		if(triangle->triangles_[2]==NULL)
		{
			DelTriangle(triangle);
			itr_triangle=triangles_.begin();
			continue;
		}
	}
}

void DynamicTriangularSurface::FillHoles()
{
	int number_of_holes_filled=0;
	TRAVERSE_TRIANGLES
	{
		DynamicVertex *v[4]={NULL,NULL,NULL,NULL};
		for(int d=0;d<3;d++)
		{
			if((*itr_triangle)->triangles_[d]==NULL)
			{
				if((*itr_triangle)->vertices_[(d+1)%3]->is_boundary_==true && (*itr_triangle)->vertices_[(d+2)%3]->is_boundary_==true)
				{
					continue;
				}
				v[0]=(*itr_triangle)->vertices_[(d+1)%3];
				v[1]=(*itr_triangle)->vertices_[(d+2)%3];
				break;
			}
			if(v[0]!=NULL && v[1]!=NULL)
			{
				break;
			}
		}
		if(v[0]==NULL && v[1]==NULL)
		{
			continue;
		}

		std::list<DynamicTriangle*>::iterator itr_triangle;
		for(itr_triangle=v[0]->triangles_.begin();itr_triangle!=v[0]->triangles_.end();itr_triangle++)
		{
			for(int d=0;d<3;d++)
			{
				if((*itr_triangle)->triangles_[d]==NULL)
				{
					if((*itr_triangle)->vertices_[(d+1)%3]->is_boundary_==true && (*itr_triangle)->vertices_[(d+2)%3]->is_boundary_==true)
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[0] && (*itr_triangle)->vertices_[(d+2)%3]==v[1])
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[1] && (*itr_triangle)->vertices_[(d+2)%3]==v[0])
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[0])
					{
						v[2]=(*itr_triangle)->vertices_[(d+2)%3];
						break;
					}
					if((*itr_triangle)->vertices_[(d+2)%3]==v[0])
					{
						v[2]=(*itr_triangle)->vertices_[(d+1)%3];
						break;
					}
					if(v[2]!=NULL)
					{
						break;
					}
				}
			}
			if(v[2]!=NULL)
			{
				break;
			}
		}
		for(itr_triangle=v[1]->triangles_.begin();itr_triangle!=v[1]->triangles_.end();itr_triangle++)
		{
			for(int d=0;d<3;d++)
			{
				if((*itr_triangle)->triangles_[d]==NULL)
				{
					if((*itr_triangle)->vertices_[(d+1)%3]->is_boundary_==true && (*itr_triangle)->vertices_[(d+2)%3]->is_boundary_==true)
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[1] || (*itr_triangle)->vertices_[(d+2)%3]==v[0])
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[0] || (*itr_triangle)->vertices_[(d+2)%3]==v[1])
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[1] || (*itr_triangle)->vertices_[(d+2)%3]==v[2])
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[2] || (*itr_triangle)->vertices_[(d+2)%3]==v[1])
					{
						continue;
					}
					if((*itr_triangle)->vertices_[(d+1)%3]==v[1])
					{
						v[3]=(*itr_triangle)->vertices_[(d+2)%3];
						break;
					}
					else
					{
						v[3]=(*itr_triangle)->vertices_[(d+1)%3];
						break;
					}
					if(v[3]!=NULL)
					{
						break;
					}
				}
			}
		}
		if(v[2]!=NULL)
		{
			AddTriangle(v[0],v[1],v[2])->CorrectCCW();
		}
		if(v[3]!=NULL)
		{
			AddTriangle(v[0],v[2],v[3])->CorrectCCW();
		}
		if(v[2]!=NULL || v[3]!=NULL)
		{
			number_of_holes_filled++;			
			ChkTrianglesNeighborConnectivity();
			std::cout<<"."<<std::flush;
		}
	}
	std::cout<<"Total holes filled "<<number_of_holes_filled<<std::endl;
}

void DynamicTriangularSurface::ChkBoundary()
{
	std::cout<<"# Boundary Checking Started"<<std::endl;
	TRAVERSE_TRIANGLES
	{
		if((*itr_triangle)->triangles_[0]==NULL)
		{
			AddEdge((*itr_triangle)->vertices_[0],(*itr_triangle)->vertices_[1]);
		}
		if((*itr_triangle)->triangles_[1]==NULL)
		{
			AddEdge((*itr_triangle)->vertices_[1], (*itr_triangle)->vertices_[2]);
		}
		if((*itr_triangle)->triangles_[2]==NULL)
		{
			AddEdge((*itr_triangle)->vertices_[2], (*itr_triangle)->vertices_[0]);
		}
	}
	std::cout<<"# Finished"<<std::endl;
	std::cout<<"# Number of Edges = "<<edges_.size()<<std::endl;

	std::list <DynamicEdge*>::iterator itr_edge = edges_.begin();
	while(itr_edge!=edges_.end())
	{
		std::cout<<"# Number of Edges = "<<edges_.size()<<std::endl;
		DynamicHole *hole=new DynamicHole;
		holes_.push_back(hole);
		hole->AddEdge((*itr_edge));
		DynamicVertex *v0=(*itr_edge)->vertices_[0];		
		DynamicTriangularSurface::edges_.erase(itr_edge);
		std::cout<<"# Number of Edges = "<<edges_.size()<<std::endl;
		std::list <DynamicEdge*>::iterator itr_edge = edges_.begin();
		while(itr_edge!=edges_.end())
		{
			if((*itr_edge)->vertices_[0]==v0 || (*itr_edge)->vertices_[1]==v0)
			{
				hole->AddEdge((*itr_edge));
				if(v0==(*itr_edge)->vertices_[0])
				{
					v0=(*itr_edge)->vertices_[1];
				}
				else
				{
					v0=(*itr_edge)->vertices_[0];
				}
				edges_.erase(itr_edge);
				itr_edge=edges_.begin();
			}
			else
			{
				itr_edge++;
			}
		}
		itr_edge=edges_.begin();
	}
	std::cout<<"# Finished"<<std::endl;
}

void DynamicTriangularSurface::MinMaxPosition(TV3& min, TV3& max)
{
	T min_compare_x = vertices_[0]->x_.values_[0];
	T min_compare_y = vertices_[0]->x_.values_[1];
	T min_compare_z = vertices_[0]->x_.values_[2];

	T max_compare_x = vertices_[0]->x_.values_[0];
	T max_compare_y = vertices_[0]->x_.values_[1];
	T max_compare_z = vertices_[0]->x_.values_[2];

	for(unsigned int i=1; i<vertices_.size(); i++)
	{
		TV3 position(vertices_[i]->x_);
		
		min_compare_x = MIN2(position.x_, min_compare_x);
		min_compare_y = MIN2(position.y_, min_compare_y);
		min_compare_z = MIN2(position.z_, min_compare_z);

		max_compare_x = MAX2(position.x_, max_compare_x);
		max_compare_y = MAX2(position.y_, max_compare_y);
		max_compare_z = MAX2(position.z_, max_compare_z);
	}

	min = TV3(min_compare_x, min_compare_y, min_compare_z);
	max = TV3(max_compare_x, max_compare_y, max_compare_z);
}


TV3 &DynamicTriangle::GetTVNormal()
{
	return n_;
}
TV3 DynamicTriangle::GetTVCenter()
{
	TV3 p0(vertices_[0]->x_);
	TV3 p1(vertices_[1]->x_);
	TV3 p2(vertices_[2]->x_);
	return (p0+p1+p2)/(T)3;
}

void DynamicTriangle::GetCenter(T* center)
{
	center[0] = vertices_[0]->x_.values_[0] + vertices_[1]->x_.values_[0] + vertices_[2]->x_.values_[0];
	center[1] = vertices_[0]->x_.values_[1] + vertices_[1]->x_.values_[1] + vertices_[2]->x_.values_[1];
	center[2] = vertices_[0]->x_.values_[2] + vertices_[1]->x_.values_[2] + vertices_[2]->x_.values_[2];

	center[0] /= (T)3; center[1] /= (T)3; center[2] /= (T)3;
}

T DynamicTriangle::SignedDistance(const TV3& location) const
{
	TV3 closest_pt = ClosestPoint(location);
	TV3 normal(n_);

	T signed_distance = (location-closest_pt).getMagnitude();

	if(dotProduct(normal, location-closest_pt) < (T)0) signed_distance = -signed_distance;

	return signed_distance;
}

TV3 DynamicTriangle::ClosestPoint(const TV3& location) const
{
	TV3 v0(vertices_[0]->x_);
	TV3 v1(vertices_[1]->x_);
	TV3 v2(vertices_[2]->x_);
	return ClosestPointFromTriangle(location, v0, v1, v2);
}

TV3 DynamicTriangle::ClosestPoint(const TV3& location, bool& on_triangle) const
{
	TV3 v0(vertices_[0]->x_);
	TV3 v1(vertices_[1]->x_);
	TV3 v2(vertices_[2]->x_);
	return ClosestPointFromTriangle(location, v0, v1, v2, on_triangle);
}

TV3 DynamicTriangle::BarycentricCoordinates(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3)
{
	TV3 u = x2-x1, v = x3-x1,w = location-x1;

	T u_dot_u = dotProduct(u,u), v_dot_v = dotProduct(v,v), u_dot_v = dotProduct(u,v),
		u_dot_w = dotProduct(u,w), v_dot_w = dotProduct(v,w);

	T denominator = u_dot_u*v_dot_v - POW2(u_dot_v), one_over_denominator;

	if(abs(denominator) > (T)1e-16)
	{
		one_over_denominator = 1/denominator;
	}
	else
	{
		one_over_denominator=(T)1e16;
	}

	T a = (v_dot_v*u_dot_w - u_dot_v*v_dot_w)*one_over_denominator, b = (u_dot_u*v_dot_w-u_dot_v*u_dot_w)*one_over_denominator;

	return TV3(1-a-b,a,b);
}

TV3 DynamicTriangle::ClosestPointFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3)
{
	TV3 nor = crossProduct(x2-x1, x3-x1); nor.normalize(); // normal unit of plane
	TV3 p = location-nor*dotProduct(nor,location-x1);

	TV3 weights = BarycentricCoordinates(p, x1, x2, x3);// can be reused for velocity calculation

	// project closest point to the triangle if it's not already inside it

	if(weights.y_<(T)0 && weights.z_<(T)0)
	{
		return x1; //return x1 
	}
	else if(weights.x_<(T)0 && weights.z_<(T)0)
	{
		return x2; //return x2 
	}
	else if(weights.x_<(T)0 && weights.y_<(T)0)
	{
		return x3; //return x3 
	}

	if(weights.x_ < (T)0) // Closest point is on edge x2--x3
	{
		return ClosestPointFromLine(p, x2, x3);
	}
	else if(weights.y_ < (T)0) // Closest point is on edge x1--x3
	{
		return ClosestPointFromLine(p, x1, x3);
	}
	else if(weights.z_ < (T)0) // Closest point is on edge x1--x2
	{
		return ClosestPointFromLine(p, x1, x2);
	}

	return p;
}


TV3 DynamicTriangle::ClosestPointFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3, bool& on_triangle)
{
	TV3 nor = crossProduct(x2-x1, x3-x1); nor.normalize(); // normal unit of plane
	TV3 p = location-nor*dotProduct(nor,location-x1);

	TV3 weights = BarycentricCoordinates(p, x1, x2, x3);// can be reused for velocity calculation

	// project closest point to the triangle if it's not already inside it
	on_triangle = false;
	if(weights.y_<(T)0 && weights.z_<(T)0)
	{
		return x1; //return x1 
	}
	else if(weights.x_<(T)0 && weights.z_<(T)0)
	{
		return x2; //return x2 
	}
	else if(weights.x_<(T)0 && weights.y_<(T)0)
	{
		return x3; //return x3 
	}

	if(weights.x_ < (T)0) // Closest point is on edge x2--x3
	{
		return ClosestPointFromLine(p, x2, x3);
	}
	else if(weights.y_ < (T)0) // Closest point is on edge x1--x3
	{
		return ClosestPointFromLine(p, x1, x3);
	}
	else if(weights.z_ < (T)0) // Closest point is on edge x1--x2
	{
		return ClosestPointFromLine(p, x1, x2);
	}

	on_triangle = true;
	return p;
}

T DynamicTriangle::SignedDistanceFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3)
{
	TV3 closest_pt = ClosestPointFromTriangle(location, x1, x2, x3);
	TV3 normal = crossProduct(x2-x1, x3-x1); normal.normalize(); // normal unit of plane
	T min_dist = (location-closest_pt).getMagnitude();

	if(dotProduct(normal, location-closest_pt) < (T)0) return -min_dist;
	else return min_dist;		
}

T DynamicTriangle::SignedDistanceFromTriangle(const TV3& location, const TV3& x1, const TV3& x2, const TV3& x3, bool& on_triangle)
{
	TV3 closest_pt = ClosestPointFromTriangle(location, x1, x2, x3, on_triangle);
	TV3 normal = crossProduct(x2-x1, x3-x1); normal.normalize(); // normal unit of plane
	T min_dist = (location - closest_pt).getMagnitude();

	if(dotProduct(normal, location-closest_pt) < (T)0) return -min_dist;
	else return min_dist;		
}

TV3 DynamicTriangle::ClosestPointFromLine(const TV3& location, const TV3& x1, const TV3& x2)
{
	T p = dotProduct(location-x1, x2-x1)/dotProduct(x2-x1, x2-x1);

	if(p < (T)0)
	{
		return x1;
	}
	else if(p > (T)1)
	{
		return x2;
	}

	return x1+(x2-x1)*p;
}

bool DynamicTriangle::IntersectTriangle(const TV3& i0, const TV3& i1, const DynamicTriangle& triangle, TV3& intersection_position)
{
	TV3 p0(triangle.vertices_[0]->x_);
	TV3 p1(triangle.vertices_[1]->x_);
	TV3 p2(triangle.vertices_[2]->x_);

	TV3 normal(triangle.n_); normal.normalize();

	TV3 q0 = i0 - p0;
	TV3 q1 = i1 - p0;

	T d0 = dotProduct(normal, q0);
	T d1 = dotProduct(normal, q1);

	TV3 closest_position = i0 - normal*d0;

	TV3 line_direction = i1 - i0;
	line_direction.normalize();

	intersection_position = i0 + line_direction*dotProduct(line_direction, closest_position - i0);

	TV3 weights = BarycentricCoordinates(intersection_position, p0, p1, p2);

	if(weights.x_ >= (T)0 && weights.y_ >= (T)0 && weights.z_ >= (T)0 && d0*d1 < (T)0) return true;

	return false;
}

bool DynamicTriangle::RayThruTriangle(const DynamicTriangle& triangle, const TV3& R1, const TV3& R2,  TV3& intersection_position_out, const T& offset)
{
	TV3 P1(triangle.vertices_[0]->x_);
	TV3 P2(triangle.vertices_[1]->x_);
	TV3 P3(triangle.vertices_[2]->x_);

	// Find Triangle Normal
	TV3 normal(triangle.n_);
	normal.normalize();

	TV3 intersect_position;
	// Find distance from LP1 and LP2 to the plane defined by the triangle
	T dist1 = dotProduct((R1-P1), normal);
	T dist2 = dotProduct((R2-P1), normal);

	if(dist1 * dist2 >= (T)0 || dist1 == dist2) return false;

	// Find point on the line that intersects with the plane
	intersect_position = R1 + (R2-R1) * ( -dist1/(dist2-dist1) );

	// Find if the interesection point lies inside the triangle by testing it against all edges
	TV3 vTest;

	vTest = crossProduct(normal, P2-P1); vTest.normalize();
	if(dotProduct(vTest, intersect_position - P1) < 0) return false;

	vTest = crossProduct(normal, P3-P2); vTest.normalize();
	if(dotProduct(vTest, intersect_position - P2) < 0) return false;

	vTest = crossProduct(normal, P1-P3); vTest.normalize();
	if(dotProduct(vTest, intersect_position - P3) < 0) return false;

	intersection_position_out = intersect_position;

	return true;
}

TV3 DynamicTriangularSurface::Center()
{
	T num = 0;
	TV3 center = TV3();
	TRAVERSE_VERTICES
	{
		num += (T)1;		
		center += (*itr_vertex)->x_;
	}

	if(num > 1e-08) center = center/num;
	return center;
}

void DynamicTriangularSurface::InertiaTensor(MATRIX_3X3& inertia_mat, const T& mass)
{
	TV3 center = Center();
	T dm = mass / (T)vertices_.size();

	// compute inertia Tensor
	T i_xx(0), i_yy(0), i_zz(0), i_xy(0), i_xz(0), i_yz(0);
	TRAVERSE_VERTICES
	{
		TV3 p = (*itr_vertex)->x_;
		p = p - center;
		
		i_xx += dm*(POW2(p.y_)+POW2(p.z_));
		i_yy += dm*(POW2(p.z_)+POW2(p.x_));
		i_zz += dm*(POW2(p.x_)+POW2(p.y_));

		i_xy += dm*(p.x_*p.y_);
		i_xz += dm*(p.x_*p.z_);
		i_yz += dm*(p.y_*p.z_);
	}

	if(i_xx < 1e-08) i_xx = (T)0;
	if(i_yy < 1e-08) i_yy = (T)0;
	if(i_zz < 1e-08) i_zz = (T)0;
	if(i_xy < 1e-08) i_xy = (T)0;
	if(i_xz < 1e-08) i_xz = (T)0;
	if(i_yz < 1e-08) i_yz = (T)0;

	inertia_mat = MATRIX_3X3(i_xx, -i_xy, -i_xz, -i_xy, i_yy, -i_yz, -i_xz, -i_yz, i_zz);
}

void DynamicTriangularSurface::SamplingPointsOnSurface(std::vector<TV3>& arr, const T dw)
{
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle& triangle = *(*itr_triangle);

		TV3& v0 = triangle.vertices_[0]->x_;
		TV3& v1 = triangle.vertices_[1]->x_;
		TV3& v2 = triangle.vertices_[2]->x_;

		T d0 = dw / (v1 - v0).getMagnitude();

		TV3 dir0 = v1-v0;

		for(T c0=(T)0; c0<=(T)1; c0+=d0)
		{
			TV3 q = v0 + dir0*c0;
			TV3 dir1 = v2-q;

			T d1 = dw / (v2 - q).getMagnitude();

			for(T c1=(T)0; c1<=(T)1; c1+=d1)
			{
				TV3 p = q + dir1*c1;			
				arr.push_back(p);
			}
		}
	}
}

void DynamicTriangularSurface::SamplingPointsOnSurface(std::vector<TV3>& vert_arr, std::vector<TV3>& nor_arr, const T dw)
{
	TRAVERSE_TRIANGLES
	{
		DynamicTriangle& triangle = *(*itr_triangle);

		TV3& v0 = triangle.vertices_[0]->x_;
		TV3& v1 = triangle.vertices_[1]->x_;
		TV3& v2 = triangle.vertices_[2]->x_;

		TV3 nor = triangle.n_.normalized();

		T h = dw*(T)4;

		T d0 = dw / (v1 - v0).getMagnitude();

		TV3 dir0 = v1 - v0;

		for (T c0 = (T)0; c0 <= (T)1; c0 += d0)
		{
			TV3 q = v0 + dir0*c0;
			TV3 dir1 = v2 - q;

			T d1 = dw / (v2 - q).getMagnitude();

			for (T c1 = (T)0; c1 <= (T)1; c1 += d1)
			{
				TV3 p = q + dir1*c1;
				TV3 n = TV3();

				n = nor;

				vert_arr.push_back(p);				
				nor_arr.push_back(n);
			}
		}
	}
}