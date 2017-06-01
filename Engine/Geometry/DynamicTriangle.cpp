#include <assert.h>

#include "CONVENTIONAL_MACROS.h"
#include "DynamicTriangle.h"
#include "DynamicEdge.h"
#include "DynamicVertex.h"

DynamicTriangle::DynamicTriangle()
{
	vertices_[0] = NULL; vertices_[1] = NULL; vertices_[2] = NULL;
	edge_vertex_[0] = NULL; edge_vertex_[1] = NULL; edge_vertex_[2] = NULL;

	is_old_ = false; wrong_ = false; is_surface_ = true;
}

DynamicTriangle::DynamicTriangle(DynamicVertex* v0,DynamicVertex* v1,DynamicVertex* v2)
{
	vertices_[0]=v0;vertices_[1]=v1;vertices_[2]=v2;
	edge_vertex_[0]=NULL;edge_vertex_[1]=NULL;edge_vertex_[2]=NULL;

	is_old_ = false; wrong_ = false; is_surface_ = true;
}

DynamicTriangle::~DynamicTriangle()
{

}

void DynamicTriangle::DelTriangle(DynamicTriangle* triangle)
{
	for(int d=0;d<3;d++)
	{
		if(triangles_[d]==triangle)
		{
			triangles_[d]=NULL;
			return;
		}
	}
}

T* DynamicTriangle::GetNormal()
{
	return n_.values_;
}

T DynamicTriangle::GetOppositeEdgeLength(DynamicVertex* v)
{
	int d;
	for(d=0;d<3;d++)
	{
		if(DynamicTriangle::vertices_[d]==v)
		{
			break;
		}
	}

	TV3 length = vertices_[(d+1)%3]->x_ - vertices_[(d+2)%3]->x_;

//	ARRAY_VECTOR3::set<T>(length,vertices_[(d+1)%3]->GetPosition());
//	ARRAY_VECTOR3::sub<T>(length,vertices_[(d+2)%3]->GetPosition(),length);

	return length.getMagnitude();
}

void DynamicTriangle::AddLocalCurvatureNormal()
{
	for(int d=0;d<3;d++)
	{
		TV3 l0,l1;

		l0 = vertices_[(d+1)%3]->x_ - vertices_[d]->x_;
		l1 = vertices_[(d+2)%3]->x_ - vertices_[d]->x_;

//		ARRAY_VECTOR3::set<T>(l0,vertices_[(d+1)%3]->GetPosition());
//		ARRAY_VECTOR3::sub<T>(l0,vertices_[d]->GetPosition(),l0);
//		ARRAY_VECTOR3::set<T>(l1,vertices_[(d+2)%3]->GetPosition());
//		ARRAY_VECTOR3::sub<T>(l1,vertices_[d]->GetPosition(),l1);

		TV3 length;

		length = vertices_[(d+1)%3]->x_ - vertices_[(d+2)%3]->x_;
//		ARRAY_VECTOR3::set<T>(length,vertices_[(d+1)%3]->GetPosition());
//		ARRAY_VECTOR3::sub<T>(length,vertices_[(d+2)%3]->GetPosition(),length);
		T weight=length.getMagnitude();
//		ARRAY_VECTOR3::normalize(l0);
//		ARRAY_VECTOR3::normalize(l1);

		l0 *= weight;
		l1 *= weight;

		vertices_[d]->curvature_normal_ += l0 + l1;

//		ARRAY_VECTOR3::mul<T>(l0,weight);
//		ARRAY_VECTOR3::mul<T>(l1,weight);
//		ARRAY_VECTOR3::add<T>(vertices_[d]->curvature_normal_.values_,l0,vertices_[d]->curvature_normal_.values_);
//		ARRAY_VECTOR3::add<T>(vertices_[d]->curvature_normal_.values_,l1,vertices_[d]->curvature_normal_.values_);
		vertices_[d]->normalizer_+=weight;
	}
}

int DynamicTriangle::GetNeighborIndex(DynamicTriangle* triangle)
{
	for(int d=0;d<3;d++)
	{
		if(triangles_[d]==triangle)
		{
			return d;
		}
	}
	std::cout<<"TRIANGLE::GetNeighborIndex error!"<<std::endl;
	return 0;
}

int DynamicTriangle::GetVertexIndex(DynamicVertex* v)
{
	for(int d=0;d<3;d++)
	{
		if(vertices_[d]==v)
		{
			return d;
		}
	}
	std::cout<<"TRIANGLE::getvertexIndex error!"<<std::endl;
	return -1;
}

bool DynamicTriangle::IsInside(DynamicVertex* v)
{
	if(vertices_[0]==v)
	{
		return true;
	}
	if(vertices_[1]==v)
	{
		return true;
	}
	if(vertices_[2]==v)
	{
		return true;
	}
	return false;
}

void DynamicTriangle::Flip(DynamicVertex* v0,DynamicVertex* v1)
{
	/*
	T center[2];
	ARRAY_VECTOR3::set(center, v0->x_);
	ARRAY_VECTOR3::set(center, v1->x_);
	ARRAY_VECTOR3::div(center, 2.0f);

	VERTEX*  v_after = new VERTEX(center, TRIANGLE::n);
	v0->Replace(v_after);
	v1->Replace(v_after);

	TRIANGLE*  n = NULL;
	for(int i = 0; i < 3; i ++)
	{
		if(triangles[i] == NULL)
			return;

		if(triangles[i]->IsInside(v0) && triangles[i]->IsInside(v1))
		{
			n = triangles[i];
		}
	}

	if(n != NULL)
	{
		TRIANGLE*  n0, n1;		// ?̿????? ???? ?ƴ? ?ٸ? ?̿????? ?????͸? ã?´?.
		for(int j = 0; j < 3; j ++)
		{
			if(n->triangles_[j] == this)
				break;
		}
		n0 = n->triangles_[(j+1)%3];
		n1 = n->triangles_[(j+2)%3];

		if(n0 != NULL)
			for(int k = 0; k < 3; k ++)
			{
				if(n0->triangles_[k] == triangles[i])
					n0->triangles_[k] = n1;
			}
		if(n1 != NULL)
			for(int k = 0; k < 3; k ++)
			{
				if(n1->triangles_[k] == triangles[i])
					n1->triangles_[k] = n0;
			}
	}
	*/
}

void DynamicTriangle::Flip()
{
	/*
	T center[3];
	ARRAY_VECTOR3::set(center, TRIANGLE::vertices[0]->x_);
	ARRAY_VECTOR3::add(center, TRIANGLE::vertices[1]->x_, center);
	ARRAY_VECTOR3::add(center, TRIANGLE::vertices[2]->x_, center);
	ARRAY_VECTOR3::div(center, 3.0f);

	VERTEX*  v_after = new VERTEX(center, TRIANGLE::n);
	for(int i = 0; i < 3; i ++)
	{
		VERTEX*  vertex = TRIANGLE::vertices[i];
		TRIANGLE::vertices[i]->Replace(v_after);		
	}
	
	// ?̿????? ?̿????? ?̿? �???? ??????Ʈ ???ش?.
	for(int i = 0; i < 3; i ++)
	{
		if(triangles[i] == NULL)
			continue;

		TRIANGLE*  n0, n1;		// ?̿????? ???? ?ƴ? ?ٸ? ?̿????? ?????͸? ã?´?.
		for(int j = 0; j < 3; j ++)
		{
			if(triangles[i]->triangles_[j] == this)
				break;
		}
		n0 = triangles[i]->triangles_[(j+1)%3];
		n1 = triangles[i]->triangles_[(j+2)%3];

		if(n0 != NULL)
			for(int k = 0; k < 3; k ++)
			{
				if(n0->triangles_[k] == triangles[i])
					n0->triangles_[k] = n1;
			}
		if(n1 != NULL)
			for(int k = 0; k < 3; k ++)
			{
				if(n1->triangles_[k] == triangles[i])
					n1->triangles_[k] = n0;
			}
	}

	for(i = 0; i < 3; i ++)
	{
		TRIANGLE::vertices[i]->triangles_.remove(this);
		v_after->triangles_.remove(this);
		if(triangles[0] != NULL)
		{
			v_after->triangles_.remove(triangles[0]);
			TRIANGLE::vertices[i]->triangles_.remove(triangles[0]);
		}
		if(triangles[1] != NULL)
		{
			v_after->triangles_.remove(triangles[1]);
			TRIANGLE::vertices[i]->triangles_.remove(triangles[1]);
		}
		if(triangles[2] != NULL)
		{
			v_after->triangles_.remove(triangles[2]);
			TRIANGLE::vertices[i]->triangles_.remove(triangles[2]);
		}
	}
	*/
}

void DynamicTriangle::SetNormal(T *n_input)
{
	n_.x_ = n_input[0];
	n_.y_ = n_input[1];
	n_.z_ = n_input[2];
}

int DynamicTriangle::CountFeatureVertex()
{
	int number=0;
	for(int i=0;i<3;i++)
	{
		if(vertices_[i]->feature_==true)
		{
			number++;
		}
	}
	return number;
}

void DynamicTriangle::ChkNeighborConnectivity()
{
	for(int d = 0; d < 3; d++)
	{
		DynamicVertex *v0 = vertices_[(d+1)%3], *v1 = vertices_[(d+2)%3];

		std::list <DynamicTriangle*>::iterator itr;

		// find a triangle that contains v0 and v1 simultaneously.
		for(itr = v0->triangles_.begin(); itr != v0->triangles_.end(); itr++)
		{
			if(*itr == this) continue;

			if((*itr)->vertices_[0] == v1 || (*itr)->vertices_[1] == v1 || (*itr)->vertices_[2] == v1)
			{
				triangles_[d] = *itr;
				break;
			}
		}

		assert(triangles_[d] != NULL);	//TODO: temporary. valid when surface is closed
	}
}

void DynamicTriangle::CorrectCCW()
{
	TV3 l0, l1;
	l0 = DynamicTriangle::vertices_[2]->x_ - DynamicTriangle::vertices_[0]->x_;
	l1 = DynamicTriangle::vertices_[1]->x_ - DynamicTriangle::vertices_[0]->x_;
//	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices_[2]->GetPosition(), TRIANGLE::vertices_[0]->GetPosition(), l0);
//	ARRAY_VECTOR3::sub<T>(TRIANGLE::vertices_[1]->GetPosition(), TRIANGLE::vertices_[0]->GetPosition(), l1);
	
	TV3 c = crossProduct(l1, l0);
//	ARRAY_VECTOR3::cross<T>(l1, l0, c);					// determine face normal
	
	DynamicTriangle::area_ = c.getMagnitude()/(T)2;		// calculate triangle area

	TV3 n = DynamicTriangle::vertices_[0]->n_ + DynamicTriangle::vertices_[1]->n_ + DynamicTriangle::vertices_[2]->GetNormal();

	if((c.x_*n.x_ + c.y_*n.y_ + c.z_*n.z_) < (T)0)
	{
		DynamicVertex* temp = DynamicTriangle::vertices_[1];
		DynamicTriangle::vertices_[1] = DynamicTriangle::vertices_[2];
		DynamicTriangle::vertices_[2] = temp;		
	}
}

/*
void DynamicTriangle::Draw()
{
	glBegin(GL_TRIANGLES);
	{
		#ifdef USE_FLOAT_T		
			glNormal3fv(vertices_[0]->n_.values_);
			glTexCoord2f(uv_[0].x_, uv_[0].y_); 
			glVertex3fv(vertices_[0]->x_.values_);
			glNormal3fv(vertices_[1]->n_.values_);
			glTexCoord2f(uv_[1].x_, uv_[1].y_); 
			glVertex3fv(vertices_[1]->x_.values_);
			glNormal3fv(vertices_[2]->n_.values_);
			glTexCoord2f(uv_[2].x_, uv_[2].y_); 
			glVertex3fv(vertices_[2]->x_.values_);
		#else
			glNormal3dv(vertices_[0]->GetNormal());
			glTexCoord2d(uv_[0].x_, uv_[0].y_); 
			glVertex3dv(vertices_[0]->GetPosition());
			glNormal3dv(vertices_[1]->GetNormal());
			glTexCoord2d(uv_[1].x_, uv_[1].y_); 
			glVertex3dv(vertices_[1]->GetPosition());
			glNormal3dv(vertices_[2]->GetNormal());
			glTexCoord2d(uv_[2].x_, uv_[2].y_); 
			glVertex3dv(vertices_[2]->GetPosition());
		#endif		
	}
	glEnd();
}
*/

/*
void DynamicTriangle::DrawBack()
{
	glBegin(GL_TRIANGLES);
	{
		#ifdef USE_FLOAT_T				
			glNormal3f(-vertices_[0]->GetNormal()[0],-vertices_[0]->GetNormal()[1],-vertices_[0]->GetNormal()[2]);
			glTexCoord2f(uv_[0].x_, uv_[0].y_); 
			glVertex3fv((float*)vertices_[0]->GetPosition());
			glNormal3f(-vertices_[2]->GetNormal()[0],-vertices_[2]->GetNormal()[1],-vertices_[2]->GetNormal()[2]);
			glTexCoord2f(uv_[2].x_, uv_[2].y_); 
			glVertex3fv((float*)vertices_[2]->GetPosition());	
			glNormal3f(-vertices_[1]->GetNormal()[0],-vertices_[1]->GetNormal()[1],-vertices_[1]->GetNormal()[2]);
			glTexCoord2f(uv_[1].x_, uv_[1].y_); 
			glVertex3fv((float*)vertices_[1]->GetPosition());
		#else
			glNormal3d(-vertices_[0]->GetNormal()[0],-vertices_[0]->GetNormal()[1],-vertices_[0]->GetNormal()[2]);
			glTexCoord2d(uv_[0].x_, uv_[0].y_); 
			glVertex3dv((double*)vertices_[0]->GetPosition());
			glNormal3d(-vertices_[2]->GetNormal()[0],-vertices_[2]->GetNormal()[1],-vertices_[2]->GetNormal()[2]);
			glTexCoord2d(uv_[2].x_, uv_[2].y_); 
			glVertex3dv((double*)vertices_[2]->GetPosition());	
			glNormal3d(-vertices_[1]->GetNormal()[0],-vertices_[1]->GetNormal()[1],-vertices_[1]->GetNormal()[2]);
			glTexCoord2d(uv_[1].x_, uv_[1].y_); 
			glVertex3dv((double*)vertices_[1]->GetPosition());
		#endif
	}
	glEnd();	
}
*/

/*
void DynamicTriangle::DrawEdges()
{
//	if(TRIANGLE::triangles[0] == NULL)
	{
		glBegin(GL_LINES);
		#ifdef USE_FLOAT_T		
			glVertex3fv(DynamicTriangle::vertices_[1]->GetPosition());
			glVertex3fv(DynamicTriangle::vertices_[2]->GetPosition());
		#else
			glVertex3dv(DynamicTriangle::vertices_[1]->GetPosition());
			glVertex3dv(DynamicTriangle::vertices_[2]->GetPosition());
		#endif
		
		glEnd();
	}

//	if(TRIANGLE::triangles[1] == NULL)
	{
		glBegin(GL_LINES);
		#ifdef USE_FLOAT_T		
			glVertex3fv(DynamicTriangle::vertices_[2]->GetPosition());
			glVertex3fv(DynamicTriangle::vertices_[0]->GetPosition());
		#else
			glVertex3dv(DynamicTriangle::vertices_[2]->GetPosition());
			glVertex3dv(DynamicTriangle::vertices_[0]->GetPosition());
		#endif		
		glEnd();
	}

//	if(TRIANGLE::triangles[2] == NULL)
	{
		glBegin(GL_LINES);
		#ifdef USE_FLOAT_T		
			glVertex3fv(DynamicTriangle::vertices_[0]->GetPosition());
			glVertex3fv(DynamicTriangle::vertices_[1]->GetPosition());
		#else
			glVertex3dv(DynamicTriangle::vertices_[0]->GetPosition());
			glVertex3dv(DynamicTriangle::vertices_[1]->GetPosition());
		#endif		
		glEnd();
	}	
}
*/

/*
void DynamicTriangle::DrawCenter()
{
	TV3 x_center = (vertices_[0]->x_ + vertices_[1]->x_ + vertices_[2]->x_)/(T)3;

	glBegin(GL_POINTS);
		glVertex3fv((float*)x_center.values_);
	glEnd();
}
*/

DynamicVertex*  DynamicTriangle::FindAnotherVertex(DynamicVertex*  v0, DynamicVertex*  v1)
{
	for(int d = 0; d < 3; d ++)
	{
		if(vertices_[d] != v0 && vertices_[d] != v1)
		{
			return vertices_[d];
		}
	}
	
	return NULL;
}

/*
void DynamicTriangle::DrawNeighborConnectivity()
{
	TV3 center = (vertices_[0]->x_ + vertices_[1]->x_ + vertices_[2]->x_)/(T)3;
	
	for(int i = 0; i < 3; i ++)
	{
		if(DynamicTriangle::triangles_[i] == NULL)
		{
			std::cout << "Null triangle connectivity" << std::endl;
			continue;
		}

		TV3 neighbor = (triangles_[i]->vertices_[0]->x_+triangles_[i]->vertices_[1]->x_+triangles_[i]->vertices_[2]->x_)/(T)3;
		
		neighbor = (center + neighbor)*(T)0.5;

		glBegin(GL_LINES);
			glVertex3fv((float*)center.values_);
			glVertex3fv((float*)neighbor.values_);
		glEnd();

	}
}
*/

/*
void DynamicTriangle::DrawNormal()
{
	TV3 p = (vertices_[0]->x_ + vertices_[1]->x_ + vertices_[2]->x_)/(T)3;

	glPushMatrix();		
	glTranslatef(p.x_, p.y_, p.z_);
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3fv((float*)DynamicTriangle::GetNormal());	
	glEnd();
	glPopMatrix();
}
*/

void DynamicTriangle::DetermineNormal()
{
	TV3 l0 = vertices_[2]->x_ - vertices_[0]->x_;
	TV3 l1 = vertices_[1]->x_ - vertices_[0]->x_;
	n_ = crossProduct(l1, l0);

	DynamicTriangle::area_ = n_.getMagnitude()*(T)0.5;
	
	n_.normalize();
//	ARRAY_VECTOR3::normalize<T>(n_.values_);	
}

void DynamicTriangle::DetermineAverageNormal()
{
//	T ave_normal[3] = {(T)0, (T)0, (T)0};

	TV3 center = GetTVCenter();
//	GetCenter(center);

	TV3 normal_weighted = TV3((T)0, (T)0, (T)0);

	for(int i=0; i < 3; ++i)
	{
//		T dist = ARRAY_VECTOR3::dist(center, vertices_[i]->x_.values_);
		T dist = (center-vertices_[i]->x_).getMagnitude();

		normal_weighted += vertices_[i]->n_ * (T)1/(dist+(T)FLT_EPSILON);		
	}

	normal_weighted.normalize();

	n_ = (n_+normal_weighted)*(T)0.5;

	//ARRAY_VECTOR3::add<T>(TRIANGLE::vertices_[0]->GetNormal(), TRIANGLE::vertices_[1]->GetNormal(), ave_normal);
	//ARRAY_VECTOR3::add<T>(ave_normal, TRIANGLE::vertices_[2]->GetNormal(), ave_normal);

	//n_[0] = n_[0] + ave_normal[0] / (T)3;
	//n_[1] = n_[1] + ave_normal[1] / (T)3;
	//n_[2] = n_[2] + ave_normal[2] / (T)3;

//	ARRAY_VECTOR3::normalize<T>(n_.values_);
	n_.normalize();
}

void DynamicTriangle::FindIntersection(int d, T *p_input, T *p_intersection, T *u)
{
	TV3 l0 = DynamicTriangle::vertices_[1]->x_ - DynamicTriangle::vertices_[0]->x_.values_;
	TV3 l1 = DynamicTriangle::vertices_[2]->x_ - DynamicTriangle::vertices_[0]->x_.values_;

	T l0x = l0.values_[(d+1)%3];
	T l0y = l0.values_[(d+2)%3];
	T l1x = l1.values_[(d+1)%3];
	T l1y = l1.values_[(d+2)%3];

	T p[2];
	p[0] = p_input[(d+1)%3] - vertices_[0]->x_.values_[(d+1)%3];
	p[1] = p_input[(d+2)%3] - vertices_[0]->x_.values_[(d+2)%3];

	T det = l0x*l1y - l1x*l0y;	// det == 0?? ???쿡?? ?ذ??? ?ȵ?
	det = 1.0f/det;
	u[0] = det*(l1y*p[0]-l1x*p[1]);
	u[1] = det*(-l0y*p[0]+l0x*p[1]);
	
	p_intersection[0] = vertices_[0]->x_.values_[0];
	p_intersection[1] = vertices_[1]->x_.values_[1];
	p_intersection[2] = vertices_[2]->x_.values_[2];
//	ARRAY_VECTOR3::set<T>(p_intersection, vertices_[0]->x_.values_);

	p_intersection[0] += l0.values_[0]*u[0] + l1.values_[0]*u[1];
	p_intersection[1] += l0.values_[1]*u[0] + l1.values_[1]*u[1];
	p_intersection[2] += l0.values_[2]*u[0] + l1.values_[2]*u[1];

	return;
}

void DynamicTriangle::SetEdge(DynamicEdge *e0, DynamicEdge *e1, DynamicEdge *e2)
{
	edges_[0] = e0;
	edges_[1] = e1;
	edges_[2] = e2;
}
