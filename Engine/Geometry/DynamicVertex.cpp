#include "CONVENTIONAL_MACROS.h"

#include "DynamicEdge.h"
#include "DynamicTriangle.h"
#include "DynamicVertex.h"

DynamicVertex::DynamicVertex()
{
	feature_ = false;
	is_mc_vertex_ = false;
	is_boundary_ = false;

	x_ = TV3((T)0, (T)0, (T)0);
	curvature_normal_ = TV3((T)0, (T)0, (T)0);
	deviation_ = TV3((T)0, (T)0, (T)0);
	velocity_ = TV3((T)0, (T)0, (T)0);

	normalizer_ = (T)0;
	normal_deviation_ = (T)0;

	phi_ = 0;
	grain_boundary_number_ = 0;
}

DynamicVertex::DynamicVertex(const T *x)
{
	feature_=false;
	is_mc_vertex_=false;
	is_boundary_=false;
	
	x_.x_ = x[0];
	x_.y_ = x[1];
	x_.z_ = x[2];

	curvature_normal_ = TV3((T)0,(T)0,(T)0);
	deviation_ = TV3((T)0,(T)0,(T)0);
	velocity_ = TV3((T)0,(T)0,(T)0);
	normalizer_=(T)0;
	normal_deviation_=(T)0;

	phi_ = 0;
	grain_boundary_number_ = 0;
}

DynamicVertex::DynamicVertex( TV3* x )
{
	feature_=false;
	is_mc_vertex_=false;
	is_boundary_=false;
	x_ = *x;
	curvature_normal_ = TV3((T)0,(T)0,(T)0);
	deviation_ = TV3((T)0,(T)0,(T)0);
	velocity_ = TV3((T)0,(T)0,(T)0);
	normalizer_=(T)0;
	normal_deviation_=(T)0;

	phi_ = 0;
	grain_boundary_number_ = 0;
}

DynamicVertex::DynamicVertex(T* x,T* n)
{
	feature_=false;
	is_mc_vertex_=false;
	is_boundary_=false;

	x_.x_ = x[0];
	x_.y_ = x[1];
	x_.z_ = x[2];

	n_.x_ = n[0];
	n_.y_ = n[1];
	n_.z_ = n[2];
	
	curvature_normal_ = TV3((T)0,(T)0,(T)0);
	deviation_ = TV3((T)0,(T)0,(T)0);
	velocity_ = TV3((T)0,(T)0,(T)0);
	normalizer_=(T)0;
	normal_deviation_=(T)0;

	phi_ = 0;
	grain_boundary_number_ = 0;
}

DynamicVertex::DynamicVertex( TV3* x, TV3* n )
{
	feature_=false;
	is_mc_vertex_=false;
	is_boundary_=false;

	x_ = *x;
	n_ = *n;

	curvature_normal_ = TV3((T)0,(T)0,(T)0);
	deviation_ = TV3((T)0,(T)0,(T)0);
	velocity_ = TV3((T)0,(T)0,(T)0);
	normalizer_=(T)0;
	normal_deviation_=(T)0;

	phi_ = 0;
	grain_boundary_number_ = 0;
}

DynamicVertex::DynamicVertex(const T x, const T y, const T z, const T nx, const T ny, const T nz)
{
	x_.values_[0] = x;
	x_.values_[1] = y;
	x_.values_[2] = z;

	n_.values_[0] = nx;
	n_.values_[1] = ny;
	n_.values_[2] = nz;
	
	feature_=false;
	is_mc_vertex_=false;
	is_boundary_=false;

	curvature_normal_ = TV3((T)0,(T)0,(T)0);
	deviation_ = TV3((T)0,(T)0,(T)0);
	velocity_ = TV3((T)0,(T)0,(T)0);

	normalizer_=(T)0;
	normal_deviation_=(T)0;

	phi_ = 0;
	grain_boundary_number_ = 0;
}

DynamicVertex::~DynamicVertex()
{
	Reset();
}

void DynamicVertex::Reset()
{
	std::list <DynamicEdge*>::iterator itr_edge;
	for (itr_edge = edges_.begin(); itr_edge != edges_.end(); itr_edge++)
	{
		SAFE_DELETE(*itr_edge);
	}
	edges_.clear();
}

void DynamicVertex::SetPosition(T* x_input)
{
	x_.x_ = x_input[0];
	x_.y_ = x_input[1];
	x_.z_ = x_input[2];
}

void DynamicVertex::SetNormal(T* n_input)
{
	n_.x_ = n_input[0];
	n_.y_ = n_input[1];
	n_.z_ = n_input[2];
}

T* DynamicVertex::GetPosition()
{
	return x_.values_; 
}

T* DynamicVertex::GetNormal()
{
	return n_.values_; 
}

void DynamicVertex::AddEdge(DynamicEdge* edge)
{
	//	triangles.remove(triangle);
	edges_.push_back(edge);
}

void DynamicVertex::DelEdge(DynamicEdge* edge)
{
	edges_.remove(edge);
}

void DynamicVertex::AddTriangle(DynamicTriangle* triangle)
{
	//	triangles.remove(triangle);
	triangles_.push_back(triangle);
}

void DynamicVertex::DelTriangle(DynamicTriangle* triangle)
{
	triangles_.remove(triangle);
}

void DynamicVertex::DetermineNormal()
{
	if((int)triangles_.size() == 0)
	{
		n_.values_[0] = (T)0;
		n_.values_[1] = (T)0;
		n_.values_[2] = (T)1;
		return;
	}

	 n_ = TV3((T)0,(T)0,(T)0);

	TRAVERSE_TRIANGLES
	{
		 n_ = (*itr_triangle)->n_ + n_;
	}

	n_.normalize();

	curvature_normal_ = TV3((T)0,(T)0,(T)0);
	normalizer_=(T)0;
}

void DynamicVertex::DetermineAverageNormal()
{
	if((int)triangles_.size() == 0) return;
/*
	{
		n_.values_[0] = (T)0; 
		n_.values_[1] = (T)0; 
		n_.values_[2] = (T)1; 
		return;
	}
*/
//	T normal_weighted[3] = {(T)0, (T)0, (T)0};
	TV3 normal_weighted = TV3((T)0, (T)0, (T)0);
//	T temp_normal[3] = {(T)0, (T)0, (T)0};

	TRAVERSE_TRIANGLES
	{
//		ARRAY_VECTOR3::add<T>(temp_normal,(*itr_triangle)->GetNormal(),temp_normal);

		TV3 center = (*itr_triangle)->GetTVCenter();

		T dist = (x_-center).getMagnitude();

		TV3 nor = (*itr_triangle)->n_;

		normal_weighted += nor * (T)1/(dist+(T)FLT_EPSILON);
	}

	//	ARRAY_VECTOR3::normalize<T>(temp_normal);

	//n_[0] += (n_[0] + temp_normal[0]) * (T)0.5;
	//n_[1] += (n_[1] + temp_normal[1]) * (T)0.5;
	//n_[2] += (n_[2] + temp_normal[2]) * (T)0.5;

	normal_weighted.normalize();

	n_ = (n_+normal_weighted) * (T)0.5;

	n_.normalize();
}

void DynamicVertex::DetCurvatureNormal()
{
	T A(0);
	TRAVERSE_TRIANGLES
	{
		A+=(*itr_triangle)->area_;
	}
	curvature_normal_ /= (T)4*A;
	// ARRAY_VECTOR3::div<T>(curvature_normal_.values_,(T)4*A); //	TRIANGLE::AddLocalCurvatureNormal() ȣ?? ?Ŀ? ?ҷ??? ?Ѵ?.
	// ARRAY_VECTOR3::div<T>(curvature_normal_,normalizer_);
}

/*
void DynamicVertex::DrawNormal()
{
	glPushMatrix();		
	glTranslatef(x_.values_[0],x_.values_[1],x_.values_[2]); 
	glBegin(GL_LINES);	
#ifdef USE_FLOAT_T
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3fv((float*)GetNormal());
#else
	glVertex3d(0.0,0.0,0.0);
	glVertex3dv(GetNormal());
#endif
	glEnd();
	glPopMatrix();	
}
*/

/*
void DynamicVertex::DrawVelocity()
{
	glPushMatrix();		
	glTranslatef(x_.values_[0],x_.values_[1],x_.values_[2]); 
	glBegin(GL_LINES);
#ifdef USE_FLOAT_T		
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3fv((float*)velocity_.values_); 
#else
	glVertex3d(0.0,0.0,0.0);
	glVertex3dv(velocity_);
#endif

	glEnd();
	glPopMatrix();
}
*/

/*
void DynamicVertex::DrawDeviation()
{
	if((int)triangles_.size()==0)
	{
		return;
	}
	glPushMatrix();
	glTranslatef(x_.values_[0],x_.values_[1],x_.values_[2]); 
	glBegin(GL_LINES);
#ifdef USE_FLOAT_T		
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3fv((float*)deviation_.values_); 
#else
	glVertex3d(0.0,0.0,0.0);
	glVertex3dv(deviation_);	
#endif	
	glEnd();
	glPopMatrix();
}
*/

/*
void DynamicVertex::DrawCurvatureNormal(const T& scalar)
{
	glPushMatrix();		
	glTranslatef(x_.values_[0],x_.values_[1],x_.values_[2]); 
	TV3 n;
	n = curvature_normal_ * scalar;
//	ARRAY_VECTOR3::set<T>(n,curvature_normal_.values_);
//	ARRAY_VECTOR3::mul<T>(n,scalar);
	glBegin(GL_LINES);
#ifdef USE_FLOAT_T		
	glVertex3f(0.0f,0.0f,0.0f);
	glVertex3fv((float*)n.values_);
#else
	glVertex3d(0.0,0.0,0.0);
	glVertex3dv(n);
#endif

	glEnd();
	glPopMatrix();
}
*/

/*
void DynamicVertex::DrawNeighborConnectivity()
{
	TRAVERSE_TRIANGLES
	{
		TV3 v_center;
		v_center = ((*itr_triangle)->vertices_[0]->x_+(*itr_triangle)->vertices_[1]->x_+(*itr_triangle)->vertices_[2]->x_)/(T)3;

		glBegin(GL_LINES);
#ifdef USE_FLOAT_T		
		glVertex3fv((float*)GetPosition());
		glVertex3fv((float*)v_center.values_);
#else
		glVertex3dv(GetPosition());
		glVertex3dv(v_center);
#endif			
		glEnd();
	}
}
*/

void DynamicVertex::Replace(DynamicVertex* v_after)
{
	if(v_after==this)
	{
		std::cout<<"Replace to the same vertex?"<<std::endl;
		return;
	}
	std::cout<<"Vertex replacing";
	TRAVERSE_TRIANGLES
	{		
		DynamicTriangle* triangle=*itr_triangle;
		for(int i=0;i<3;i++)
		{
			if(triangle->vertices_[i]==this)
			{
				triangle->vertices_[i]=v_after;
				v_after->AddTriangle(triangle);
			}
		}
	}
	std::cout<<"End"<<std::endl;
	v_after->triangles_.unique();
}
