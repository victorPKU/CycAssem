#include "geometry.h"
#include <math.h>

void getdx(float x1[3],float x2[3],float dx[3])
{
	dx[0]=x1[0]-x2[0];
	dx[1]=x1[1]-x2[1];
	dx[2]=x1[2]-x2[2];
}
float get_len(float x[3])
{
	float r;

	r=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

	return sqrt(r);
}
float distance(float x1[3], float x2[3])
{
	float x12[3];

	getdx(x1,x2,x12);
	return get_len(x12);
}
float distance2(float x1[3], float x2[3])
{
	return (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]);
}
void cpnx(float xi[][3],float xo[][3], int n)
{
	int ia,i;

	for(ia=0;ia<n;ia++){
		for(i=0;i<3;i++)
			xo[ia][i]=xi[ia][i];
	}
}
void cpx(float xi[3],float xo[3])
{
	int i;

	for(i=0;i<3;i++)
		xo[i]=xi[i];
}


float iprod(float a[3],float b[3])
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
float cos_angle(const float a[3],const float b[3])
{
  /* 
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  float   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; 
  
  ip=ipa=ipb=0.0;
  for(m=0; m<3; m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}
float cal_angle(float x1[3], float x2[3], float x3[3])
{
	float x12[3], x32[3];
	float cosa;
	float a;

	getdx(x1,x2,x12);
	getdx(x3,x2,x32);
	cosa=cos_angle(x12,x32);
	a=acos(cosa);

	return TODEG(a);
}
void oprod(const float a[3],const float b[3],float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

float cal_dih(float x1[3], float x2[3], float x3[3], float x4[3])
{
	float x12[3],x32[3],x34[3], m[3], n[3];
	float ipr,phi,cos_phi,sign;
	int i;

	getdx(x1,x2,x12);  
	getdx(x3,x2,x32);	
	getdx(x3,x4,x34);	

	oprod(x12,x32,m); 
	oprod(x32,x34,n);
  	cos_phi=cos_angle(m,n);
  	phi=acos(cos_phi);
	ipr=iprod(x12,n);
	sign=(ipr<0.0)?-1.0:1.0;
	phi=sign*phi; 

  	return TODEG(phi);
}
void get_center(int n, float x[][3], float center[3])
{
	int i,j;
	for(j=0;j<3;j++){
		center[j]=0.0;
	}
	for(i=0;i<n;i++){
		for(j=0;j<3;j++)
			center[j]+=x[i][j];
	}
	for(j=0;j<3;j++){
		center[j]/=n;
	}
}
void direction( float a[3], float b[3], float c[3])
{
	int i;
	float r;
	
	getdx(a,b,c);
	r=get_len(c);
	for(i=0;i<3;i++) 
		c[i]=c[i]/r;
}
void translate(int n, float x0[][3],float x[][3],float v[3])
{
	int i,j;

	for(i=0;i<n;i++){
		for(j=0;j<3;j++){
			x[i][j]=x0[i][j]+v[j];
		}
	}
}
void irtranslate(int n, float x0[][3],float x[][3],float v[3])
{
	int i,j;

	for(i=0;i<n;i++){
		for(j=0;j<3;j++){
			x[i][j]=x0[i][j]-v[j];
		}
	}
}

