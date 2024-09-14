
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TODEG(A)     ((A)*57.295779513)


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
void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+30+8*i);
}

#define LINELEN 81
void readpepterminal(char *pn, float N1[3],float C1[3],float CA1[3],float CA2[3],float CAn[3], float Nn[3], float Cn[3], int len)
{
	FILE *pf;
	char line[LINELEN];
	char res[10]="*********";
	int rn=0;
	
	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",pn);
		exit(0);
	}

	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			if(strncmp(line+17,res,9)){
				rn++;
				strncpy(res,line+17,9);
				res[9]='\0';
			}
			if(!strncmp(line+13,"N  ",3)){
				if(rn==1) getxyz(line,N1);
				if(rn==len) getxyz(line,Nn);
			}
			else if(!strncmp(line+13,"CA ",3)){
				if(rn==1) getxyz(line, CA1);
				else if(rn==2) getxyz(line, CA2);
				else if(rn==len) getxyz(line, CAn);
			}
			else if(!strncmp(line+13,"C  ",3)){
				if(rn==1) getxyz(line,C1);
				if(rn==len) getxyz(line,Cn);
			}
		}
	}
	fclose(pf);
}


int main(int argc, char *argv[])
{
	float N1[3],C1[3], CA1[3],CA2[3], CAn[3], Cn[3],Nn[3];
	float CAdis,N1CACA,CACACn,dih1,NnCACA,CACAC1,dih2;
	int len=atoi(argv[2]);

	readpepterminal(argv[1],N1,C1,CA1,CA2,CAn,Nn,Cn,len);
	CAdis=distance(CA1,CAn);

	N1CACA=cal_angle(N1,CA1,CAn);
	CACACn=cal_angle(CA1,CAn,Cn);
	dih1=cal_dih(N1,CA1,CAn,Cn);

	NnCACA=cal_angle(Nn,CAn,CA1);
	CACAC1=cal_angle(CAn,CA1,C1);
	dih2=cal_dih(Nn,CAn,CA1,C1);

	printf("%s CA-CA %9.3f\n", argv[1],CAdis);
	printf("%s N1-Cn %9.3f%9.3f%9.3f\n", argv[1],N1CACA,CACACn,dih1);
	printf("%s C1-Nn %9.3f%9.3f%9.3f\n", argv[1],NnCACA,CACAC1,dih2);
}


