#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float distance(float a[3], float b[3])
{
	return (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]);
}

float iprod(float a[3],float b[3])
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

void oprod(const float a[3],const float b[3],float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

void rot_axis(int an, float xi[][3],float xo[][3], float ax[3], float angle)
{
	int ia,i;
	float cosa,sina;
	float dot;
	float cross[3];

	for(ia=0;ia<an;ia++){
		cosa=cos(angle);
		sina=sin(angle);
		dot=iprod(ax,xi[ia]);
		oprod(ax,xi[ia],cross);
		for(i=0;i<3;i++)
			xo[ia][i]=cosa*xi[ia][i]+sina*cross[i]+dot*(1-cosa)*ax[i];
	}
}
float get_len(float x[3])
{
	float r;

	r=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

	return sqrt(r);
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


void getxAxis(float xyz[3][3],float ax[3])
{
	int i;
	float l;


	for(i=0;i<3;i++) { ax[i]=xyz[1][i]-xyz[0][i];}
	l=get_len(ax);
	for(i=0;i<3;i++) 
		 ax[i]=ax[i]/l;
}
void getzAxis(float xyz[3][3],float az[3])
{
	int i;
	float l;
	float ax[3],ay[3];

	for(i=0;i<3;i++) { ax[i]=xyz[1][i];}
	l=get_len(ax);
	for(i=0;i<3;i++) 
		 ax[i]=ax[i]/l;
	for(i=0;i<3;i++) { ay[i]=xyz[2][i];}
	l=get_len(ay);
	for(i=0;i<3;i++) 
		 ay[i]=ay[i]/l;

	oprod(ax,ay,az);
	l=get_len(az);
	for(i=0;i<3;i++) 
		 az[i]=az[i]/l;
}

void rotpro(int an,float xyz[][3],float txyz[3][3])
{
	int i,j;
	float center[3];
	float ax[3],az[3];
	float xaxis[3]={1.0,0.0,0.0};
	float zaxis[3]={0.0,0.0,1.0};
	float rax[3],xang;
	float raz[3],zang;
	float l;

	for(j=0;j<3;j++){
		center[j]=(txyz[0][j]+txyz[1][j])/2.0;
	}
	for(i=0;i<an;i++){
		for(j=0;j<3;j++){
			xyz[i][j]=xyz[i][j]-center[j];
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			txyz[i][j]=txyz[i][j]-center[j];
		}
	}
	getxAxis(txyz,ax);
	oprod(xaxis,ax,rax);
	l=get_len(rax);
	for(i=0;i<3;i++) 
		 rax[i]=rax[i]/l;
	xang=(acos(cos_angle(xaxis,ax)));
	rot_axis(an, xyz,xyz, rax, -xang);	
	rot_axis(3, txyz,txyz, rax, -xang);

	
	getzAxis(txyz,az);
	oprod(zaxis,az,raz);
	l=get_len(raz);
	for(i=0;i<3;i++) 
		 raz[i]=raz[i]/l;
	zang=(acos(cos_angle(zaxis,az)));
	rot_axis(an, xyz,xyz, raz, -zang);	
	for(i=0;i<an;i++){
		for(j=0;j<3;j++){
			xyz[i][j]=xyz[i][j];
		}
	}
	rot_axis(3, txyz,txyz, raz, -zang);

}
void printrotted(char *fn, int an,char nam[][14],float xyz[][3])
{
	FILE *pdbfile;
	int i,n=1;

	pdbfile=fopen(fn,"w");
	for(i=0;i<an;i++){
		fprintf(pdbfile,"ATOM  %5d  %s    %8.3f%8.3f%8.3f\n",n,nam[i],xyz[i][0],xyz[i][1],xyz[i][2]);
		n++;
	}
	fclose(pdbfile);
}
void printboxcenter(float txyz[3][3])
{
	printf("box1       %8.3f%8.3f%8.3f\n",txyz[0][0],txyz[0][1],txyz[0][2]);
	printf("box2       %8.3f%8.3f%8.3f\n",txyz[1][0],txyz[1][1],txyz[1][2]);
	printf("3atm       %8.3f%8.3f%8.3f\n",txyz[2][0],txyz[2][1],txyz[2][2]);

	printf("box2-box1  %8.3f%8.3f%8.3f\n",txyz[1][0]-txyz[0][0],txyz[1][1]-txyz[0][1],txyz[1][2]-txyz[0][2]);
}

#define LINELEN 100
int main(int argc, char *argv[])
{
	FILE *pf;
	char line[LINELEN];
	int an,i;
	float xyz[5000][3];
	char nam[5000][14];
	float txyz[3][3];

	if((pf=fopen(argv[1],"r"))==NULL){
		printf("ERROR: Can not open structure file %s!\n",argv[1]);
		exit(0) ;
	}
	an=0;
	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			strncpy(nam[an], line+13,13);
			nam[an][13]='\0';
			for(i=0;i<3;i++)
				xyz[an][i]=atof(line+30+8*i);
			an++;
		}
	}
	fclose(pf);	
	if((pf=fopen(argv[2],"r"))==NULL){
		printf("ERROR: Can not open structure file %s!\n",argv[2]);
		exit(0) ;
	}
	fgets(line,LINELEN,pf);
	for(i=0;i<3;i++) txyz[0][i]=atof(line+30+8*i);
	fgets(line,LINELEN,pf);
	for(i=0;i<3;i++) txyz[1][i]=atof(line+30+8*i);
	fgets(line,LINELEN,pf);
	for(i=0;i<3;i++) txyz[2][i]=atof(line+30+8*i);

	fclose(pf);

	rotpro(an, xyz, txyz);
	printrotted(argv[3], an, nam,xyz);
	printboxcenter(txyz);
}

