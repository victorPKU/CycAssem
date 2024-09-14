
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846
#define TODEG(A)     ((A)*57.295779513)
#define ROTFLINELEN   100
typedef struct{
	int an;
	float xyz[100][3];
	char name[100];
	float score;
	int N;
	int NC;
	int C;
	int CN;
	int CAN;
	int CAC;
}PROTEIN;



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
void read_structure(char *pn, PROTEIN *p)
{
	FILE *pf;
	char line[LINELEN];
	int n=0;
	int rn=1;
	int l;
	
	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",pn);
		exit(0);
	}
	strcpy(p->name, pn);
	l=strlen(pn);
	p->name[l]='\0';
	p->score=0.0;
	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			getxyz(line,p->xyz[n]);
			if(!strncmp(line+13,"N  ",3)){
				if(rn==1) p->N=n;
				p->CN=n;
			}
			else if(!strncmp(line+13,"CA ",3)){
				if(rn==1) p->CAN=n;
				p->CAC=n;
			}
			else if(!strncmp(line+13,"C  ",3)){
				if(rn==1){
					p->NC=n;
					rn++;
				}
				p->C=n;
			}
			n++;
		}
		else if(!strncmp(line,"REMARK    score",15)){
			p->score=atof(line+15);
		}
	}
	fclose(pf);
	p->an=n;
}
int checkcollision(PROTEIN *p1,PROTEIN *p2)
{
	int ia,ja;
	float dis;

	for(ia=0;ia<p1->an;ia++){
		for(ja=0;ja<p2->an;ja++){
			dis=distance2(p1->xyz[ia], p2->xyz[ja]);
			if(dis<2.5*2.5){
				/*printf("%5d%5d %8.3f\n",ia,ja,sqrt(dis));*/
				if(ia==p1->C&&ja==p2->N) continue;
				if(ja==p2->C&&ia==p1->N) continue;
				if(ia==p1->CAC&&ja==p2->N) continue;
				if(ja==p2->CAC&&ia==p1->N) continue;
				if(ia==p1->C&&ja==p2->CAN) continue;
				if(ja==p2->C&&ia==p1->CAN) continue;
				return 0;
			}
		}
	}
	return 1;
}

#define MAXCADIS (13.5*13.5)
#define MINCADIS (4.7*4.7)
#define AMMINCADIS (3.4*3.4)

#define TERMINALDISCAMAX3   (7.5*7.5)
#define TERMINALDISCAMAX4   (10.5*10.5)

#define MAXFRAGNUM 2000

int main(int argc, char *argv[])
{
	FILE *f1;
	char line0[100];
	char xyzf[100];
	PROTEIN p1[MAXFRAGNUM],p2[MAXFRAGNUM];
	int N1,N2;
	int i;
	int i1,i2;
	int col;
	float dis,w,ang1,ang2;
	float dis2,w2,ang12,ang22;
	float CAdis,NCACA,CACAC,dih;
	float CCACA,CACAN,Dih;
	float CAdis2,NCACA2,CACAC2,dih2;
	float CCACA2,CACAN2,Dih2;
	int isAM1,isAM2;
	int resn;
	float maxtotscore;
	int maxlen;

	if((f1=fopen(argv[1],"r"))==NULL){
		printf("ERROR: Can not open file %s.\n", argv[1]);
		exit(0);
	}
	N1=0;
	while(fgets(line0,100,f1)){
		i=0;while(line0[i]!='.') {xyzf[i]=line0[i]; i++;}
		xyzf[i]='.';xyzf[i+1]='p';xyzf[i+2]='d';xyzf[i+3]='b';xyzf[i+4]='\0';
		read_structure(xyzf, p1+N1);
		N1++;
		if(N1>=MAXFRAGNUM){
			printf("Warning max fragment: %d\n",MAXFRAGNUM);
			break;
		}
	}
	fclose(f1);

	if((f1=fopen(argv[2],"r"))==NULL){
		printf("ERROR: Can not open file %s.\n", argv[2]);
		exit(0);
	}
	N2=0;
	while(fgets(line0,100,f1)){
		i=0;while(line0[i]!='.') {xyzf[i]=line0[i]; i++;}
		xyzf[i]='.';xyzf[i+1]='p';xyzf[i+2]='d';xyzf[i+3]='b';xyzf[i+4]='\0';
		read_structure(xyzf, p2+N2);
		N2++;
		if(N2>=MAXFRAGNUM){
			printf("Warning max fragment: %d\n",MAXFRAGNUM);
			break;
		}
	}
	fclose(f1);

	maxtotscore=atof(argv[3]);
	maxlen=atoi(argv[4]);  /*maximum length of the two linking fragments*/
	for(i1=0;i1<N1;i1++) for(i2=0;i2<N2;i2++){
		if(p1[i1].score+p2[i2].score> maxtotscore) continue;

		col=checkcollision(p1+i1,p2+i2);
		if(col==1){
			CAdis=distance2(p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC]);
			CAdis2=distance2(p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC]);
			if(CAdis>MAXCADIS||CAdis2>MAXCADIS||CAdis<AMMINCADIS||CAdis2<AMMINCADIS){
				continue;
			}
			if(CAdis<MINCADIS){
				isAM1=0;
				dis=distance2(p2[i2].xyz[p2[i2].N], p1[i1].xyz[p1[i1].C]);
				if(dis>1.0*1.0&&dis<1.8*1.8){
					isAM1=1;
					w=cal_dih(p2[i2].xyz[p2[i2].CAN],p2[i2].xyz[p2[i2].N], p1[i1].xyz[p1[i1].C],p1[i1].xyz[p1[i1].CAC]);
					if(w<0) w=-w;
					if(w<135){
						continue;
					}
					ang1=cal_angle(p2[i2].xyz[p2[i2].CAN],p2[i2].xyz[p2[i2].N], p1[i1].xyz[p1[i1].C]);
					ang2=cal_angle(p2[i2].xyz[p2[i2].N], p1[i1].xyz[p1[i1].C],p1[i1].xyz[p1[i1].CAC]);
					if(ang1<90||ang1>150||ang2<90||ang2>150){
						continue;
					}
				}
				else continue;
			}
			if(CAdis2<MINCADIS){
				isAM2=0;
				dis2=distance2(p1[i1].xyz[p1[i1].N], p2[i2].xyz[p2[i2].C]);
				if(dis>1.0*1.0&&dis2<1.8*1.8){
					isAM2=1;
					w2=cal_dih(p1[i1].xyz[p1[i1].CAN],p1[i1].xyz[p1[i1].N], p2[i2].xyz[p2[i2].C],p2[i2].xyz[p2[i2].CAC]);
					if(w2<0) w2=-w2;
					if(w2<135){
						continue;
					}
					ang12=cal_angle(p1[i1].xyz[p1[i1].CAN],p1[i1].xyz[p1[i1].N], p2[i2].xyz[p2[i2].C]);
					ang22=cal_angle(p1[i1].xyz[p1[i1].N], p2[i2].xyz[p2[i2].C],p2[i2].xyz[p2[i2].CAC]);
					if(ang12<90||ang12>150||ang22<90||ang22>150){
						continue;
					}
				}
				else continue;
			}
			resn=0;
			if(CAdis>=MINCADIS){
				if(CAdis<TERMINALDISCAMAX3) resn=1;
				else if(CAdis<TERMINALDISCAMAX4) resn=2;
				else resn=3;
			}
			if(CAdis2>=MINCADIS){
				if(CAdis2<TERMINALDISCAMAX3) resn+=1;
				else if(CAdis2<TERMINALDISCAMAX4) resn+=2;
				else resn+=3;
			}
			if(resn>(maxlen-4)||resn<1) continue;

			printf("%40s %40s %8.3f%8.3f  ",p1[i1].name,p2[i2].name,p1[i1].score,p2[i2].score);
			if(CAdis>=MINCADIS){
				NCACA=cal_angle(p2[i2].xyz[p2[i2].N],p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC]);
				CACAC=cal_angle(p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC],p1[i1].xyz[p1[i1].C]);
				dih=cal_dih(p2[i2].xyz[p2[i2].N],p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC],p1[i1].xyz[p1[i1].C]);

				CCACA=cal_angle(p2[i2].xyz[p2[i2].NC],p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC]);
				CACAN=cal_angle(p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC],p1[i1].xyz[p1[i1].CN]);
				Dih=cal_dih(p2[i2].xyz[p2[i2].NC],p2[i2].xyz[p2[i2].CAN],p1[i1].xyz[p1[i1].CAC],p1[i1].xyz[p1[i1].CN]);
				printf(" 1: %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f ",sqrt(CAdis),NCACA,CACAC,dih,CCACA,CACAN,Dih);
			}
			else{
				if(isAM1==1){			
					printf(" 0: %8.3f%8.3f%8.3f%8.3f%8.3f                 ",sqrt(CAdis),sqrt(dis),w,ang1,ang2);
				}
			}
			if(CAdis2>=MINCADIS){
				NCACA2=cal_angle(p1[i1].xyz[p1[i1].N],p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC]);
				CACAC2=cal_angle(p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC],p2[i2].xyz[p2[i2].C]);
				dih2=cal_dih(p1[i1].xyz[p1[i1].N],p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC],p2[i2].xyz[p2[i2].C]);

				CCACA2=cal_angle(p1[i1].xyz[p1[i1].NC],p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC]);
				CACAN2=cal_angle(p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC],p2[i2].xyz[p2[i2].CN]);
				Dih2=cal_dih(p1[i1].xyz[p1[i1].NC],p1[i1].xyz[p1[i1].CAN],p2[i2].xyz[p2[i2].CAC],p2[i2].xyz[p2[i2].CN]);
				printf("1: %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f \n",sqrt(CAdis2),NCACA2,CACAC2,dih2,CCACA2,CACAN2,Dih2);
			}
			else{
				if(isAM2==1){			
					printf("0: %8.3f%8.3f%8.3f%8.3f%8.3f                 \n",sqrt(CAdis2),sqrt(dis2),w2,ang12,ang22);
				}
			}
		}
	}
}


