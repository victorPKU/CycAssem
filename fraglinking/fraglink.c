
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mem.h"
#include "geometry.h"
#include "superpose.h"

#define PI 3.14159265358979323846
#define TODEG(A)     ((A)*57.295779513)
typedef struct{
	int start;
	int CA;
	int C;
	int N;
	int end;
	char nam[10];
}RES;
typedef struct{
	int an;
	float xyz[100][3];
	char nam[100][4];
	int resn;
	RES res[5];
	int N;
	int NC;
	int C;
	int CN;
	int CAN;
	int CAC;
}FRAGMENT;


void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+30+8*i);
}
#define LINELEN 81
void readfrag(char *pn, FRAGMENT *p)
{
	FILE *pf;
	char line[LINELEN];
	int n=0;
	int rn=0;
	
	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",pn);
		exit(0);
	}
	
	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			if(rn==0||strncmp(line+17,p->res[rn-1].nam,9)){
				strncpy(p->res[rn].nam,line+17,9);
				p->res[rn].nam[9]='\0';
				p->res[rn].start=n;
				if(rn!=0) p->res[rn-1].end=n-1;
				rn++;
			}
			getxyz(line,p->xyz[n]);
			strncpy(p->nam[n],line+13,3);
			p->nam[n][3]='\0';
			if(!strncmp(line+13,"N  ",3)){
				if(rn==1) p->N=n;
				p->res[rn-1].N=n;
				p->CN=n;
			}
			else if(!strncmp(line+13,"CA ",3)){
				if(rn==1) p->CAN=n;
				p->res[rn-1].CA=n;
				p->CAC=n;
			}
			else if(!strncmp(line+13,"C  ",3)){
				if(rn==1){
					p->NC=n;
				}
				p->res[rn-1].C=n;
				p->C=n;
			}
			n++;
		}
	}
	fclose(pf);
	p->an=n;
	p->resn=rn;
	p->res[rn-1].end=n-1;
}
int readtarget(float target[][3], char *fn)
{
	FILE *pf;
	char line[LINELEN];
	int n=0;
	
	if((pf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",fn);
		exit(0);
	}

	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			getxyz(line,target[n]);
			n++;
		}
	}
	fclose(pf);
	return n;
}
#define TERMINALDISCAMIN   (4.70*4.70)
#define TERMINALDISCAMAX3   (7.5*7.5)
#define TERMINALDISCAMAX4   (10.5*10.5)
#define TERMINALDISCAMAX5   (13.5*13.5)
void getgapgeo(FRAGMENT *p1,FRAGMENT *p2,float *tmdis,float *ncaca,float *cacac,float *dih, float *ncaca2,float *cacac2,float *dih2,float terxyz[2][2][4][3])
{
	int i;

	for(i=0;i<3;i++){
		terxyz[0][0][0][i]=p1->xyz[p1->N][i];
		terxyz[0][0][1][i]=p1->xyz[p1->CAN][i];
		terxyz[0][0][2][i]=p2->xyz[p2->CAC][i];
		terxyz[0][0][3][i]=p2->xyz[p2->C][i];
		terxyz[1][0][0][i]=p2->xyz[p2->N][i];
		terxyz[1][0][1][i]=p2->xyz[p2->CAN][i];
		terxyz[1][0][2][i]=p1->xyz[p1->CAC][i];
		terxyz[1][0][3][i]=p1->xyz[p1->C][i];

		terxyz[0][1][0][i]=p2->xyz[p2->CN][i];
		terxyz[0][1][1][i]=p2->xyz[p2->CAC][i];
		terxyz[0][1][2][i]=p1->xyz[p1->CAN][i];
		terxyz[0][1][3][i]=p1->xyz[p1->NC][i];
		terxyz[1][1][0][i]=p1->xyz[p1->CN][i];
		terxyz[1][1][1][i]=p1->xyz[p1->CAC][i];
		terxyz[1][1][2][i]=p2->xyz[p2->CAN][i];
		terxyz[1][1][3][i]=p2->xyz[p2->NC][i];
	}
	tmdis[0]=distance(p1->xyz[p1->CAN], p2->xyz[p2->CAC]);
	ncaca[0]=cal_angle(p1->xyz[p1->N],p1->xyz[p1->CAN],p2->xyz[p2->CAC]);
	cacac[0]=cal_angle(p1->xyz[p1->CAN],p2->xyz[p2->CAC],p2->xyz[p2->C]);
	dih[0]=cal_dih(p1->xyz[p1->N],p1->xyz[p1->CAN],p2->xyz[p2->CAC],p2->xyz[p2->C] );
	tmdis[1]=distance(p2->xyz[p2->CAN], p1->xyz[p1->CAC]);
	ncaca[1]=cal_angle(p2->xyz[p2->N],p2->xyz[p2->CAN],p1->xyz[p1->CAC]);
	cacac[1]=cal_angle(p2->xyz[p2->CAN],p1->xyz[p1->CAC],p1->xyz[p1->C]);
	dih[1]=cal_dih(p2->xyz[p2->N],p2->xyz[p2->CAN],p1->xyz[p1->CAC],p1->xyz[p1->C] );

	ncaca2[0]=cal_angle(p2->xyz[p2->CN],p2->xyz[p2->CAC],p1->xyz[p1->CAN]);
	cacac2[0]=cal_angle(p2->xyz[p2->CAC],p1->xyz[p1->CAN],p1->xyz[p1->NC]);
	dih2[0]=cal_dih(p2->xyz[p2->CN],p2->xyz[p2->CAC],p1->xyz[p1->CAN],p1->xyz[p1->NC] );
	ncaca2[1]=cal_angle(p1->xyz[p1->CN],p1->xyz[p1->CAC],p2->xyz[p2->CAN]);
	cacac2[1]=cal_angle(p1->xyz[p1->CAC],p2->xyz[p2->CAN],p2->xyz[p2->NC]);
	dih2[1]=cal_dih(p1->xyz[p1->CN],p1->xyz[p1->CAC],p2->xyz[p2->CAN],p2->xyz[p2->NC] );
}

void centerterminal(float terxyz0[2][4][3],float tcenter[2][3],float terxyz[2][4][3])
{
	int i,j,k;

	for(k=0;k<2;k++){
		for(j=0;j<3;j++) tcenter[k][j]=0.0;

		for(i=0;i<4;i++){
			for(j=0;j<3;j++) tcenter[k][j]+=terxyz0[k][i][j];
		}
		for(j=0;j<3;j++) tcenter[k][j]/=4;
		for(i=0;i<4;i++){
			for(j=0;j<3;j++) terxyz[k][i][j]=terxyz0[k][i][j]-tcenter[k][j];
		}
	}
}

#define ROUTELEN 41
void printpepstructure(int id, char *frag1name, char *frag2name,char *fragdir,char *pepname,float rmsd,float steric,float  docking,float tersteric, float score,char pepfname[100], FRAGMENT *p, int patchi, int typ, int ter[2],float fragscore[2])
{
	FILE *pf;
	int atno=1;
	int ir,i;

	if((pf=fopen(pepfname,"a"))==NULL){
		printf("ERROR: Can not create the file %s\n",pepfname);
		exit(0);
	}
	fprintf(pf,"TITLE   %2d %2d %2d\n",patchi+1, id, typ+1);
	fprintf(pf,"REMARK  dock1 %2d %8.3f %s\n",ter[0],fragscore[0], frag1name);
	fprintf(pf,"REMARK  dock2 %2d %8.3f %s\n",ter[1],fragscore[1], frag2name);
	fprintf(pf,"REMARK  link %s/%s\n", fragdir,pepname);
	fprintf(pf,"REMARK  score %8.3f %8.3f%8.3f%8.3f%8.3f\n",score,rmsd,steric,docking,tersteric);

	/*printf(" %s %d\n", pepname, p->resn);*/
	for(ir=0;ir<p->resn;ir++){
		for(i=p->res[ir].start;i<p->res[ir].start+4;i++){
			fprintf(pf,"ATOM %6d  %s %s    %8.3f%8.3f%8.3f\n",atno,p->nam[i],p->res[ir].nam,p->xyz[i][0],p->xyz[i][1],p->xyz[i][2]);
			atno++;
		}
		if(ir!=0&&ir!=p->resn-1){
			for(i=p->res[ir].start+4;i<=p->res[ir].end;i++){		
				fprintf(pf,"ATOM %6d  %s %s    %8.3f%8.3f%8.3f\n",atno,p->nam[i],p->res[ir].nam,p->xyz[i][0],p->xyz[i][1],p->xyz[i][2]);
				atno++;
			}
		}
	}
	fprintf(pf,"END\n");
	fclose(pf);

}

void buildcomplex(float lpcenter[3],float trans[3], float rotR[3][3],FRAGMENT *p)
{
	int ia;
	float transxyz[100][3];
	float fitxyz[100][3];

	irtranslate(p->an,p->xyz ,transxyz,trans);
	do_rot(p->an,transxyz,fitxyz,rotR);
	translate(p->an,fitxyz,p->xyz,lpcenter);
}
void getCAxyz(FRAGMENT *p,float linkerCA[][3])
{
	int ir;
	int i;

	for(ir=0;ir<p->resn;ir++){
		for(i=0;i<3;i++) linkerCA[ir][i]=p->xyz[p->res[ir].CA][i];		
	}
}
void printcoordinate(float xyz1[3], float xyz2[3])
{
	printf("%8.3f%8.3f%8.3f\n",xyz1[0],xyz1[1],xyz1[2]);
	printf("%8.3f%8.3f%8.3f\n",xyz2[0],xyz2[1],xyz2[2]);
}
#define ACCEPTSTERIC 4.0
#define ACCEPTTERSTERIC 1.5
float checksteric(FRAGMENT *p,FRAGMENT *p1,FRAGMENT *p2, int ter[2])
{
	int ia,ja;
	int s,e;
	float d2,d;
	float steric=0.0;
	float tersteric=0.0;

	for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
		for(ja=0;ja<p2->res[p2->resn-1].start;ja++){
			d2=distance2(p->xyz[ia],p2->xyz[ja]);
			if(d2>(3.0*3.0)) steric+=0.0;
			else if(d2<(2.4*2.4)){ steric+=(ACCEPTSTERIC*2);}
			else{
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;/*printcoordinate(p->xyz[ia],p2->xyz[ja]);*/
			}
			if(steric> ACCEPTSTERIC) return steric;
		}
		if(steric> ACCEPTSTERIC) return steric;
	}
	for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
		for(ja=p1->res[1].start;ja<p1->an;ja++){
			d2=distance2(p->xyz[ia],p1->xyz[ja]);
			if(d2>(3.0*3.0)) steric+=0.0;
			else if(d2<(2.4*2.4)){ steric+=(ACCEPTSTERIC*2);}
			else{
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;/*printcoordinate(p->xyz[ia],p2->xyz[ja]);*/
			}
			if(steric> ACCEPTSTERIC) return steric;
		}
		if(steric> ACCEPTSTERIC) return steric;
	}
	tersteric=0.0;ter[0]=1;
	for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
		for(ja=4;ja<p1->res[1].start;ja++){
			d2=distance2(p->xyz[ia],p1->xyz[ja]);
			if(d2>(2.8*2.8)) tersteric+=0.0;
			else if(d2<(2.2*2.2)){ ter[0]=0; break;}
			else{
				d=sqrt(d2);
				tersteric+=1.0-(d-2.2)/0.6;
			}
			if(tersteric> ACCEPTTERSTERIC) { ter[0]=0; break;}
		}
		if(ter[0]==0) break;
	}
	tersteric=0.0;ter[1]=1;
	for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
		for(ja=p2->res[p2->resn-1].start+4;ja<p2->an;ja++){
			d2=distance2(p->xyz[ia],p2->xyz[ja]);
			if(d2>(2.8*2.8)) tersteric+=0.0;
			else if(d2<(2.2*2.2)){ ter[1]=0; break;}
			else{
				d=sqrt(d2);
				tersteric+=1.0-(d-2.2)/0.6;
			}
			if(tersteric> ACCEPTTERSTERIC) { ter[1]=0; break;}
		}
		if(ter[1]==0) break;
	}
	return steric;
}
float checksteric2(FRAGMENT *p,FRAGMENT *p1,FRAGMENT *p2, int ter[2])
{
	int ia,ja;
	int s,e;
	float d2,d;
	float R1,R2;
	float steric=0.0;
	float tersteric=0.0;

	for(ia=p->NC;ia<=p->CN;ia++){
		for(ja=0;ja<p2->res[p2->resn-1].start;ja++){
			if(ja==p2->res[p2->resn-2].C) continue;
			d2=distance2(p->xyz[ia],p2->xyz[ja]);
			if(d2>(3.0*3.0)) steric+=0.0;
			else if(d2<(2.4*2.4)){ steric+=(ACCEPTSTERIC*2);}
			else{
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;/*printcoordinate(p->xyz[ia],p2->xyz[ja]);*/
			}
			if(steric> ACCEPTSTERIC) return steric;
		}
		if(steric> ACCEPTSTERIC) return steric;
	}
	for(ia=p->NC;ia<=p->CN;ia++){
		for(ja=p1->res[1].start+1;ja<p1->an;ja++){
			d2=distance2(p->xyz[ia],p1->xyz[ja]);
			if(d2>(3.0*3.0)) steric+=0.0;
			else if(d2<(2.4*2.4)){ steric+=(ACCEPTSTERIC*2);}
			else{
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;/*printcoordinate(p->xyz[ia],p2->xyz[ja]);*/
			}
			if(steric> ACCEPTSTERIC) return steric;
		}
		if(steric> ACCEPTSTERIC) return steric;
	}
	tersteric=0.0;ter[0]=1;
	d=sqrt(distance2(p->xyz[p->res[p->resn-1].N],p1->xyz[p1->res[0].N]));
	/*printf("%8.3f%8.3f%8.3f %d %8.3f%8.3f%8.3f %d d1=%8.3f \n", p->xyz[p->res[p->resn-1].N][0],p->xyz[p->res[p->resn-1].N][1],p->xyz[p->res[p->resn-1].N][2],p->res[p->resn-1].N, p1->xyz[p1->res[0].N][0],p1->xyz[p1->res[0].N][1],p1->xyz[p1->res[0].N][2], p1->res[0].N, d);*/
	if(d<0.8){
	/*printf("d1=%8.3f \n", d);*/
		for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
			for(ja=4;ja<p1->res[1].start;ja++){
				d2=distance2(p->xyz[ia],p1->xyz[ja]);
				if(d2>(2.8*2.8)) tersteric+=0.0;
				else if(d2<(2.2*2.2)){ ter[0]=0; break;}
				else{
					d=sqrt(d2);
					tersteric+=1.0-(d-2.2)/0.6;
				}
				if(tersteric> ACCEPTTERSTERIC) { ter[0]=0; break;}
			}
			if(ter[0]==0) break;
		}
	}
	else{
		ter[0]=0;
	}
	tersteric=0.0;ter[1]=1;
	d=sqrt(distance2(p->xyz[p->res[0].C],p1->xyz[p1->res[p2->resn-1].C]));
	if(d<0.8){
	/*printf("d2=%8.3f \n", d);*/
		for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
			for(ja=p2->res[p2->resn-1].start+4;ja<p2->an;ja++){
				d2=distance2(p->xyz[ia],p2->xyz[ja]);
				if(d2>(2.8*2.8)) tersteric+=0.0;
				else if(d2<(2.2*2.2)){ ter[1]=0; break;}
				else{
					d=sqrt(d2);
					tersteric+=1.0-(d-2.2)/0.6;
				}
				if(tersteric> ACCEPTTERSTERIC) { ter[1]=0; break;}
			}
			if(ter[1]==0) break;
		}
	}
	else{
		ter[1]=0;
	}
	return steric;
}


#define ACCEPTDOCKSTERIC 5.0

float dock2target(int an,float target[][3],FRAGMENT *p)
{
	int ia,ja;
	float d2,d;
	float steric=0.0;
	float attr=0.0;

	for(ia=p->res[1].start;ia<p->res[p->resn-1].start;ia++){
		for(ja=0;ja<an;ja++){
			d2=distance2(p->xyz[ia],target[ja]);
			/*if(d2>(3.2*3.2)&&d2<4.2*4.2) attr+=1.0;
			else*/ if(d2<3.0*3.0){
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;
			}
			if(steric> ACCEPTDOCKSTERIC) return steric;
		}
	}
	return steric/*-attr*0.1*/;
}

/*#define DISDIFF 0.1
#define ANG1DIFF 10
#define ANG2DIFF 10
#define DIHDIFF  16*/
#define BESTPATCHNUM 500
int add2best2(float score,float *bestscore,int *index, int c)
{
        int i,j;
        int w;

        if(c<BESTPATCHNUM){
                i=c;w=c;
        }
        else{
                if(score>bestscore[index[BESTPATCHNUM-1]]) return -1;
                i=BESTPATCHNUM;w=index[BESTPATCHNUM-1];
        }
        while(i>=1&&score<bestscore[index[i-1]]){
                index[i]=index[i-1];
                i--;
        }
        index[i]=w;
        return w;
}

#define DISDIFF 0.5
#define ANG1DIFF 15
#define ANG2DIFF 15
#define DIHDIFF  35

#define LIBFLINELEN   120

int readfraglib(char *libname, int typ, char patchname[BESTPATCHNUM][56], float patchgeo[BESTPATCHNUM][4],float tmdis,float ncaca,float cacac,float dih)
{
	FILE *gf;
	char disgeofn[100];
	char tn[100];
	char line[LINELEN];
	float d,a1,a2,h;
	
	int dish;
	int N=0;

	float score;
	float bestscore[BESTPATCHNUM];
	int index[BESTPATCHNUM];
	int id,ii;

      for(ii=0;ii<BESTPATCHNUM;ii++){
            index[ii]=ii;
		bestscore[ii]=4;
       }

	for(dish=0;dish<(DISDIFF*2000+1);dish++){
		d=tmdis+dish*0.001-DISDIFF;
		if(typ==0) sprintf(disgeofn,"%s/D%.3f",libname,d);
		else sprintf(disgeofn,"%s/d%.3f",libname,d);
		if((gf=fopen(disgeofn,"r"))==NULL) continue; /*{printf("ERROR: Can not open file %s\n",disgeofn); break;}*/
		while(fgets(line, LIBFLINELEN,gf)){
			a1=atof(line+27);
			if(a1-ncaca<-ANG1DIFF||a1-ncaca>ANG1DIFF) continue;
			a2=atof(line+36);
			if(a2-cacac<-ANG2DIFF||a2-cacac>ANG2DIFF) continue;
			h=atof(line+45);
			if(h-dih<-180) h=h+360;
			else if(h-dih>180) h=h-360;
			if(h-dih<-DIHDIFF||h-dih>DIHDIFF) continue;
			score=((d-tmdis)/DISDIFF)*((d-tmdis)/DISDIFF)+((a1-ncaca)/ANG1DIFF)*((a1-ncaca)/ANG1DIFF)+((a2-cacac)/ANG2DIFF)*((a2-cacac)/ANG2DIFF)+((h-dih)/DIHDIFF)*((h-dih)/DIHDIFF);
			id=add2best2(score,bestscore,index,N);
			if(id!=-1){
				strncpy(patchname[id],line,17);
				patchname[id][17]='\0';
				patchgeo[id][0]=d;patchgeo[id][1]=a1;patchgeo[id][2]=a2;patchgeo[id][3]=h;
			}
			N++;
		}
		fclose(gf);
	}

	if(N>BESTPATCHNUM) return BESTPATCHNUM;
	return N;
}
void readpatchxyz(char *dir,char pnam[56],FRAGMENT *p,float terxyz[4][3],int typ)
{
        FILE *pf;
        char line[LINELEN];
        char res[10]="*********";
	int n=0;
	int rn=0;
	int i;
	char nam[100];
	
	sprintf(nam, "%s/%s",dir,pnam);

	if((pf=fopen(nam,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",nam);
		exit(0);
	}

	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){			
			if(!strncmp(line+13,"N  ",3)||!strncmp(line+13,"CA ",3)||!strncmp(line+13,"C  ",3)||!strncmp(line+13,"O  ",3)){
			if(rn==0||strncmp(line+17,p->res[rn-1].nam,9)){
				strncpy(p->res[rn].nam,line+17,9);
				p->res[rn].nam[9]='\0';
				p->res[rn].start=n;
				if(rn!=0) p->res[rn-1].end=n-1;
				rn++;
			}
			getxyz(line,p->xyz[n]);
			strncpy(p->nam[n],line+13,3);
			p->nam[n][3]='\0';
			if(!strncmp(line+13,"N  ",3)){
				if(rn==1) p->N=n;
				p->res[rn-1].N=n;
				p->CN=n;
			}
			else if(!strncmp(line+13,"CA ",3)){
				if(rn==1) p->CAN=n;
				p->res[rn-1].CA=n;
				p->CAC=n;
			}
			else if(!strncmp(line+13,"C  ",3)){
				if(rn==1){
					p->NC=n;
				}
				p->res[rn-1].C=n;
				p->C=n;
			}
			n++;
			}
		}
	}
	fclose(pf);
	p->an=n;
	p->resn=rn;
	p->res[rn-1].end=n-1;

	if(typ==0){
		for(i=0;i<3;i++){
			terxyz[0][i]=p->xyz[p->CN][i];
			terxyz[1][i]=p->xyz[p->CAC][i];
			terxyz[2][i]=p->xyz[p->CAN][i];
			terxyz[3][i]=p->xyz[p->NC][i];
		}
	}
	else{
		for(i=0;i<3;i++){
			terxyz[0][i]=p->xyz[p->N][i];
			terxyz[1][i]=p->xyz[p->CAN][i];
			terxyz[2][i]=p->xyz[p->CAC][i];
			terxyz[3][i]=p->xyz[p->C][i];
		}
	}
}
#define BESTNUM 30
#define BESTCLUSTERNUM 3
#define CULSTERRMSD 0.6
int CAcompare(int resn,float linkerCA[5][3], int N, int bestCAn[BESTCLUSTERNUM], float bestCA[BESTCLUSTERNUM][5][3])
{
	int clusterid;

	int ic,ir,i;
	float rmsd;

	for(ic=0;ic<N;ic++){
		rmsd=0.0;
		if(resn!=bestCAn[ic]) continue;
		for(ir=0;ir<resn;ir++){
			rmsd+=(bestCA[ic][ir][0]-linkerCA[ir][0])*(bestCA[ic][ir][0]-linkerCA[ir][0])+(bestCA[ic][ir][1]-linkerCA[ir][1])*(bestCA[ic][ir][1]-linkerCA[ir][1])+(bestCA[ic][ir][2]-linkerCA[ir][2])*(bestCA[ic][ir][2]-linkerCA[ir][2]);
		}
		rmsd=sqrt(rmsd/resn);
		if(rmsd<CULSTERRMSD) return ic+1;
	}
	return ic+1;
}
void cp2clusterCA(int resn,float linkerCA[5][3], float bestCA[5][3])
{
	int ir,i;

	for(ir=0;ir<resn;ir++) for(i=0;i<3;i++) bestCA[ir][i]=linkerCA[ir][i];
}


int add2best(float score,float *bestscore,int *index, int c)
{
        int i,j;
        int w;

        if(c<BESTNUM){
                i=c;w=c;
        }
        else{
                if(score>bestscore[index[BESTNUM-1]]) return -1;
                i=BESTNUM;w=index[BESTNUM-1];
        }
        while(i>=1&&score<bestscore[index[i-1]]){
                index[i]=index[i-1];
                i--;
        }
        index[i]=w;
        return w;
}
void pepterdock(int an,float target[][3],FRAGMENT *p1,FRAGMENT *p2, float score[4])
{
	int ia,ja;
	float d2,d;
	float steric=0.0;
	float attr=0.0;

	for(ia=4;ia<p1->res[1].start;ia++){
		for(ja=0;ja<an;ja++){
			d2=distance2(p1->xyz[ia],target[ja]);
			if(d2>(3.2*3.2)&&d2<4.2*4.2) attr+=1.0;
			else if(d2<3.0*3.0){
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;
			}
		}
	}
	score[0]=steric*0.4-attr*0.1;
	if(score[0]>0.0) score[0]=0.0;

	steric=0.0;
	attr=0.0;
	for(ia=p2->res[p2->resn-1].start+4;ia<p2->an;ia++){
		for(ja=0;ja<an;ja++){
			d2=distance2(p2->xyz[ia],target[ja]);
			if(d2>(3.2*3.2)&&d2<4.2*4.2) attr+=1.0;
			else if(d2<3.0*3.0){
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;
			}
		}
	}
	score[3]=steric*0.4-attr*0.1;
	if(score[3]>0.0) score[3]=0.0;

	steric=0.0;
	attr=0.0;
	for(ia=4;ia<p2->res[1].start;ia++){
		for(ja=0;ja<an;ja++){
			d2=distance2(p2->xyz[ia],target[ja]);
			if(d2>(3.2*3.2)&&d2<4.2*4.2) attr+=1.0;
			else if(d2<3.0*3.0){
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;
			}
		}
	}
	score[2]=steric*0.4-attr*0.1;
	if(score[2]>0.0) score[2]=0.0;

	steric=0.0;
	attr=0.0;
	for(ia=p1->res[p1->resn-1].start+4;ia<p1->an;ia++){
		for(ja=0;ja<an;ja++){
			d2=distance2(p1->xyz[ia],target[ja]);
			if(d2>(3.2*3.2)&&d2<4.2*4.2) attr+=1.0;
			else if(d2<3.0*3.0){
				d=sqrt(d2);
				steric+=1.0-(d-2.4)/0.6;
			}
		}
	}
	score[1]=steric*0.4-attr*0.1;
	if(score[1]>0.0) score[1]=0.0;
}

#define ACCEPTRMSD 0.8 
#define ACCEPTDOCKING 2.0

int main(int argc, char *argv[])
{
	FRAGMENT p1,p2;
	int an;
	float target[2000][3];
	float tmdis[2];
	float ncaca[2];
	float cacac[2];
	float dih[2];
	float ncaca2[2];
	float cacac2[2];
	float dih2[2];

	char patchname[2][2][3][BESTPATCHNUM][56];
	float patchgeo[2][2][3][BESTPATCHNUM][4];
	int matchN[2][2][3];
	int im,il;

	FRAGMENT patch;

	float terxyz0[2][2][4][3];
	float terxyz[2][2][4][3];
	float tcenter[2][2][3];
	float pterxyz[4][3];

	float pepterdockscore[4];
	float rmsd, steric,docking;
	int tersteric[2];
	float terminalscore;
	float score;

	float trans[3];
	float rotR[3][3];

	int ii,jj,kk,N[2],pt[2][BESTNUM][3];
	float best[2][BESTNUM];
	float linkerCA[5][3];
	int bestCAn[2][BESTNUM];
	float bestCA[2][BESTNUM][5][3];
	int clusterid,clusterN;
	int index[2][BESTNUM];
	int id;

	FILE *f1;
	char line0[300];
	float fragscore[2];
	char frg1nm[100];
	char frg2nm[100];
	char outname[100];
	int i,j;

	int NN=0;

	an=readtarget(target, argv[2]);
	if((f1=fopen(argv[1],"r"))==NULL){
		printf("ERROR: Can not open file %s.\n", argv[1]);
		exit(0);
	}
	while(fgets(line0,300,f1)){
		NN++;
		sprintf(outname, "%s_%d",argv[3],NN);
		i=0;while(line0[i]==' ') i++;
		j=0;while(line0[i]!=' ') {frg1nm[j]=line0[i];i++;j++;}
		frg1nm[j]='\0';
		
		while(line0[i]==' ') i++;
		j=0;while(line0[i]!=' ') {frg2nm[j]=line0[i];i++;j++;}
		frg2nm[j]='\0';

		fragscore[0]=atof(line0+i+1);
		fragscore[1]=atof(line0+i+9);

	readfrag(frg1nm, &p1);
	readfrag(frg2nm, &p2);

	getgapgeo(&p1,&p2,tmdis,ncaca,cacac,dih,ncaca2,cacac2,dih2,terxyz0);
	centerterminal(terxyz0[0],tcenter[0],terxyz[0]);
	centerterminal(terxyz0[1],tcenter[1],terxyz[1]);
	pepterdock(an,target,&p1,&p2, pepterdockscore);

	for(ii=0;ii<2;ii++){
		if(tmdis[ii]<7.5&&tmdis[ii]>4.5){
			matchN[ii][0][0]=readfraglib(argv[4], 0,patchname[ii][0][0], patchgeo[ii][0][0],tmdis[ii],ncaca[ii],cacac[ii],dih[ii]);
		}
		else
			matchN[ii][0][0]=0;
		if(tmdis[ii]<10.5&&tmdis[ii]>5.5){
			matchN[ii][0][1]=readfraglib(argv[5], 0,patchname[ii][0][1], patchgeo[ii][0][1],tmdis[ii],ncaca[ii],cacac[ii],dih[ii]);
		}
		else
			matchN[ii][0][1]=0;
		if(tmdis[ii]<13.5&&tmdis[ii]>6.5){
			matchN[ii][0][2]=readfraglib(argv[6], 0,patchname[ii][0][2], patchgeo[ii][0][2],tmdis[ii],ncaca[ii],cacac[ii],dih[ii]);
		}
		else
			matchN[ii][0][2]=0;
		if(tmdis[ii]<7.5&&tmdis[ii]>4.5){
			matchN[ii][1][0]=readfraglib(argv[4], 1,patchname[ii][1][0], patchgeo[ii][1][0],tmdis[ii],ncaca2[ii],cacac2[ii],dih2[ii]);
		}
		else
			matchN[ii][1][0]=0;
		if(tmdis[ii]<10.5&&tmdis[ii]>5.5){
			matchN[ii][1][1]=readfraglib(argv[5], 1,patchname[ii][1][1], patchgeo[ii][1][1],tmdis[ii],ncaca2[ii],cacac2[ii],dih2[ii]);
		}
		else
			matchN[ii][1][1]=0;
		if(tmdis[ii]<13.5&&tmdis[ii]>6.5){
			matchN[ii][1][2]=readfraglib(argv[6], 1,patchname[ii][1][2], patchgeo[ii][1][2],tmdis[ii],ncaca2[ii],cacac2[ii],dih2[ii]);
		}
		else
			matchN[ii][1][2]=0;
		if((matchN[ii][0][0]+matchN[ii][0][1]+matchN[ii][0][2]+matchN[ii][1][0]+matchN[ii][1][1]+matchN[ii][1][2])==0) break;
	}
	if(ii<2) continue;
	
      for(ii=0;ii<BESTNUM;ii++){
            index[0][ii]=ii;
            index[1][ii]=ii;
		best[0][ii]=ACCEPTDOCKING+ACCEPTSTERIC+ACCEPTRMSD*4+5;
 		best[1][ii]=ACCEPTDOCKING+ACCEPTSTERIC+ACCEPTRMSD*4+5;
       }
	for(ii=0;ii<2;ii++){
		N[ii]=0;
		for(jj=0;jj<2;jj++) for(il=0;il<3;il++){
			for(im=0;im<matchN[ii][jj][il];im++){
				readpatchxyz(argv[7],patchname[ii][jj][il][im],&patch,pterxyz,jj);
				rmsd=fit_terminal(terxyz[ii][jj],pterxyz,trans, rotR);
				/*printf("%s %8.3f\n", patchname[il][im],rmsd);*/
				if(rmsd<ACCEPTRMSD){
					buildcomplex(tcenter[ii][jj],trans, rotR,&patch);
					if(ii==0&&jj==0) steric=checksteric( &patch,&p1,&p2,tersteric);
					else if(ii==1&&jj==0)  steric=checksteric( &patch,&p2,&p1,tersteric);
					else if(ii==0&&jj==1)  steric=checksteric2( &patch,&p1,&p2,tersteric);
					else if(ii==1&&jj==1)  steric=checksteric2( &patch,&p2,&p1,tersteric);
					/*printf("%s %8.3f\n", patchname[il][im],steric);*/
					if(steric<ACCEPTSTERIC) {
						docking=dock2target(an,target,&patch);
						if(docking < ACCEPTDOCKING){
							terminalscore=0.0;
							if(ii==0){
								if(tersteric[0]==1) terminalscore+=pepterdockscore[0];
								if(tersteric[1]==1) terminalscore+=pepterdockscore[3];
							}
							else{
								if(tersteric[0]==1) terminalscore+=pepterdockscore[2];
								if(tersteric[1]==1) terminalscore+=pepterdockscore[1];
							}
							score=(rmsd)*4+(steric)+(docking)+il*2+terminalscore;
							id=add2best(score,best[ii],index[ii],N[ii]);
							if(id!=-1){ pt[ii][id][0]=il;pt[ii][id][1]=im;pt[ii][id][2]=jj;best[ii][id]=score;}
							N[ii]++;
						}
					}
				}
			}
		}
		if(N[ii]==0) break;
		else{
			if(N[ii]>BESTNUM) N[ii]=BESTNUM;
		}
	}
	if(ii<2) continue;
	
	printf("%5d %s %s\n",NN,frg1nm,frg2nm);		
	for(ii=0;ii<2;ii++){
		clusterN=0;
		for(kk=0;kk<N[ii];kk++){
			id=index[ii][kk];
			il=pt[ii][id][0];
			im=pt[ii][id][1];
			jj=pt[ii][id][2];
			if(best[ii][id]-best[ii][index[ii][0]]>2) break;
			readpatchxyz(argv[7],patchname[ii][jj][il][im],&patch,pterxyz,jj);
			rmsd=fit_terminal(terxyz[ii][jj],pterxyz,trans, rotR);
			buildcomplex(tcenter[ii][jj],trans, rotR,&patch);
			if(ii==0&&jj==0) steric=checksteric( &patch,&p1,&p2,tersteric);
			else if(ii==1&&jj==0)  steric=checksteric( &patch,&p2,&p1,tersteric);
			else if(ii==0&&jj==1)  steric=checksteric2( &patch,&p1,&p2,tersteric);
			else if(ii==1&&jj==1)  steric=checksteric2( &patch,&p2,&p1,tersteric);
			docking=dock2target(an,target,&patch);
			terminalscore=0.0;
			if(ii==0){
				if(tersteric[0]==1) terminalscore+=pepterdockscore[0];
				if(tersteric[1]==1) terminalscore+=pepterdockscore[3];
			}
			else{
				if(tersteric[0]==1) terminalscore+=pepterdockscore[2];
				if(tersteric[1]==1) terminalscore+=pepterdockscore[1];
			}
			score=(rmsd)*4+(steric)+(docking)+il*2+terminalscore;
			getCAxyz(&patch,linkerCA);
			clusterid=CAcompare(patch.resn,linkerCA,clusterN,bestCAn[ii], bestCA[ii]);
			/*printf("%5d%5d %8.3f%8.3f %d %d\n",ii,kk, score, best[ii][id], clusterid, clusterN);*/
			if(clusterid>clusterN){
				printpepstructure(clusterid, frg1nm, frg2nm,argv[7],patchname[ii][jj][il][im],rmsd,steric,docking, terminalscore, score, outname, &patch,ii,jj,tersteric,fragscore);
				bestCAn[ii][clusterN]=patch.resn;
				cp2clusterCA(patch.resn,linkerCA, bestCA[ii][clusterN]);
				clusterN++;	
				if(clusterN>=BESTCLUSTERNUM) break;
			}
		}
	}
	}
	fclose(f1);

	return 1;
}
