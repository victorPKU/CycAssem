/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * watgrid.c: Detect the solvent geometry and chemical environment on the input protein surface by 3D grids.
 ***************************************************************************/

#include "watgrid.h"
#include "geometry.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mem.h"

void pro2grid(float p[3], int g[3])
{
	int i;

	for(i=0;i<3;i++){
		g[i]=(int)(p[i]/GRIDSPACESIZE+GRIDSIZE[i]/2);
	}
}
void grid2xyz(int gx,int gy,int gz, float x[3])
{
	x[0]=(gx-GRIDSIZE[0]/2)*GRIDSPACESIZE;
	x[1]=(gy-GRIDSIZE[1]/2)*GRIDSPACESIZE;
	x[2]=(gz-GRIDSIZE[2]/2)*GRIDSPACESIZE;
}
int getgindex(int gx,int gy,int gz)
{
	return gx*GRIDSIZE[1]*GRIDSIZE[2]+gy*GRIDSIZE[2]+gz;
}

void init_watenegrid(float *grid, int N)
{
	int i;
	for(i=0;i<N;i++)
		grid[i]=WATBULKENE;
}

#define WATSPACE 0
#define PROTEINSPACE 1
void init_mattergrid(int *grid, int N)
{
	int i;
	for(i=0;i<N;i++)
		grid[i]=WATSPACE;
}
void gen_mattergrid(STRUCTURE *pro, int *matbulk)
{
	int ia;
	int i,j,k;
	float rad,rad2;
	int g,gi[3];
	float gxyz[3];
	float dis2;
	int n;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==SURFATOM)
			rad=(pro->a[ia].sol_rad+WATERVDWRAD)*REPSURFVDWR;
		else 
			rad=(pro->a[ia].sol_rad+WATERVDWRAD)*REPCOREVDWR;
		rad2=rad*rad;
		n=(int)(rad/GRIDSPACESIZE);
		pro2grid(pro->a[ia].xyz,gi);
		for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
			for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
				for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
					g=getgindex(i,j,k);
					grid2xyz(i,j,k,gxyz);
					dis2=distance2(pro->a[ia].xyz,gxyz);
					if(dis2<rad2){
						matbulk[g]=PROTEINSPACE;
					}
				}
	}
}

float hbene(float dis, float I, float R)
{
	float r2,r4,r8;

	if(dis<R) return I;
	r2=(R*R)/(dis*dis);
	r4=r2*r2;
	r8=r4*r4;

	return I*r8;
}
void getshift( int i, float shift[3])
{
	shift[0]=(i/9-1)*GRIDSPACESIZE/3.0;
	shift[1]=(i%9/3-1)*GRIDSPACESIZE/3.0;
	shift[2]=(i%3-1)*GRIDSPACESIZE/3.0;
}
#define HBMAXDIS 3.50
#define PROGRID -2
#define BULKGRID -1
#define NOHBBOND 0.0
void init_HBgrid(float *grid, int N)
{
	int i;
	for(i=0;i<N;i++)
		grid[i]=NOHBBOND;
}
float getwatshellthick(float fluc0, float fluc)
{
	float k;

	k=fluc/fluc0;
	if(k<0.6) k=0.6;

	k=exp(1-k);
	
	return WATSHELLTHICK*k;
}
void gen_HBgrid(STRUCTURE *pro, float *HBgrid)
{
	int ia;
	int i,j,k;
	float rad,rad2,surfrad,surfrad2;
	int g,gi[3];
	float gxyz[3];
	float shift[3];
	float dis2,dis;
	int n;
	float *subgrid;
	int *sublab;
	int I,ii;

	snew(subgrid, GRID3SIZE);
	snew(sublab, GRID3SIZE);
	for(i=0;i<GRID3SIZE;i++)
		sublab[i]=0;
	for(I=0;I<27;I++){
		getshift(I,shift);
		for(i=0;i<GRID3SIZE;i++)
			subgrid[i]=BULKGRID;
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->a[ia].surf_core==SURFATOM){
				rad=(pro->a[ia].sol_rad+WATERVDWRAD)*REPSURFVDWR;
				surfrad=(pro->a[ia].sol_rad+getwatshellthick(pro->fluc0,pro->a[ia].fluc));
			}
			else {
				rad=(pro->a[ia].sol_rad+WATERVDWRAD)*REPCOREVDWR;
				surfrad=rad;
			}
			rad2=rad*rad;
			surfrad2=surfrad*surfrad;
			n=(int)(surfrad/GRIDSPACESIZE);
			pro2grid(pro->a[ia].xyz,gi);
			for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
				for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
					for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
						g=getgindex(i,j,k);
						grid2xyz(i,j,k,gxyz);
						for(ii=0;ii<3;ii++) gxyz[ii]+=shift[ii];
						dis2=distance2(pro->a[ia].xyz,gxyz);
						if(dis2<rad2){
							subgrid[g]=PROGRID;
						}
						else if(dis2<surfrad2){
							subgrid[g]=NOHBBOND;
						}
					}
		}
		for(ia=0;ia<pro->atmN;ia++){
			if(pro->a[ia].surf_core==SURFATOM&&pro->a[ia].HBstrength>0.01){
			rad=pro->a[ia].sol_rad+WATERVDWRAD+WATCAGE;
			rad2=rad*rad;
			n=(int)(rad/GRIDSPACESIZE);
			pro2grid(pro->a[ia].xyz,gi);
			for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
				for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
					for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
						g=getgindex(i,j,k);
						if(subgrid[g]<0.0) continue;
						grid2xyz(i,j,k,gxyz);
						for(ii=0;ii<3;ii++) gxyz[ii]+=shift[ii];
						dis2=distance2(pro->a[ia].xyz,gxyz);
						if(dis2<rad2){
							dis=sqrt(dis2);
							subgrid[g]+=hbene(dis,pro->a[ia].HBstrength,HBMAXDIS);
						}
					}
			}
		}
		for(i=0;i<GRID3SIZE;i++){
			if(subgrid[i]>HBgrid[i]){
				HBgrid[i]=subgrid[i];
				sublab[i]=I+1;
			}
		}
	}

	for(i=0;i<GRID3SIZE;i++){
		if(sublab[i]>0){
			HBgrid[i]=sublab[i]*100+HBgrid[i];
		}
	}

	sfree(subgrid);
	sfree(sublab);
}
float proRad(STRUCTURE *p, float span[3])
{
	int ia;
	float R=0.0,d;
	int i;

	for(i=0;i<3;i++) span[i]=0.0;
	for(ia=0;ia<p->atmN;ia++){
		d=(p->a[ia].xyzc[0])*(p->a[ia].xyzc[0])+(p->a[ia].xyzc[1])*(p->a[ia].xyzc[1])+(p->a[ia].xyzc[2])*(p->a[ia].xyzc[2]);
		if(d>R)
			R=d;
		for(i=0;i<3;i++)
			if(fabs(p->a[ia].xyzc[i])>span[i]) span[i]=fabs(p->a[ia].xyzc[i]);
	}
	return sqrt(R);
}

int ps2wGS(float ps)
{
	int n;

	n=(int)((ps+WATERVDWRAD)*2/GRIDSPACESIZE)+1;

	return n;
}

void decide_gridsize(STRUCTURE *p1)
{
	float span1[3];
	int GS1[3];
	int i;

	proRad(p1,span1);
	for(i=0;i<3;i++){
		GS1[i]=ps2wGS(span1[i]+WATERVDWRAD+WATCAGE);
	}
	for(i=0;i<3;i++) GRIDSIZE[i]=GS1[i];
	GRID3SIZE=GS1[0]*GS1[1]*GS1[2];
}
int inpocket(int a,int b,float m)
{
	int ai,aj,ak;
	int bi,bj,bk;
	int x,y,z,d;

	ai=a/(GRIDSIZE[1]*GRIDSIZE[2]);
	bi=b/(GRIDSIZE[1]*GRIDSIZE[2]);
	x=ai-bi;
	if(x<0) x=-x;
	if(x>m) return 0;
	ak=a%(GRIDSIZE[2]);
	bk=b%(GRIDSIZE[2]);
	z=ak-bk;
	if(z<0) z=-z;
	if(z>m) return 0;
	aj=(a-ak)/GRIDSIZE[2]%GRIDSIZE[1];
	bj=(b-bk)/GRIDSIZE[2]%GRIDSIZE[1];
	y=aj-bj;
	if(y<0) y=-y;
	if(y>m) return 0;

	d=x*x+z*z;
	if(d>m*m) return 0;
	d=d+y*y;
	if(d>m*m) return 0;

	return 1;	
}
#define PROPOINT 10
#define OUTPOINT 0
#define POCKETPOINT 3
#define INPOINT   2
#define EDGEPOINT 1
float gdistance2(int i,int j,int k)
{
	return (i*i+j*j+k*k)*GRIDSPACESIZE*GRIDSPACESIZE;
}
#define POCKETDEP  0.7
void getpocket(STRUCTURE *pro, int *mattergrid, int *grid, float RAD)
{
	int ia;
	int i,j,k;
	int ii,jj,kk,gg;
	float rad,rad2;
	int g,gi[3];
	float gxyz[3];
	float dis2;
	int n;
	float gsn;
	int GSN;

	for(i=0;i<GRID3SIZE;i++){
		if(mattergrid[i]==PROTEINSPACE) grid[i]=PROPOINT;
		else grid[i]=OUTPOINT;
	}
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==SURFATOM){
		rad=pro->a[ia].sol_rad+RAD;
		rad2=rad*rad;
		n=(int)(rad/GRIDSPACESIZE);
		pro2grid(pro->a[ia].xyz,gi);
		for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
			for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
				for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
					g=getgindex(i,j,k);
					if(mattergrid[g]==PROTEINSPACE) continue;
					grid2xyz(i,j,k,gxyz);
					dis2=distance2(pro->a[ia].xyz,gxyz);
					if(dis2<=rad2){
						grid[g]=POCKETPOINT;
					}
				}
		}
	}
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==SURFATOM){
		rad=pro->a[ia].sol_rad+RAD+GRIDSPACESIZE*1.75;
		rad2=rad*rad;
		n=(int)(rad/GRIDSPACESIZE);
		pro2grid(pro->a[ia].xyz,gi);
		for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
			for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
				for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
					g=getgindex(i,j,k);
					if(mattergrid[g]==PROTEINSPACE) continue;
					grid2xyz(i,j,k,gxyz);
					dis2=distance2(pro->a[ia].xyz,gxyz);
					if(dis2<=rad2){
						if(grid[g]==OUTPOINT) grid[g]=EDGEPOINT;
					}
				}
		}
	}
	gsn=(RAD)/GRIDSPACESIZE;
	GSN=(int)gsn;
	
	for(i=0;i<GRIDSIZE[0];i++)for(j=0;j<GRIDSIZE[1];j++)for(k=0;k<GRIDSIZE[2];k++){
		g=getgindex(i,j,k);
		if(grid[g]==EDGEPOINT){
			for(ii=BSMAX(0,i-GSN);ii<=BSMIN(GRIDSIZE[0]-1,i+GSN);ii++)
			for(jj=BSMAX(0,j-GSN);jj<=BSMIN(GRIDSIZE[1]-1,j+GSN);jj++)
			for(kk=BSMAX(0,k-GSN);kk<=BSMIN(GRIDSIZE[2]-1,k+GSN);kk++){
				gg=getgindex(ii,jj,kk);
				if(grid[gg]==POCKETPOINT){
					dis2=gdistance2(ii-i,jj-j,kk-k);
					if(dis2<(RAD+POCKETDEP)*(RAD+POCKETDEP)) grid[gg]=INPOINT;
				}
			}
		}
	}
}

void init_pocketgrid(int *mattergrid, int *grid, int N)
{
	int i;
	for(i=0;i<N;i++){
		if(mattergrid[i]==PROTEINSPACE) grid[i]=PROPOINT;
		else grid[i]=WATBULKGEO;
	}
}

void gen_pocketgrid(STRUCTURE *pro,int *pocketgrid, int *mattergrid)
{
	int *grid;
	int i;
	int ia;
	int j,k;
	float rad,rad2;
	int g,gi[3];
	float gxyz[3];
	float dis2;
	int n;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==SURFATOM){
			rad=(pro->a[ia].sol_rad+getwatshellthick(pro->fluc0,pro->a[ia].fluc));
			rad2=rad*rad;
			n=(int)(rad/GRIDSPACESIZE);
			pro2grid(pro->a[ia].xyz,gi);
			for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
				for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
					for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
						g=getgindex(i,j,k);
						if(mattergrid[g]==PROTEINSPACE) continue;
						grid2xyz(i,j,k,gxyz);
						dis2=distance2(pro->a[ia].xyz,gxyz);
						if(dis2<rad2){
							if(pocketgrid[g]==WATBULKGEO) pocketgrid[g]=WATSURFGEO;
						}
					}
		}
	}
	snew(grid, GRID3SIZE);

	getpocket(pro, mattergrid, grid, WATERVDWRAD+WATBULK);	
	for(i=0;i<GRID3SIZE;i++){
		if(grid[i]==POCKETPOINT&&pocketgrid[i]<WATCAGEGEO&&pocketgrid[i]>WATBULKGEO) pocketgrid[i]=WATCAGEGEO;
	}
	getpocket(pro, mattergrid, grid, WATERVDWRAD+WATCAGE);	
	for(i=0;i<GRID3SIZE;i++){
		if(grid[i]==POCKETPOINT&&pocketgrid[i]<WAT10DGEO&&pocketgrid[i]>WATBULKGEO) pocketgrid[i]=WAT10DGEO;
	}
	getpocket(pro, mattergrid, grid, WATERVDWRAD+WAT10D);	
	for(i=0;i<GRID3SIZE;i++){
		if(grid[i]==POCKETPOINT&&pocketgrid[i]<WATPRISMGEO&&pocketgrid[i]>WATBULKGEO) pocketgrid[i]=WATPRISMGEO;
	}
	getpocket(pro, mattergrid, grid, WATERVDWRAD+WATPRISM);	
	for(i=0;i<GRID3SIZE;i++){
		if(grid[i]==POCKETPOINT&&pocketgrid[i]<WATRINGGEO&&pocketgrid[i]>WATBULKGEO) pocketgrid[i]=WATRINGGEO;
	}
	getpocket(pro, mattergrid, grid, WATERVDWRAD+WATRING);	

	sfree(grid);
}
float getwatene(float whb, float hb)
{
	float hbn;
	float t;

	hbn=whb+hb;
	t=-(hbn+1)*(hbn+1)/WATHB_LAM;

	return exp(t)*WATHB_K-WATHB_C;
}
float gethbstrength(float p)
{
	int g;
	
	g=(int)(p/100);
	return p-g*100;
}
int getshiftnum(float p)
{
	int g;

	g=(int)(p/100);

	return g-1;	
}
#define watneighbordis  0.80
int getneighbor(int n[][3])
{
	int i,j,k;
	int MAXN;
	float dis;
	int N;

	MAXN=(int)(watneighbordis/GRIDSPACESIZE);
	N=0;
	for(i=-MAXN;i<MAXN;i++)
	for(j=-MAXN;j<MAXN;j++)
	for(k=-MAXN;k<MAXN;k++){
		dis=(i*i+j*j+k*k)*GRIDSPACESIZE*GRIDSPACESIZE;
		if(dis<=watneighbordis*watneighbordis){
			n[N][0]=i;n[N][1]=j;n[N][2]=k;
			N++;
		}
	}
    return N;
}
float smoothenemap(float *enegrid,int *mattergrid, int *pocketgrid, int i,int j,int k,int neiN,int neighbor[][3])
{
	int g;
	int ni;
	float pmax,nmax;
	
	pmax=0.0; nmax=0.0;

	for(ni=0;ni<neiN;ni++){
		g=getgindex(i+neighbor[ni][0],j+neighbor[ni][1],k+neighbor[ni][2]);
		if(mattergrid[g]==PROTEINSPACE) continue;
		if(pocketgrid[g]==WATBULKGEO) continue;
		if(enegrid[g]>0.0){
			if(enegrid[g]>pmax) pmax=enegrid[g];
		}
		else{
			if((-enegrid[g])>nmax) nmax=-enegrid[g];
		}
	}
	return pmax-nmax;
}
void gen_watenegrid(float *watenegrid, int *mattergrid, int *pocketgrid, float *HBgrid)
{
	int g;
	int i,j,k;
	float *enegrid;
	int neighbor[1000][3];
	int neiN;

	snew(enegrid, GRID3SIZE);
	for(g=0;g<GRID3SIZE;g++){
        enegrid[g]=watenegrid[g];
		if(mattergrid[g]==PROTEINSPACE) continue;
		if(pocketgrid[g]==WATRINGGEO) enegrid[g]=getwatene(WATRINGHBN, gethbstrength(HBgrid[g]));
		else if(pocketgrid[g]==WATPRISMGEO) enegrid[g]=getwatene(WATPRISMHBN, gethbstrength(HBgrid[g]));
		else if(pocketgrid[g]==WAT10DGEO) enegrid[g]=getwatene(WAT10DHBN, gethbstrength(HBgrid[g]));
		else if(pocketgrid[g]==WATCAGEGEO) enegrid[g]=getwatene(WATCAGEHBN, gethbstrength(HBgrid[g]));
		else if(pocketgrid[g]==WATSURFGEO) enegrid[g]=getwatene(WATSURFHBN, gethbstrength(HBgrid[g]));
	}
	neiN=getneighbor(neighbor);
	for(i=0;i<GRIDSIZE[0];i++)
	for(j=0;j<GRIDSIZE[1];j++)
	for(k=0;k<GRIDSIZE[2];k++){
		g=getgindex(i,j,k);
        watenegrid[g]=enegrid[g];
		if(mattergrid[g]==PROTEINSPACE) continue;
		if(pocketgrid[g]==WATBULKGEO) continue;
		watenegrid[g]=smoothenemap(enegrid,mattergrid,pocketgrid,i,j,k,neiN,neighbor);
	}
	sfree(enegrid);
}
void sortenepoint(int N, float *ene,int *index)
{
	int i,j,flag=1;
	int gap=sqrt(N)*2;
	int tmp;

	for(i=0;i<N;i++) index[i]=i;
	while(gap>1){
		gap=gap/2;
		do{
			flag=0;
			for(i=0;i<N-gap;i++){
				j=i+gap;
				if(ene[index[i]]<ene[index[j]]){
					tmp=index[i];
					index[i]=index[j];
					index[j]=tmp;
					flag=1;
				}
			}
		}while(flag!=0);
	}
}
#define SAVEWATENE   0.1

int enegrid2watpositin(int *mattergrid, int *pocketgrid,float *enegrid, float *HBgrid,float enewatxyz[][4], int start, int type)
{
	int pn;
	float (*p)[3];
	float *ene;
	int *index;
	int g;
	int i,j,k;
	int ii;
	int sn;
	float shift[3];
	float gxyz[3];
	int wn;

	pn=0;
	for(g=0;g<GRID3SIZE;g++){
		if(mattergrid[g]!=PROTEINSPACE&&pocketgrid[g]>=WATSURFGEO){
			if(type==1){
				if(enegrid[g]>SAVEWATENE) pn++;
			}
			else if(type==-1){
				if(enegrid[g]<-SAVEWATENE) pn++;
			}
		}
	}
	snew(p,pn);
	snew(ene,pn);
	snew(index,pn);
	pn=0;
	for(i=0;i<GRIDSIZE[0];i++)
	for(j=0;j<GRIDSIZE[1];j++)
	for(k=0;k<GRIDSIZE[2];k++){
		g=getgindex(i,j,k);
		if(mattergrid[g]==PROTEINSPACE) continue;
		if(type==1){
				if(enegrid[g]<SAVEWATENE) continue;
		}
		else if(type==-1){
			if(enegrid[g]>-SAVEWATENE) continue;
		}
		if(pocketgrid[g]>=WATSURFGEO){
			grid2xyz(i,j,k,gxyz);
			sn=getshiftnum(HBgrid[g]);
			getshift(sn,shift);
			for(ii=0;ii<3;ii++) p[pn][ii]=gxyz[ii]+shift[ii];
			ene[pn]=fabs(enegrid[g]);
			pn++;
		}
	}
	sortenepoint(pn,ene,index);

	wn=0;
	for(g=0;g<pn;g++){
		for(k=0;k<wn;k++){
			if(distance2(p[index[g]],enewatxyz[k+start])<WATWATDIS)
				break;
		}
		if(k>=wn){
			for(i=0;i<3;i++) enewatxyz[wn+start][i]=p[index[g]][i];
			enewatxyz[wn+start][3]=ene[index[g]]*type;
			wn++;
		}
	}

	sfree(p);
	sfree(ene);
	sfree(index);

	return wn;
}
void printwatxyz(char *inputnam,char *gfnam,int wn,float watxyz[][4])
{
	FILE *gf;
	int i;

	if((gf=fopen(gfnam,"w"))==NULL){
		printf("ERROR: Can not open water map file %s\n",gfnam);
		exit(0);
	}
	fprintf(gf,"HEADER    water map file for %s\n", inputnam);

	for(i=0;i<wn;i++){
		fprintf(gf,"HETATM%5d  O   HOH A%4d    %8.3f%8.3f%8.3f%6.2f%6.3f\n",
			i+1,i+1,watxyz[i][0],watxyz[i][1],watxyz[i][2],1.0,watxyz[i][3]);
	}
	fclose(gf);
}


