/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * grid.c: Map the proteins to interaction grids.
 ***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "grid.h"
#include "geometry.h"

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

/*#define MAXNEIGHB    12
#define NEIGHBGRID   1.548
#define EXPCOEF      1.925
#define MAXNEIGHB    10
#define NEIGHBGRID   1.469
#define EXPCOEF      2.135*/

/*#define MAXNEIGHB    10
#define NEIGHBGRID   1.40
#define EXPCOEF      2.35*/

#define MAXNEIGHB   8
#define NEIGHBGRID  1.19
#define EXPCOEF     3.25

int getnearest(float x[3],int gn[MAXNEIGHB],float coe[MAXNEIGHB])
{
	float gxyz[3];
	float dis2;
	int gi[3];
	int i,j,k,n;
	int N=0;
	
	n=(int)(NEIGHBGRID/GRIDSPACESIZE);
	pro2grid(x,gi);
	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				grid2xyz(i,j,k,gxyz);
				dis2=distance2(x,gxyz);
				if(dis2<NEIGHBGRID*NEIGHBGRID){
					if(N<MAXNEIGHB){
						gn[N]=getgindex(i,j,k);
						coe[N]=exp(-EXPCOEF*dis2);
						N++;
					}
				}
			}

	return N;
}
void init_rgrid(float *grid, int N)
{
	int i;
	for(i=0;i<N;i++)
		grid[i]=0.0;
}
void init_grid(fftwf_complex *grid, int N)
{
	int i;
	for(i=0;i<N;i++){
		grid[i][0]=0.0;
		grid[i][1]=0.0;
	}
}
void init_matter(int *grid, int N)
{
	int i;
	for(i=0;i<N;i++)
		grid[i]=0;
}
void gen_matter_region(STRUCTURE *pro, int *matter)
{
	int ia;
	int i,j,k;
	float rad,rad2;
	int g,gi[3];
	float gxyz[3];
	float dis2;
	int n;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==COREATOM){
			rad=(pro->a[ia].vdw_rad*REPCOREVDWR+AVERAGEVDWRAD*REPSURFVDWR);
		}
		else{
			rad=(pro->a[ia].vdw_rad+AVERAGEVDWRAD)*REPSURFVDWR;
		}
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
						if(matter[g]==0)
							matter[g]=1;
					}
				}
	}
}
void atom_rep_grid(float xyz[3],float rad,float score,float *rep)
{
	int i,j,k;
	float rad2;
	int g,gi[3];
	float gxyz[3];
	int n;
	
	rad2=rad*rad;
	n=(int)(rad/GRIDSPACESIZE);
	pro2grid(xyz,gi);
	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				g=getgindex(i,j,k);
				if(rep[g]>score)
					continue;
				grid2xyz(i,j,k,gxyz);
				if(distance2(xyz,gxyz)<rad2)
					rep[g]=score;
			}
}

void gen_rep_grid(STRUCTURE *pro, float *rep)
{
	int ia;
	float score;
	float rad;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==COREATOM){
			rad=pro->a[ia].vdw_rad*REPCOREVDWR;
			score=3.0;
		}
		else{
			rad=pro->a[ia].vdw_rad*REPSURFVDWR;
			score=1.0;
		}
		atom_rep_grid(pro->a[ia].xyz, rad, score,rep);
	}
}

void gen_matter_grid(STRUCTURE *pro, fftwf_complex *vdwgrid, fftwf_complex *elegrid, fftwf_complex *solgrid1, fftwf_complex *solgrid2, fftwf_complex *solgrid3)
{
	int ia;
	int g[MAXNEIGHB];
	int ii;
	float coe[MAXNEIGHB];
	int nn;
	float coet;

	for(ia=0;ia<pro->atmN;ia++){
		nn=getnearest(pro->a[ia].xyz,g,coe);
		coet=0.0;
		for(ii=0;ii<nn;ii++) coet+=coe[ii];
		for(ii=0;ii<nn;ii++){
			vdwgrid[g[ii]][0]+=pro->a[ia].sq_vdw_eps*coe[ii]/coet;
			if(pro->a[ia].charge>0.01||pro->a[ia].charge<-0.01)
				elegrid[g[ii]][0]+=pro->a[ia].charge*coe[ii]/coet;
			solgrid1[g[ii]][0]+=pro->a[ia].volume*coe[ii]/coet;
			if(pro->a[ia].soltype==1||pro->a[ia].soltype==2){
				solgrid2[g[ii]][0]+=pro->a[ia].gfree/3.5*coe[ii]/coet;
			}
			else if(pro->a[ia].soltype==3)
				solgrid3[g[ii]][0]+=pro->a[ia].gfree/6.0*coe[ii]/coet;
		}
	}
}
#define MAXDIS 7.2
void gen_atom_vdw_field(float xyz[3],float rad, float sq_vdw_eps, int *matter,fftwf_complex *vdwgrid)
{
	int i,j,k;
	int n;
	float Fdis2,Sdis2,dis2;
	int g,gi[3];
	float gxyz[3];

	Fdis2=VDWFIRSTSHELL*VDWFIRSTSHELL*(rad+AVERAGEVDWRAD)*(rad+AVERAGEVDWRAD);
	Sdis2=VDWSECONDSHELL*VDWSECONDSHELL*(rad+AVERAGEVDWRAD)*(rad+AVERAGEVDWRAD);
	n=(int)(VDWSECONDSHELL*(rad+AVERAGEVDWRAD)/GRIDSPACESIZE);

	pro2grid(xyz, gi);

	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				g=getgindex(i,j,k);
				if(matter[g]!=0)
					continue;		
				grid2xyz(i,j,k,gxyz);
				dis2=distance2(xyz,gxyz);
				if(dis2<Fdis2){
					vdwgrid[g][1]+=(-ATTRSHORT)*sq_vdw_eps;
				}
				else if(dis2<Sdis2){
					vdwgrid[g][1]+=(-ATTRLONG)*sq_vdw_eps;
				}
			}
}
void gen_atom_ele_field(float xyz[3],float occupy,float charge,int *matter,fftwf_complex *elegrid)
{
	int i,j,k;
	int n;
	float dis2;
	float s1,s2,s3;
	int g,gi[3];
	float gxyz[3];
	
	s1=ELE1SHELL*ELE1SHELL;
	s2=ELE2SHELL*ELE2SHELL;
	s3=ELE3SHELL*ELE3SHELL;
	n=(int)(ELE3SHELL/GRIDSPACESIZE);
	pro2grid(xyz, gi);
	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				g=getgindex(i,j,k);
				if(matter[g]!=0)
					continue;		
				grid2xyz(i,j,k,gxyz);
				dis2=distance2(xyz,gxyz);
				if(dis2<s1){
					elegrid[g][1]+=ELE_1*occupy*charge;
				}
				else if(dis2<s2){
					elegrid[g][1]+=ELE_2*occupy*charge;
				}
				else if(dis2<s3){
					elegrid[g][1]+=ELE_3*occupy*charge;
				}
			}
}
void gen_atom_gfree_field(float xyz[3],int soltype, float gfree, int *matter,fftwf_complex *solgrid1)
{
	int i,j,k;
	int n;
	float Fdis2,Sdis2,dis2;
	int g,gi[3];
	float gxyz[3];

	if(soltype==1||soltype==3){
		Fdis2=SOL1FIRSTSHELL*SOL1FIRSTSHELL;
		Sdis2=SOL1SECONDSHELL*SOL1SECONDSHELL;
		n=(int)(SOL1SECONDSHELL/GRIDSPACESIZE);
	}
	else {
		Fdis2=SOL2FIRSTSHELL*SOL2FIRSTSHELL;
		Sdis2=SOL2SECONDSHELL*SOL2SECONDSHELL;
		n=(int)(SOL2SECONDSHELL/GRIDSPACESIZE);
	}

	pro2grid(xyz, gi);

	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				g=getgindex(i,j,k);
				if(matter[g]!=0)
					continue;		
				grid2xyz(i,j,k,gxyz);
				dis2=distance2(xyz,gxyz);
				if(dis2<Fdis2){
					if(soltype==1){
						solgrid1[g][1]+=(-SOL1SHORT_3p5/3.5)*gfree;
					}
					else if(soltype==2){
						solgrid1[g][1]+=(-SOL2SHORT_3p5/3.5)*gfree;
					}
					else if(soltype==3){
						solgrid1[g][1]+=(-SOL1SHORT_6p0/6.0)*gfree;
					}
				}
				else if(dis2<Sdis2){
					if(soltype==1){
						solgrid1[g][1]+=(-SOL1LONG_3p5/3.5)*gfree;
					}
					else if(soltype==2){
						solgrid1[g][1]+=(-SOL2LONG_3p5/3.5)*gfree;
					}
					else if(soltype==3){
						solgrid1[g][1]+=(-SOL1LONG_6p0/6.0)*gfree;
					}
				}
			}
}
void gen_atom_volum_field(float xyz[3],int soltype, float volum, int *matter,fftwf_complex *solgrid2,fftwf_complex *solgrid3)
{
	int i,j,k;
	int n;
	float Fdis2,Sdis2,dis2;
	int g,gi[3];
	float gxyz[3];

	if(soltype==1||soltype==3){
		Fdis2=SOL1FIRSTSHELL*SOL1FIRSTSHELL;
		Sdis2=SOL1SECONDSHELL*SOL1SECONDSHELL;
		n=(int)(SOL1SECONDSHELL/GRIDSPACESIZE);
	}
	else{
		Fdis2=SOL2FIRSTSHELL*SOL2FIRSTSHELL;
		Sdis2=SOL2SECONDSHELL*SOL2SECONDSHELL;
		n=(int)(SOL2SECONDSHELL/GRIDSPACESIZE);
	}

	pro2grid(xyz, gi);

	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				g=getgindex(i,j,k);
				if(matter[g]!=0)
					continue;		
				grid2xyz(i,j,k,gxyz);
				dis2=distance2(xyz,gxyz);
				if(dis2<Fdis2){
					if(soltype==1||soltype==3){
						solgrid2[g][1]+=(-SOL1SHORT_3p5)*volum;
						solgrid3[g][1]+=(-SOL1SHORT_6p0)*volum;
					}
					else if(soltype==2){
						solgrid2[g][1]+=(-SOL2SHORT_3p5)*volum;
						solgrid3[g][1]+=(-SOL2SHORT_6p0)*volum;
					}
				}
				else if(dis2<Sdis2){
					if(soltype==1||soltype==3){
						solgrid2[g][1]+=(-SOL1LONG_3p5)*volum;
						solgrid3[g][1]+=(-SOL1LONG_6p0)*volum;
					}
					else if(soltype==2){
						solgrid2[g][1]+=(-SOL2LONG_3p5)*volum;
						solgrid3[g][1]+=(-SOL2LONG_6p0)*volum;
					}
				}
			}
}
void gen_field_grid(STRUCTURE *pro, int *matter,fftwf_complex *vdwgrid, fftwf_complex *elegrid, fftwf_complex *solgrid1, fftwf_complex *solgrid2,fftwf_complex *solgrid3)
{
	int ia;

	for(ia=0;ia<pro->atmN;ia++){
		gen_atom_vdw_field(pro->a[ia].xyz,pro->a[ia].vdw_rad,pro->a[ia].sq_vdw_eps, matter,vdwgrid);
		if(pro->a[ia].charge>0.01||pro->a[ia].charge<-0.01)
			gen_atom_ele_field(pro->a[ia].xyz,1.0,pro->a[ia].charge,matter,elegrid);
		gen_atom_volum_field(pro->a[ia].xyz,pro->a[ia].soltype, pro->a[ia].volume, matter,solgrid2,solgrid3);
		if(pro->a[ia].gfree>0.1||pro->a[ia].gfree<-0.1)
			gen_atom_gfree_field(pro->a[ia].xyz,pro->a[ia].soltype, pro->a[ia].gfree, matter,solgrid1);
	}
}

#define WATEXPCOEF      0.50
#define WATMAXGRIDP    20

float den_dis(float dis2)
{
	return exp(-WATEXPCOEF*dis2);
}
int gen_atomden_grid(float xyz[3],float R,int p[],float per[])
{
	int i,j,k;
	int n;
	float dis2;
	float R2;
	int g,gi[3];
	float gxyz[3];
	int pn;
	float den[WATMAXGRIDP];
	float tden;

	R2=R*R;	
	n=(int)(R/GRIDSPACESIZE);
	pro2grid(xyz, gi);
	pn=0;
	tden=0.0;
	for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
		for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
			for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
				g=getgindex(i,j,k);
				grid2xyz(i,j,k,gxyz);
				dis2=distance2(xyz,gxyz);
				if(dis2<R2){
					den[pn]=den_dis(dis2);
					tden+=den[pn];
					p[pn]=g;
					pn++;
				}
			}
	for(n=0;n<pn;n++){
		per[n]=den[n]/tden;
	}

	return pn;
}
void gen_watden_grid(SURFWAT *wat, fftwf_complex *watgrid)
{
	int iw;
	float k[WATMAXGRIDP];
	int p[WATMAXGRIDP];
	int pn;
	int i;

	for(iw=0;iw<wat->wn;iw++)
	{
		pn=gen_atomden_grid(wat->w[iw].xyz,SURFWATRAD,p,k);
		for(i=0;i<pn;i++){
			watgrid[p[i]][1]+=(wat->w[iw].ene*k[i]);
		}
	}
}
void gen_proden_grid(STRUCTURE *pro, fftwf_complex *watgrid)
{
      int ia;
      int i,j,k;
      int n;
      float dis2;
      int g,gi[3];
      float gxyz[3];
	float R;
	float CORE=1.25;
	float FIRSTSHELLRAD=1.4;

        for(ia=0;ia<pro->atmN;ia++)
        {
		if(pro->a[ia].surf_core==SURFATOM) R=pro->a[ia].sol_rad;
		else R=pro->a[ia].sol_rad*CORE;
                n=(int)(R/GRIDSPACESIZE);
                pro2grid(pro->a[ia].xyz, gi);
                for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
                        for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
                                for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
                                        g=getgindex(i,j,k);
                                        grid2xyz(i,j,k,gxyz);
                                        dis2=distance2(pro->a[ia].xyz,gxyz);
                                        if(dis2<(R*R)){
                                                if(dis2<FIRSTSHELLRAD*FIRSTSHELLRAD&&watgrid[g][0]<1.0) watgrid[g][0]=1.0;
								else if(watgrid[g][0]<0.5)watgrid[g][0]=0.5;
                                        	    }
                                           }
        }
}
void gen_water_grid(STRUCTURE *pro, SURFWAT *wat, fftwf_complex *watgrid)
{
	gen_proden_grid(pro, watgrid);	
	gen_watden_grid(wat, watgrid);
}
void gen_piele_grid(STRUCTURE *pro, int *matter,fftwf_complex *indgrid, float shellR)
{
	int ia;
	int i,j,k;
	int n;
	float dis2;
	int g,gi[3];
	float gxyz[3];

	n=(int)(shellR/GRIDSPACESIZE);
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].piele<0.0001) continue;
		pro2grid(pro->a[ia].xyz, gi);
		for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
			for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
				for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
					g=getgindex(i,j,k);
					if(matter[g]!=0)
						continue;		
					grid2xyz(i,j,k,gxyz);
					dis2=distance2(pro->a[ia].xyz,gxyz);
					if(dis2<(shellR*shellR)){
						indgrid[g][1]+=(-pro->a[ia].piele*INDUCESTRENGTH);
					}
				}
	}
}

void gen_cation_grid(STRUCTURE *pro, fftwf_complex *indgrid)
{
	int ia;
	int g[MAXNEIGHB];
	int ii;
	float coe[MAXNEIGHB];
	int nn;
	float coet;

	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].cation<0.0001&&pro->a[ia].piele<0.0001) continue;
		nn=getnearest(pro->a[ia].xyz,g,coe);
		coet=0.0;
		for(ii=0;ii<nn;ii++) coet+=coe[ii];
		if(pro->a[ia].cation>=0.0001){
			for(ii=0;ii<nn;ii++){
				indgrid[g[ii]][0]+=pro->a[ia].cation*coe[ii]/coet;
			}
		}
		else{
			for(ii=0;ii<nn;ii++){
				indgrid[g[ii]][0]+=pro->a[ia].piele*coe[ii]/coet;
			}
		}
	}
}
void gen_induction_grid(STRUCTURE *pro, int *matter,fftwf_complex *indgrid)
{
	gen_piele_grid(pro, matter, indgrid,INDUCESHELL);	
	gen_cation_grid(pro,indgrid);
}

#define HBNEGATOM  1
#define HBPOSATOM  2
#define HBACCATOM  3
#define HBDONATOM  4
#define HBBBNATOM  5
#define HBBBOATOM  6
#define BBNHBSTRENGHT  0.15
#define BBOHBSTRENGHT  0.25

void gen_HBfield_grid(STRUCTURE *pro, int *matter,fftwf_complex *hbdgrid, float shellR, int type)
{
	int ia;
	int i,j,k;
	int n;
	float dis2;
	int g,gi[3];
	float gxyz[3];
	float v;

	n=(int)(shellR/GRIDSPACESIZE);
	for(ia=0;ia<pro->atmN;ia++){
		v=0.0;
		if(type==HBBBNATOM&&pro->a[ia].bbNO=='N') v=BBNHBSTRENGHT;
		else if(type==HBBBOATOM&&pro->a[ia].bbNO=='O') v=-BBOHBSTRENGHT;
		if(fabs(v)<0.0001) continue;
		pro2grid(pro->a[ia].xyz, gi);
		for(i=BSMAX(0,gi[0]-n);i<=BSMIN(GRIDSIZE[0]-1,gi[0]+n+1);i++)
			for(j=BSMAX(0,gi[1]-n);j<=BSMIN(GRIDSIZE[1]-1,gi[1]+n+1);j++)
				for(k=BSMAX(0,gi[2]-n);k<=BSMIN(GRIDSIZE[2]-1,gi[2]+n+1);k++){
					g=getgindex(i,j,k);
					if(matter[g]!=0)
						continue;		
					grid2xyz(i,j,k,gxyz);
					dis2=distance2(pro->a[ia].xyz,gxyz);
					if(dis2<(shellR*shellR)){
						hbdgrid[g][1]+=(v*HBONDSTRENGTH);
					}
				}
	}
}

void gen_HB_grid(STRUCTURE *pro, fftwf_complex *hbdgrid, int type)
{
	int ia;
	int g[MAXNEIGHB];
	int ii;
	float coe[MAXNEIGHB];
	int nn;
	float coet;
	float v;

	for(ia=0;ia<pro->atmN;ia++){
		v=0.0;
		if(type==HBACCATOM) v=-pro->a[ia].oxy;
		else if(type==HBDONATOM) v=pro->a[ia].hyd;
		if(fabs(v)<0.0001) continue;
		nn=getnearest(pro->a[ia].xyz,g,coe);
		coet=0.0;
		for(ii=0;ii<nn;ii++) coet+=coe[ii];
		for(ii=0;ii<nn;ii++){
			hbdgrid[g[ii]][0]+=v*coe[ii]/coet;
		}
	}
}

void gen_hbond_grid(STRUCTURE *pro, int *matter,fftwf_complex *Nhbgrid,fftwf_complex *Ohbgrid)
{
	gen_HBfield_grid(pro, matter, Nhbgrid,HBONDSHELL,HBBBNATOM);	
	gen_HBfield_grid(pro, matter, Ohbgrid,HBONDSHELL,HBBBOATOM);	
	gen_HB_grid(pro,Nhbgrid,HBACCATOM);
	gen_HB_grid(pro,Ohbgrid,HBDONATOM);
}

void print_grid(STRUCTURE *pro,fftwf_complex *grid,char *fn,float min1,float max1,float min2,float max2)
{
	int i,j,k,n=0;
	FILE *gridf;
	int g;
	float x[3];
	int ai,atno=1;
	
	gridf=fopen(fn,"w");
	
	for(ai=0;ai<pro->atmN;ai++){
		fprintf(gridf,"ATOM %6d  %s  %8.3f%8.3f%8.3f\n",
			   atno,pro->a[ai].name,
			   pro->a[ai].xyz[0],pro->a[ai].xyz[1],pro->a[ai].xyz[2]);
		atno++;
	}
	for(i=0;i<GRIDSIZE[0];i+=1)
		for(j=0;j<GRIDSIZE[1];j+=1){
			for(k=0;k<GRIDSIZE[2];k+=1){
				g=getgindex(i,j,k);
				grid2xyz(i,j,k, x);
				if(fabs(grid[g][1])>min1&&fabs(grid[g][1])<max1){
					n++;
					fprintf(gridf,"HETATM%5d  O   HOH B%4d   %8.3f%8.3f%8.3f  1.00%6.2f\n",n,n,x[0],x[1],x[2],grid[g][1]);
				}
			}
		}
	for(i=0;i<GRIDSIZE[0];i+=1)
		for(j=0;j<GRIDSIZE[1];j+=1){
			for(k=0;k<GRIDSIZE[2];k+=1){
				g=getgindex(i,j,k);
				grid2xyz(i,j,k, x);
				if(fabs(grid[g][0])>min2&&fabs(grid[g][0])<max2){
					n++;
					fprintf(gridf,"HETATM%5d  C   ALA A%4d   %8.3f%8.3f%8.3f  1.00%6.2f\n",n,n,x[0],x[1],x[2],grid[g][0]);
				}
			}
		}
		
	fclose(gridf);
}

void gen_receptor_grid(STRUCTURE *pro,SURFWAT *wat, int wmod)
{
	init_rgrid(Rrepgrid, GRID3SIZE);
	init_grid(Rvdwgrid, GRID3SIZE);
	init_grid(Relegrid, GRID3SIZE);
	init_grid(Rsolgrid1, GRID3SIZE);
	init_grid(Rsolgrid2, GRID3SIZE);
	init_grid(Rsolgrid3, GRID3SIZE);
	if(wmod==WATMODE) init_grid(Rwatgrid, GRID3SIZE);
	init_grid(Rindgrid, GRID3SIZE);
	init_grid(RNhbgrid, GRID3SIZE);
	init_grid(ROhbgrid, GRID3SIZE);
	
	init_matter(mattergrid, GRID3SIZE);
	gen_matter_region(pro, mattergrid);
	
	gen_rep_grid(pro, Rrepgrid);
	gen_matter_grid(pro, Rvdwgrid,Relegrid,Rsolgrid1,Rsolgrid2,Rsolgrid3);
	gen_field_grid(pro, mattergrid,Rvdwgrid,Relegrid,Rsolgrid1,Rsolgrid2,Rsolgrid3);
	if(wmod==WATMODE) gen_water_grid(pro, wat, Rwatgrid);
	gen_induction_grid(pro, mattergrid,Rindgrid);
	gen_hbond_grid(pro, mattergrid,RNhbgrid,ROhbgrid);
}

void gen_ligand_grid(STRUCTURE *pro,SURFWAT *wat, int wmod)
{
	init_rgrid(Lrepgrid, GRID3SIZE);
	init_grid(Lvdwgrid, GRID3SIZE);
	init_grid(Lelegrid, GRID3SIZE);
	init_grid(Lsolgrid1, GRID3SIZE);
	init_grid(Lsolgrid2, GRID3SIZE);
	init_grid(Lsolgrid3, GRID3SIZE);
	if(wmod==WATMODE) init_grid(Lwatgrid, GRID3SIZE);
	init_grid(Lindgrid, GRID3SIZE);
	init_grid(LNhbgrid, GRID3SIZE);
	init_grid(LOhbgrid, GRID3SIZE);
	
	init_matter(mattergrid, GRID3SIZE);
	gen_matter_region(pro, mattergrid);
	
	gen_rep_grid(pro, Lrepgrid);
	gen_matter_grid(pro, Lvdwgrid,Lelegrid,Lsolgrid1,Lsolgrid2,Lsolgrid3);
	gen_field_grid(pro, mattergrid,Lvdwgrid,Lelegrid,Lsolgrid1,Lsolgrid2,Lsolgrid3);
	if(wmod==WATMODE) gen_water_grid(pro, wat, Lwatgrid);
	gen_induction_grid(pro, mattergrid,Lindgrid);
	gen_hbond_grid(pro, mattergrid,LNhbgrid,LOhbgrid);
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

static int prime[164]={
	11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,
	109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,
	227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,
	347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,
	461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,
	599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,
	727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,
	859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997
};

int ps2GS(float ps)
{
	int n;
	int i;
	int isok=0;
	
	n=(int)((ps+AVERAGEVDWRAD)*2/GRIDSPACESIZE)+1;
	
	while(!isok){
		isok=1;
		if(n%4!=0){
			n++;
			isok=0;
			continue;
		}
		for(i=0;i<164;i++){
			if(prime[i]>n) break;
			if(n%prime[i]==0){
				isok=0;
				n=n+1;
				break;
			}
		}
	}
	
	return n;	
}

int decide_gridsize(STRUCTURE *p1,STRUCTURE *p2)
{
	float Rad1,Rad2,span1[3],span2[3];
	int GS1[3],GS2[3];
	int i;

	Rad1=proRad(p1,span1);
	Rad2=proRad(p2,span2);

	for(i=0;i<3;i++){
		GS1[i]=ps2GS(span1[i]+Rad2);
		GS2[i]=ps2GS(span2[i]+Rad1);	
	}
	if(GS1[0]*GS1[1]*GS1[2]<=GS2[0]*GS2[1]*GS2[2]){
		for(i=0;i<3;i++) GRIDSIZE[i]=GS1[i];
		GRID3SIZE=GS1[0]*GS1[1]*GS1[2];
		return 0;
	}
	else{
		for(i=0;i<3;i++) GRIDSIZE[i]=GS2[i];
		GRID3SIZE=GS2[0]*GS2[1]*GS2[2];
		return 1;
	}
}

