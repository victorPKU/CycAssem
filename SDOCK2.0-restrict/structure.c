/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * structure.c : read the input structure file and center, rotate the proteins.
 ***************************************************************************/

#include "structure.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mem.h"
#include "geometry.h"

void fixpro(STRUCTURE *pro)
{
	int ia;

	for(ia=0;ia<pro->atmN;ia++){
		cpx(pro->a[ia].xyzc, pro->a[ia].xyz);
	}
}
void fixwat(SURFWAT *wat)
{
	int iw;

	for(iw=0;iw<wat->wn;iw++){
		cpx(wat->w[iw].xyzc, wat->w[iw].xyz);
	}
}
void rotligand(STRUCTURE *pro,float rm[3][3])
{
	int ai;
	
	for(ai=0;ai<pro->atmN;ai++){
		rotatm(pro->a[ai].xyzc,pro->a[ai].xyz,rm);
	}
}
void rotligwat(SURFWAT *wat,float rm[3][3])
{
	int wi;
	
	for(wi=0;wi<wat->wn;wi++){
		rotatm(wat->w[wi].xyzc,wat->w[wi].xyz,rm);
	}
}

void rotligandCa(int CaN, float Ca0[][3], float Ca[][3],float rm[3][3])
{
	int ati;
	
	for(ati=0;ati< CaN;ati++){
		rotatm(Ca0[ati],Ca[ati],rm);
	}
}

#define LINELEN 100
void read_structure(STRUCTURE *p, char *fn)
{
	FILE *pf;
	char line[LINELEN];
	int an;
	int i;
	int surfat=0;
	int caat=0;
	
	if((pf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open structure file %s!\n",fn);
		exit(0) ;
	}
	an=0;
	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			an++;
		}
	}
			
	snew(p->a,an);
	
	rewind(pf);
	an=0;
	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			strncpy(p->a[an].name, line+13,14);
			p->a[an].isCa=0;
			if(!strncmp(p->a[an].name,"CA ",3)||!strncmp(p->a[an].name,"CAT",3)){
				p->a[an].isCa=1;
				caat++;
			}
			for(i=0;i<3;i++)
				p->a[an].xyzc[i]=atof(line+30+8*i);
			p->a[an].surf_core=atoi(line+54);
			if(p->a[an].surf_core==SURFATOM) surfat++;
			p->a[an].fluc=atof(line+58);
			p->a[an].charge=atof(line+66);
			p->a[an].cation=0.0;
			p->a[an].hyd=0.0;
			p->a[an].oxy=0.0;
			p->a[an].bbNO=' ';
			an++;
		}
		else if(!strncmp(line,"FFPA",4)){
			p->a[an-1].vdw_rad=atof(line+30);
			p->a[an-1].sol_rad=atof(line+38);
			p->a[an-1].volume=atof(line+46);
			p->a[an-1].soltype=atoi(line+54);
			p->a[an-1].sq_vdw_eps=atof(line+58);
			p->a[an-1].gfree=atof(line+66);
		}
		else if(!strncmp(line,"ADPA",4)){
			p->a[an-1].HBstrength=atof(line+30);
			p->a[an-1].piele=atof(line+38);
			p->a[an-1].cation=atof(line+46);
			if(p->a[an-1].surf_core==SURFATOM){
				p->a[an-1].bbNO=line[56];
				p->a[an-1].hyd=atof(line+58);
				p->a[an-1].oxy=atof(line+66);
			}
		}
	}
	fclose(pf);	
	
	p->atmN=an;
	p->surfatN=surfat;
	p->CaN=caat;
}
void read_surfacewater(SURFWAT *wat, char *fn)
{
	FILE *wf;
	char line[LINELEN];
	int an;
	int i;
	
	if((wf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open structure file %s!\n",fn);
		exit(0) ;
	}
	an=0;
	while(fgets(line,LINELEN,wf)){
		if(!strncmp(line,"HETATM",6)){
			an++;
		}
	}
			
	snew(wat->w,an);
	
	rewind(wf);
	an=0;
	while(fgets(line,LINELEN,wf)){
		if(!strncmp(line,"HETATM",6)){
			for(i=0;i<3;i++)
				wat->w[an].xyzc[i]=atof(line+30+8*i);
			wat->w[an].ene=atof(line+60);
			an++;
		}
	}
	fclose(wf);	
	
	wat->wn=an;
}
int get_ligandca(STRUCTURE * ligand, float Ca[][3])
{
	int ia,i;
	int n=0;

	for(ia=0;ia<ligand->atmN;ia++){
		if(ligand->a[ia].isCa==1){
			for(i=0;i<3;i++)
				Ca[n][i]=ligand->a[ia].xyzc[i];
			n++;
			if(n>=ligand->CaN) break;
		}
	}

	return n;
}
