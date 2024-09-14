/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * protein.c : read force field parameter file and assign atom properties
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "protein.h"
#include "surface.h"
#include "geometry.h"
#include "mem.h"


void get_reslink(PROTEIN *pro)
{
	int ir;
	int ic;
	
	for(ic=0;ic<pro->chN;ic++){
		for(ir=0;ir<pro->c[ic].resN;ir++){
			if(ir==0){
				pro->c[ic].r[ir].pre_link=-1;
				continue;
			}
			if(ir==pro->c[ic].resN-1){
				pro->c[ic].r[ir].next_link=-1;
			}
			if(distance2(pro->c[ic].r[ir].a[0].xyz,pro->c[ic].r[ir-1].a[0].xyz)<=BLDIS2){
				pro->c[ic].r[ir].pre_link=pro->c[ic].r-pro->r+ir-1;
				pro->c[ic].r[ir-1].next_link=pro->c[ic].r-pro->r+ir;
			}
			else{
				pro->c[ic].r[ir].pre_link=-1;
				pro->c[ic].r[ir-1].next_link=-1;
			}
		}
	}
}
void centerpro(PROTEIN *pro)
{
	int ia;
	int i;
	float cx[3];

	for(i=0;i<3;i++)
		cx[i]=0.0;
	for(ia=0;ia<pro->atmN;ia++){
		for(i=0;i<3;i++)
			cx[i]+=pro->a[ia].xyz[i];
	}
	for(i=0;i<3;i++)
		cx[i]/=pro->atmN;
	for(ia=0;ia<pro->atmN;ia++){
		for(i=0;i<3;i++)
			pro->a[ia].xyzc[i]=pro->a[ia].xyz[i]-cx[i];
	}
}
void fixcenterpro(PROTEIN *pro,float cx[3])
{
	int ia;
	int i;

	for(ia=0;ia<pro->atmN;ia++){
		for(i=0;i<3;i++)
			pro->a[ia].xyzc[i]=pro->a[ia].xyz[i]-cx[i];
	}
}

#define SURFATOM  1
#define COREATOM  2

float cal_surface(PROTEIN *pro,float water_rad)
{
	float (*x)[3];
	float *s;
	float *r;
	int ia,i;
	float surf;
	int n;
	
	n=pro->atmN;
	
	snew(x,n);
	snew(s,n);
	snew(r,n);
	
	for(ia=0;ia<n;ia++){
		for(i=0;i<3;i++)
			x[ia][i]=pro->a[ia].xyzc[i];
		r[ia]=pro->a[ia].sol_rad;
	}
	
	surf=surfarea(n, x, r,s,SURF_CAL_DENSIT, water_rad);
	
	sfree(x);
	sfree(r);
	
	for(ia=0;ia<n;ia++){
		if(s[ia]>1.0)
			pro->a[ia].surf_core=SURFATOM;
		else
			pro->a[ia].surf_core=COREATOM;
	}
	sfree(s);
	return surf;
}
#define ATMFILELINE  90
void read_atm(char *atmparaf)
{
	FILE *atmf;
	char line[ATMFILELINE];
	if((atmf=fopen(atmparaf,"r"))==NULL){
		printf("ERROR: Can not open force field parameter file %s.\n", atmparaf);
		exit(0);
	}
	while(fgets(line,ATMFILELINE,atmf)){
		if(line[0]=='%')
			continue;
		strncpy(atdata[atdataN].resname,line,3);
		strncpy(atdata[atdataN].atname,line+4,3);
		atdata[atdataN].vdw_eps=sqrt(-atof(line+12));
		atdata[atdataN].vdw_rad=atof(line+19);
		atdata[atdataN].charge=atof(line+28);
		atdata[atdataN].sol_rad=atof(line+39);
		atdata[atdataN].HBstrength=atof(line+44);
		atdata[atdataN].soltype=atoi(line+48);
		atdata[atdataN].vol=atof(line+50);
		atdata[atdataN].gfree=atof(line+55);
		atdata[atdataN].piele=0.0;
		atdata[atdataN].cation=0.0;
		atdata[atdataN].hyd=0.0;
		atdata[atdataN].oxy=0.0;
		if(line[64]=='+') atdata[atdataN].cation=atof(line+66);
		if(line[64]=='*') atdata[atdataN].piele=atof(line+66);
		if(line[65]=='>'||line[65]=='x') atdata[atdataN].hyd=atof(line+66);
		if(line[65]=='o'||line[65]=='x') atdata[atdataN].oxy=atof(line+66);
		atdataN++;
	}
	fclose(atmf);
}
int lookupATMtab(char *aname,char *rname)
{
	int i;
	
	for(i=0;i<atdataN;i++){
		if(!strcmp(aname,atdata[i].atname)
			&&(!strcmp(rname,atdata[i].resname)||!strcmp(atdata[i].resname, "ARB")))
			return i;
	}
	return -1;
}

void get_atom_para(PROTEIN *pro, char *fn)
{
	int ia;
	int index;
	
	read_atm(fn);
	for(ia=0;ia<pro->atmN;ia++){
		index=lookupATMtab(pro->a[ia].abbrname,pro->a[ia].name+4);
		if(index==-1){
			printf("ERROR: Can not find atom type '%s' in the force field parameter file %s!\n",pro->a[ia].name,fn);
			exit(0);
		}
		pro->a[ia].vdw_eps=atdata[index].vdw_eps;
		pro->a[ia].vdw_rad=atdata[index].vdw_rad;
		pro->a[ia].charge=atdata[index].charge;
		pro->a[ia].sol_rad=atdata[index].sol_rad;
		pro->a[ia].HBstrength=atdata[index].HBstrength;
		pro->a[ia].soltype=atdata[index].soltype;
		pro->a[ia].volume=atdata[index].vol;
		pro->a[ia].gfree=atdata[index].gfree;
		pro->a[ia].piele=atdata[index].piele;
		pro->a[ia].cation=atdata[index].cation;
		pro->a[ia].hyd=atdata[index].hyd;
		pro->a[ia].oxy=atdata[index].oxy;
		if(!strcmp(pro->a[ia].abbrname,"N  ")&&strcmp(pro->a[ia].name+4,"PRO")) pro->a[ia].bbNO='N';
		else if(!strcmp(pro->a[ia].abbrname,"O  ")) pro->a[ia].bbNO='O';
		else pro->a[ia].bbNO=' ';
	}
}

void free_pro(PROTEIN *pro)
{
	sfree(pro->a);
	sfree(pro->r);
	sfree(pro->c);
}
