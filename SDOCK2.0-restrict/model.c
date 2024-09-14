/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * model.c : build models according to the docking records and write the model structure file.
 ***************************************************************************/

#include "build.h"
#include <stdio.h>
#include <stdlib.h>

void translateligand(STRUCTURE *pro, float x, float y, float z)
{
	int ati;
	
	for(ati=0;ati<pro->atmN;ati++){
		pro->a[ati].xyz[0]-=x;
		pro->a[ati].xyz[1]-=y;
		pro->a[ati].xyz[2]-=z;
	}
}
void build_model(STRUCTURE *receptor, STRUCTURE* ligand, int t, float m[3][3], float x, float y, float z,float shift[3])
{
	if(t==0){
		fixpro(receptor);
		rotligand(ligand,m);
		translateligand(ligand,x,y,z);
	}
	else{
		fixpro(ligand);
		rotligand(receptor,m);
		translateligand(receptor,x,y,z);
	}
	translateligand(receptor,shift[0],shift[1],shift[2]);
	translateligand(ligand,shift[0],shift[1],shift[2]);
}
void write_complex(STRUCTURE *receptor, STRUCTURE* ligand,char *fn,int cluster,int member,float score, char *rname, char *lname,int t,char *recordf,float shift[3], int liandOnly)
{
	FILE *pf;
	int ai,atno=1;

	if((pf=fopen(fn,"w"))==NULL){
		printf("ERROR: Can not create model structure file %s\n",fn);
		exit(0);
	}
	fprintf(pf,"HEADER    SDOCK model file\n");
	if(t==0){
		fprintf(pf,"REMARK    static %s\n",rname);
		fprintf(pf,"REMARK    mobile %s\n",lname);
	}
	else{
		fprintf(pf,"REMARK    mobile %s\n",rname);
		fprintf(pf,"REMARK    static %s\n",lname);
	}
	fprintf(pf,"REMARK    shift vector %8.3f%8.3f%8.3f\n",shift[0],shift[1],shift[2]);
	fprintf(pf,"REMARK    record file %s\n",recordf);
	fprintf(pf,"REMARK    cluster %5d : member %5d\n",cluster,member);
	fprintf(pf,"REMARK    score %8.3f\n",score);
	if(liandOnly==0){
		for(ai=0;ai<receptor->atmN;ai++){
			fprintf(pf,"ATOM %6d  %s   %8.3f%8.3f%8.3f\n",
				   atno,receptor->a[ai].name,
				   receptor->a[ai].xyz[0],receptor->a[ai].xyz[1],receptor->a[ai].xyz[2]);
			atno++;
		}
		fprintf(pf,"END\n");
	}
	else{
		fprintf(pf,"REMARK    LIGAND ONLY\n");
	}
	for(ai=0;ai<ligand->atmN;ai++){
		fprintf(pf,"ATOM %6d  %s   %8.3f%8.3f%8.3f\n",
			   atno,ligand->a[ai].name,
			   ligand->a[ai].xyz[0],ligand->a[ai].xyz[1],ligand->a[ai].xyz[2]);
		atno++;
	}
	fprintf(pf,"END\n");

	printf("==> The complex model %s is built!\n",fn);

	fclose(pf);
}
