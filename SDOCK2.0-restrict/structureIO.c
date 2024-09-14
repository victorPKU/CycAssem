/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * structureIO.c : read pdb file, and write the processed file..
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "protein.h"
#include "mem.h"

#define LINELEN 81

void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+8*i);
}
int IstheChn(char c,int n,char *chn)
{	
	int i;
	
	if(chn[0]=='@')
		return 1;
	for(i=0;i<n;i++){
		if(chn[i]=='_'&&c==' ')
			return 1;
		if(chn[i]==c)
			return 1;
	}
	return 0;
}
int IstheucAA(char c,char *het, int n, char name[][5])
{
	int i;
	
	for(i=0;i<n;i++){
		if(name[i][0]=='@'&&!strncmp(het,name[i]+1,3))
			return 1;
		if(c==name[i][0]&&!strncmp(het,name[i]+1,3))
			return 1;
	}
	return 0;
}
int Isthesmol(char c,char *het, int n, char name[][5])
{
	int i;
	
	for(i=0;i<n;i++){
		if(name[i][0]=='@'&&!strncmp(het,name[i]+1,3))
			return 1;
		if(c==name[i][0]&&!strncmp(het,name[i]+1,3))
			return 1;
	}
	return 0;
}

void readprotein(PROTEIN *newpro, char *pn, char *chn, int ucAAN, char ucAA[][5], int smolN, char smol[][5])
{
	FILE *pf;
	char line[LINELEN];
	int chi,chj;
	int atmi;
	int resi,resj;
	char id=0;
	int chN;
	char conform=' ';
	int MAXATMN=0;
	int MAXRESN;
	int MAXCHNN;
	int restyp=NMAA;
	
	chi=atmi=resi=0;
	
	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",pn);
		exit(0);
	}

	chN=strlen(chn);
	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)||!strncmp(line,"HETATM",6)){
			 MAXATMN++;
		}
	}
	rewind(pf);
	MAXRESN=MAXATMN;
	MAXCHNN=MAXRESN;
	
	snew(newpro->a, MAXATMN);
	snew(newpro->r, MAXRESN);
	snew(newpro->c, MAXCHNN);

	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ENDMDL",5)) break;
		restyp=-1;
		if(!strncmp(line,"ATOM",4)&&IstheChn(line[21],chN,chn)) restyp=NMAA;
		else if(!strncmp(line,"HETATM",6)&&IstheucAA(line[21],line+17, ucAAN, ucAA)) restyp=UCAA;
		else if(!strncmp(line,"HETATM",6)&&Isthesmol(line[21],line+17, smolN, smol)) restyp=SMOL;
		else continue;
		if(line[13]=='H'||line[12]=='H'){
			continue;
		}
		if(line[16]!=' '){
			if(conform==' ')
				conform=line[16];
			else if(line[16]!=conform)
				continue;
		}
		if(id!=line[21]){
			id=line[21];
			newpro->c[chi].id=line[21];
			newpro->c[chi].r=newpro->r+resi;
			chi++;
		}
		if(resi==0||strncmp(line+17,newpro->r[resi-1].name,9)){
			strncpy(newpro->r[resi].name,line+17,9);
			newpro->r[resi].name[9]='\0';
			strncpy(newpro->r[resi].abbrname,line+17,3);
			newpro->r[resi].abbrname[3]='\0';
			newpro->r[resi].ishet=restyp;
			newpro->r[resi].atch=newpro->c+chi-1;
			newpro->r[resi].a=newpro->a+atmi;
			newpro->r[resi].atmN=0;
			newpro->r[resi].CA=newpro->a+atmi;
			resi++;
		}
		newpro->r[resi-1].atmN++;
		newpro->a[atmi].atres=newpro->r+resi-1;
		getxyz(line+30,newpro->a[atmi].xyz);
		strncpy(newpro->a[atmi].name,line+13,7);
		newpro->a[atmi].name[7]='\0';
		strncpy(newpro->a[atmi].abbrname,line+13,3);
		newpro->a[atmi].abbrname[3]='\0';
		if(restyp!=NMAA&&line[12]!=' '){
			newpro->a[atmi].name[0]=line[12];
			newpro->a[atmi].name[1]=line[13];
			newpro->a[atmi].name[2]='+';
			newpro->a[atmi].abbrname[0]=line[12];
			newpro->a[atmi].abbrname[1]=line[13];
			newpro->a[atmi].abbrname[2]='+';
		}
		if(!strncmp("CA ",line+13,3)||!strncmp("CAT",line+13,3))
			newpro->r[resi-1].CA=newpro->a+atmi;
		atmi++;
	}
	newpro->atmN=atmi;
	newpro->chN=chi;
	newpro->resN=resi;

	if(chi==0){
		if(chn[0]=='@')
			printf("ERROR: No protein chains are found in the file %s!\n",pn);
		else 
			printf("ERROR: chains %s are not found in the file %s!\n",chn,pn);
		free_pro(newpro);
		exit(0);
	}
	srenew(newpro->a, newpro->atmN);
	srenew(newpro->r, newpro->resN);
	srenew(newpro->c, newpro->chN);
	
	for(chi=0;chi<(newpro->chN-1);chi++){
		newpro->c[chi].resN=newpro->c[chi+1].r-newpro->c[chi].r;
	}
	newpro->c[chi].resN=newpro->r-newpro->c[chi].r+newpro->resN;

	if(chn[0]!='@'){
		for(chi=0;chi<chN;chi++){
			for(chj=0;chj<newpro->chN;chj++){
				if(newpro->c[chj].id==chn[chi])
					break;
				else if(newpro->c[chj].id==' '&&chn[chi]=='_')
					break;
			}
			if(chj>=newpro->chN)
				printf("WARNING: the chain %c is not found in the pdb file %s!\n",chn[chi],pn);
		}
	}

	for(resi=0;resi<ucAAN;resi++){
		for(resj=0;resj<newpro->resN;resj++){
			if(ucAA[resi][0]=='@'&&ucAA[resi][1]==newpro->r[resj].name[0]&&ucAA[resi][2]==newpro->r[resj].name[1]&&ucAA[resi][3]==newpro->r[resj].name[2])
				break;
			else if(ucAA[resi][0]==newpro->r[resj].name[4]&&ucAA[resi][1]==newpro->r[resj].name[0]&&ucAA[resi][2]==newpro->r[resj].name[1]&&ucAA[resi][3]==newpro->r[resj].name[2]){
				break;
			}
		}
		if(resj>=newpro->resN)
			printf("WARNING: the uncommon amino acid %s is not found in the pdb file %s!\n",ucAA[resi],pn);
	}
	
	for(resi=0;resi<smolN;resi++){
		for(resj=0;resj<newpro->resN;resj++){
			if(smol[resi][0]=='@'&&smol[resi][1]==newpro->r[resj].name[0]&&smol[resi][2]==newpro->r[resj].name[1]&&smol[resi][3]==newpro->r[resj].name[2])
				break;
			else if(smol[resi][0]==newpro->r[resj].name[4]&&smol[resi][1]==newpro->r[resj].name[0]&&smol[resi][2]==newpro->r[resj].name[1]&&smol[resi][3]==newpro->r[resj].name[2]){
				break;
			}
		}
		if(resj>=newpro->resN)
			printf("WARNING: the cofactor %s is not found in the pdb file %s!\n",smol[resi],pn);
	}

	fclose(pf);
}

void write_structure(PROTEIN *pro, char *fn,char *argv[], int arg_chain, int arg_atmfile, int arg_ucAA, int arg_smol)
{
	int ai,atno=1;
	FILE *pf;
	RES *res;
	
	if((pf=fopen(fn,"w"))==NULL){
		printf("ERROR: Can not create the file %s\n",fn);
		exit(0);
	}
	fprintf(pf,"HEADER    SDOCK input file\n");
	fprintf(pf,"REMARK    origin pdb file: %s\n", argv[1]);
	if(arg_chain==0)
		fprintf(pf,"REMARK    chains: <ALL>\n");
	else
		fprintf(pf,"REMARK    chains: %s\n", argv[arg_chain]);
	if(arg_ucAA==0)
		fprintf(pf,"REMARK    uncommon amino acids: no\n" );
	else
		fprintf(pf,"REMARK    uncommon amino acids: %s\n", argv[arg_ucAA]);
	if(arg_smol==0)
		fprintf(pf,"REMARK    cofactors: no\n");
	else
		fprintf(pf,"REMARK    cofactors: %s\n", argv[arg_smol]);
	if(arg_atmfile==0)
		fprintf(pf,"REMARK    force field parameter file: ATM\n");
	else
		fprintf(pf,"REMARK    force field parameter file: %s\n", argv[arg_atmfile]);
	

	for(ai=0;ai<pro->atmN;ai++){
		res=pro->a[ai].atres;
		fprintf(pf,"ATOM%7d  %s%s    %8.3f%8.3f%8.3f%3d %8.3f%8.3f\n",
			   atno,pro->a[ai].name,res->name+3,
			   pro->a[ai].xyzc[0],pro->a[ai].xyzc[1],pro->a[ai].xyzc[2],pro->a[ai].surf_core,pro->a[ai].fluc*100,pro->a[ai].charge);
		fprintf(pf,"FFPA%7d  %s%s    %8.3f%8.3f%8.3f%3d %8.3f%8.3f\n",
			   atno,pro->a[ai].name,res->name+3,
			   pro->a[ai].vdw_rad,pro->a[ai].sol_rad,pro->a[ai].volume,pro->a[ai].soltype,pro->a[ai].vdw_eps,pro->a[ai].gfree);
		fprintf(pf,"ADPA%7d  %s%s    %8.3f%8.3f%8.3f%3c %8.3f%8.3f\n",
			   atno,pro->a[ai].name,res->name+3,
			   pro->a[ai].HBstrength,pro->a[ai].piele,pro->a[ai].cation,pro->a[ai].bbNO,pro->a[ai].hyd,pro->a[ai].oxy);
		atno++;
	}
	fclose(pf);

	printf("\n==> The input file %s is successfully processed to %s now!\n", argv[1],fn);
}
