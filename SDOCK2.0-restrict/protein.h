/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * protein.h : read pdb file, assign atom properties and write the processed file.
 ***************************************************************************/
#ifndef _protein_h
#define _protein_h

typedef struct{
	char name[8];
	char abbrname[4];
	void *atres;
	float vdw_rad;
	float vdw_eps;
	float charge;
	float sol_rad;
	float HBstrength;
	int soltype;
	float volume;
	float gfree;
	float piele;
	float cation;
	char bbNO;
	float hyd;
	float oxy;
	float fluc;
	int surf_core;
	float xyz[3];
	float xyzc[3];
}ATM;
#define NMAA 0
#define UCAA 1
#define SMOL 2
typedef struct{
	void *atch;
	char name[10];
	char abbrname[4];
	int ishet;
	int pre_link;
	int next_link;
      int blockN;
	int atmN;
	ATM *a;
      ATM *CA;
      int block_i[4];
      int block_at[4][10];
}RES;
typedef struct{
	int resN;
	RES *r;
	char id;
}CHN;
typedef struct{
	int chN;
	int resN;
	int atmN;
      int blockN;
      int ssbN;
	ATM *a;
	RES *r;
	CHN *c;
      int *ssb;
}PROTEIN;
#define WATERRAD  1.40
#define BLDIS2 4.0*4.0
#define SURF_CAL_DENSIT   400
#define PI 3.14159265358979323846
typedef struct{
	char resname[4];
	char atname[4];
	float vdw_eps;
	float vdw_rad;
	float charge;
	float HBstrength;
	float sol_rad;
	int soltype;
	float vol;
	float gfree;
	float piele;
	float cation;
	float hyd;
	float oxy;
}PROTEINAT;

int atdataN;
PROTEINAT atdata[800];

void readprotein(PROTEIN *newpro, char *pn, char *chn, int ucAAN, char ucAA[][5], int smolN, char smol[][5]);
void write_structure(PROTEIN *pro, char *fn,char *argv[], int arg_chain, int arg_atmfile, int arg_ucAA, int arg_smol);
void get_atom_para(PROTEIN *pro, char *fn);
float cal_surface(PROTEIN *pro,float water_rad);
void centerpro(PROTEIN *pro);
void fixcenterpro(PROTEIN *pro,float cx[3]);
void get_reslink(PROTEIN *pro);
void free_pro(PROTEIN *pro);

#endif
