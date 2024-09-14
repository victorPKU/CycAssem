/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * structure.h : read the input structure file and center, rotate the proteins.
 ***************************************************************************/

#ifndef _structure_h
#define _structure_h

typedef struct{
	char name[15];
	int isCa;
	float vdw_rad;
	float sq_vdw_eps;
	float charge;
	int soltype;
	float volume;
	float gfree;
	int surf_core;
	float sol_rad;
	float HBstrength;
	float piele;
	float cation;
	char bbNO;
	float hyd;
	float oxy;
	float fluc;
	float xyzc[3];
	float xyz[3];
}ATM;
typedef struct{
	int atmN;
	ATM *a;
	int surfatN;
	int CaN;
	float fluc0;
}STRUCTURE;

typedef struct{
	float xyzc[3];
	float xyz[3];
	float ene;
}WATATM;
typedef struct{
	int wn;
	WATATM *w;
}SURFWAT;

#define SURFATOM  1
#define COREATOM  2
#define WATMODE   1
#define NOWATMODE 0

#define PI 3.14159265358979323846

void read_structure(STRUCTURE *pro, char *fn);
void read_surfacewater(SURFWAT *wat, char *fn);
void fixpro(STRUCTURE *pro);
void fixwat(SURFWAT *wat);
void rotligand(STRUCTURE *pro,float m[3][3]);
void rotligwat(SURFWAT *wat,float rm[3][3]);
int get_ligandca(STRUCTURE * ligand, float Ca[][3]);
void rotligandCa(int CaN, float Ca0[][3], float Ca[][3],float rm[3][3]);


#endif
