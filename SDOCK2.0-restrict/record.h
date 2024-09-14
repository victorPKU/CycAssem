/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * record.h : save the docking results and sort them by docking score.
 ***************************************************************************/
#ifndef _record_h
#define _record_h

int  RESULTKEEPN;
typedef struct{
	int cluster;
	int num;
	float RMSD;	
	int roti;	
	float mx;
	float my;
	float mz;
	float score;
	float rep;
	float vdw;
	float ele;
	float sol;
	float wat;
	float ind;
	float hbN;
	float hbO;
}DOCKRECORD;

DOCKRECORD *dock_result;
DOCKRECORD *dock_resultf;
int *result_index;
struct LIST{
	int i;
	struct LIST *n;
};
struct LIST *order;
int count;
	
void init_record();
void arrange_record();
void save_dock(int roti, float trans_x,float trans_y,float trans_z,float score,float rep,float vdw,float ele, float sol, float wat, float ind, float hbN, float hbO);
void free_record();
int read_record(char *fn, int *t,int cluster, int number);

#endif

