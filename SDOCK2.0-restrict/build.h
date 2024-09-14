/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * build.h : build models according to the docking records and write the model structure file.
 ***************************************************************************/

#include "structure.h"

void translateligand(STRUCTURE *pro, float x, float y, float z);
void build_model(STRUCTURE *receptor, STRUCTURE* ligand, int t,float m[3][3], float x, float y, float z,float shift[3]);
void write_complex(STRUCTURE *receptor, STRUCTURE* ligand,char *fn,int cluster,int member,float score,char *rname, char *lname,int t,char *recordf,float shift[3], int liandOnly);
