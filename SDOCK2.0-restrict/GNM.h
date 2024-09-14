/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * GNM.h : GNM model for atom fluctuation calculation.
 ***************************************************************************/

#ifndef GNM_H_
#define GNM_H_
#include "protein.h"

void genblocks(PROTEIN *pro);
void genmatrix(PROTEIN *pro, double *matrix);
void calfluctuation(double *matrix, PROTEIN *pro);

#endif /* GNM_H_ */
