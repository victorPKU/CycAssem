/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * rotplan.h : read the rotation sampling file and translate the quaternion to rotational matrix.
 ***************************************************************************/

float (*ROTPLAN)[4];
int ROTN;

int readrotplan(char *fn);
void quat2matrix(float quat[4], float m[3][3]);
void free_rotplan();
