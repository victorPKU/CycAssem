/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * geometry.h : distance calculation.
 ***************************************************************************/
#ifndef _geometry_h
#define _geometry_h

float distance2(float a[3],float b[3]);
float distance(float a[3],float b[3]);
void cpx(float xi[3],float xo[3]);
void rotatm(float xyz0[3],float xyz[3],float rm[3][3]);

#endif
