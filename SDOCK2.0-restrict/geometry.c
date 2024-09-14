/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * geometry.c : distance calculation.
 ***************************************************************************/
#include "geometry.h"
#include <math.h>

float distance2(float a[3], float b[])
{
        return (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]);
}
float distance(float a[3], float b[3])
{
        return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]));
}
void cpx(float xi[3],float xo[3])
{
        int i;

        for(i=0;i<3;i++)
                xo[i]=xi[i];
}
void rotatm(float xi[3],float xo[3],float rm[3][3])
{
	int i,j;
	for(i=0;i<3;i++){
		xo[i]=0.0;
		for(j=0;j<3;j++)
			xo[i]+=xi[j]*rm[i][j];
	}
}

