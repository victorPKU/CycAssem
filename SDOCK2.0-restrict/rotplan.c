/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * rotplan.c : read the rotation sampling file and translate the quaternion to rotational matrix.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "rotplan.h"
#include "mem.h"

#define ROTFLINELEN   100

int readrotplan(char *fn)
{
	FILE *pf;
	int i,j;
	char line[ROTFLINELEN];

	if((pf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open the rotation sampling file %s!\n",fn);
		exit(0) ;
	}
	i=0;
	while(fgets(line,ROTFLINELEN,pf)){
		i++;
	}
	snew(ROTPLAN,i);
	i=0;
	rewind(pf);
	while(fgets(line,ROTFLINELEN,pf)){
		for(j=0;j<4;j++)
			ROTPLAN[i][j]=atof(line+15*j);
		/*printf("%5d%16.8f%16.8f%16.8f%16.8f\n",i,ROTPLAN[i][0],ROTPLAN[i][1],ROTPLAN[i][2],ROTPLAN[i][3]);*/
		i++;
	}
	fclose(pf);	
	return i;
}
void quat2matrix(float quat[4], float m[3][3])
{
	float a=quat[0],b=quat[1],c=quat[2],d=quat[3];

	m[0][0]=a*a+b*b-c*c-d*d;
	m[0][1]=2*b*c-2*a*d;
	m[0][2]=2*b*d+2*a*c;
	m[1][0]=2*b*c+2*a*d;
	m[1][1]=a*a-b*b+c*c-d*d;
	m[1][2]=2*c*d-2*a*b;
	m[2][0]=2*b*d-2*a*c;
	m[2][1]=2*c*d+2*a*b;
	m[2][2]=a*a-b*b-c*c+d*d;
}
void free_rotplan()
{
	sfree(ROTPLAN);
}
