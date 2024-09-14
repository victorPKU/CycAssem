/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * watgrid.h: Detect the solvent geometry and chemical environment on the input protein surface by 3D grids.
 ***************************************************************************/

#ifndef GRID_H_
#define GRID_H_

#include "structure.h"

#define GRIDSPACESIZE 0.3
#define WATERVDWRAD  1.40
#define REPCOREVDWR 1.10
#define REPSURFVDWR 1.00
#define AVERAGEVDWRAD  1.80

#define WATRING 1.4
#define WATPRISM  2.2
#define WAT10D 3.0
#define WATCAGE  3.8
#define WATBULK  4.6

#define WATRINGGEO   5
#define WATPRISMGEO  4
#define WAT10DGEO   3
#define WATCAGEGEO  2
#define WATSURFGEO  1
#define WATBULKGEO  0

#define WATBULKENE  0.0
#define WATSHELLTHICK 2.80
#define WATWATDIS (2.8*2.8)

#define WATRINGHBN   1.65
#define WATPRISMHBN  2.2
#define WAT10DHBN   2.65
#define WATCAGEHBN  3.0
#define WATSURFHBN  3.25
#define WATBULKHBN  0

#define BSMAX(a,b)  (((a)>(b))?(a):(b))
#define BSMIN(a,b)  (((a)<(b))?(a):(b))

#define WATHB_LAM  10.0
#define WATHB_K    10.0
#define WATHB_C     1.0

int GRIDSIZE[3];
int GRID3SIZE;

void decide_gridsize(STRUCTURE *p1);
void init_mattergrid(int *grid, int N);
void init_HBgrid(float *grid, int N);
void init_watenegrid(float *grid, int N);
void init_pocketgrid(int *mattergrid, int *grid, int N);

void gen_mattergrid(STRUCTURE *pro, int *matbulk);
void gen_HBgrid(STRUCTURE *pro, float *HBgrid);
void gen_pocketgrid(STRUCTURE *pro,int *pocketgrid, int *mattergrid);

void gen_watenegrid(float *watenegrid, int *mattergrid, int *pocketgrid, float *HBgrid);
int enegrid2watpositin(int *mattergrid, int *pocketgrid,float *enegrid, float *HBgrid,float enewatxyz[][4], int start, int type);

void printwatxyz(char *inputnam,char *gfnam,int wn,float watxyz[][4]);

#endif /* GRID_H_ */
