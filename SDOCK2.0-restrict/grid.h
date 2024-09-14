/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * grid.h: Map the proteins to interaction grids.
 ***************************************************************************/

#include <fftw3.h>
#include "structure.h"

#define GRIDSPACESIZE 1.2

#define AVERAGEVDWRAD  1.80
#define REPCOREVDWR 1.2247449
#define REPSURFVDWR 0.750

#define ATTRSHORT  0.5808
#define ATTRLONG   0.1203
#define ELE_1   2.713
#define ELE_2   1.092
#define ELE_3   0.571
#define SOL1SHORT_3p5  (0.005556)
#define SOL1LONG_3p5   (0.002140)
#define SOL2SHORT_3p5  (0.003955)
#define SOL2LONG_3p5   (0.001361)

#define SOL1SHORT_6p0  (0.005703)
#define SOL1LONG_6p0   (0.002680)
#define SOL2SHORT_6p0  (0.004099)
#define SOL2LONG_6p0   (0.001851)

#define VDWFIRSTSHELL 1.4
#define VDWSECONDSHELL 1.8
#define ELE1SHELL 5.04
#define ELE2SHELL 6.48
#define ELE3SHELL 7.92
#define SOL1FIRSTSHELL 4.76
#define SOL1SECONDSHELL 6.12
#define SOL2FIRSTSHELL 5.6
#define SOL2SECONDSHELL 7.2
#define INDUCESHELL 4.50
#define INDUCESTRENGTH 3.0
#define HBONDSHELL 4.00
#define HBONDSTRENGTH 3.0

#define SURFWATRAD 1.4

#define WEIPARAN  8

float REPWEIGHT;
float ELEWEIGHT;
float SOLWEIGHT;
float WATWEIGHT;
float INDWEIGHT;
float HBNWEIGHT;
float HBOWEIGHT;

#define BSMAX(a,b)  (((a)>(b))?(a):(b))
#define BSMIN(a,b)  (((a)<(b))?(a):(b))

int GRIDSIZE[3];
int GRID3SIZE;

fftwf_complex	*Rrepgrid2,*Relegrid,*Rvdwgrid,*Rsolgrid1,*Rsolgrid2,*Rsolgrid3,*Rwatgrid,*Rindgrid,*RNhbgrid,*ROhbgrid;
fftwf_complex	*Lrepgrid2,*Lelegrid,*Lvdwgrid,*Lsolgrid1,*Lsolgrid2,*Lsolgrid3,*Lwatgrid,*Lindgrid,*LNhbgrid,*LOhbgrid;
float	*Rrepgrid,*Lrepgrid;
int *mattergrid;
float  *repscore, *vdwscore, *elescore, *solscore,*watscore,*indscore,*Nhbscore,*Ohbscore,*totalscore;
fftwf_complex  *ffttmp;

int getgindex(int gx,int gy,int gz);
int decide_gridsize(STRUCTURE *p1,STRUCTURE *p2);
void gen_receptor_grid(STRUCTURE *pro,SURFWAT *wat, int wmod);
void gen_ligand_grid(STRUCTURE *pro,SURFWAT *wat, int wmod);

