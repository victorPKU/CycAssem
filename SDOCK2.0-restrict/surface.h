/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * surface.h : Calculate solvent accessible surface area for every protein atom.
 ***************************************************************************/
#define TEST_NSC 0

#define TEST_ARC 0
#define TEST_DOD 0
#define TEST_CUBE 0

#define UNSP_ICO_DOD      9
#define UNSP_ICO_ARC     10

#define PI 3.14159265358979323846
#define FOURPI		(4.*PI)
#define TORAD(A)     ((A)*0.017453293)
#define DP_TOL		0.001
#define MAXIMUM(a,b) (((a) > (b)) ? (a) : (b) )

float surfarea(int an, float ax[][3],float *rad,float *asurf,int densit, float solvR);
