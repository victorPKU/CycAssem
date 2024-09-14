/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * cluster.h : cluster the results by interface Calpha RMSD and output the docking results.
 ***************************************************************************/

#define MAXMEMBER  300
#define MAXCLUSTER  9999
typedef struct{
	int N;
	float rmsdcutoff;
	int id[MAXMEMBER];
	float Lrmsd[MAXMEMBER];
}CLUSTER;

int clusterN;
int max_clusterN;
CLUSTER cluster[MAXCLUSTER];
float LRMSDCUTOFF;

void cluster_result(int RCaN,float RCa[][3],int CaN,float Ca[][3], int mod);
void print_cluster(int mode,char *pA, char *pB, int t, char *fn,char *rfn, float cw, float ew, float sw, float ww, float iw, float nw, float ow);
