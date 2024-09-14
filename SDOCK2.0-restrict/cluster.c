/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * cluster.c : cluster the results by interface Calpha RMSD and output the docking results.
 ***************************************************************************/

#include "cluster.h"
#include "structure.h"
#include "record.h"
#include "mem.h"
#include "rotplan.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void new_cluster(int i, int mod)
{
	if(clusterN>=max_clusterN)
		return;
	cluster[clusterN].N=1;
	cluster[clusterN].id[0]=i;
	cluster[clusterN].Lrmsd[0]=0.0;
	if(mod==1||clusterN==0) cluster[clusterN].rmsdcutoff=LRMSDCUTOFF;
	/*else cluster[clusterN].rmsdcutoff=LRMSDCUTOFF*exp(1.0*clusterN/max_clusterN);*/
	else cluster[clusterN].rmsdcutoff=LRMSDCUTOFF+log10(1.0*clusterN);
	clusterN++;
}

void add2cluster(int i, int cluster_i,float rmsd)
{
	int k=cluster[cluster_i].N;

	if(k>=MAXMEMBER)
		return;
	cluster[cluster_i].id[k]=i;
	cluster[cluster_i].Lrmsd[k]=rmsd;
	cluster[cluster_i].N++;
}

void gen_ligand(int CaN, float Ca0[][3], float Ca[][3],DOCKRECORD r)
{
	int ia;
	float rm[3][3];

	quat2matrix(ROTPLAN[r.roti], rm);
	rotligandCa(CaN, Ca0, Ca,rm);
	for(ia=0;ia<CaN;ia++){
		Ca[ia][0]-=r.mx;
		Ca[ia][1]-=r.my;
		Ca[ia][2]-=r.mz;
	}
}

#define INTERFCADIS2   (15*15)

int IsInterfCa(int RCaN,float RCa[][3], float LCa[3])
{
	int i;
	float d;

	for(i=0;i<RCaN;i++){
		d=(RCa[i][0]-LCa[0])*(RCa[i][0]-LCa[0])+(RCa[i][1]-LCa[1])*(RCa[i][1]-LCa[1])+(RCa[i][2]-LCa[2])*(RCa[i][2]-LCa[2]);
		if(d<INTERFCADIS2)
			return 1;		
	}
	return 0;
}
float cal_Lrmsd(int RCaN,float RCa[][3],int CaN,float Ca[][3],float A_Ca[][3],float B_Ca[][3],int iA, int iB,float rmsdcutoff)
{
	int Aindex,Bindex;
	int ia,i;
	float Lrmsd=0;
	float d;
	int an;
	float INTERFRMSDIS2=(rmsdcutoff*rmsdcutoff*1.5*1.5);

	Aindex=result_index[iA];
	Bindex=result_index[iB];
	gen_ligand(CaN,Ca,A_Ca,dock_result[Aindex]);
	gen_ligand(CaN,Ca,B_Ca,dock_result[Bindex]);

	an=0;
	for(ia=0;ia<CaN;ia++){
		if(IsInterfCa(RCaN,RCa,B_Ca[ia])==0)
			continue;
		d=0.0;
		for(i=0;i<3;i++)
			d+=(A_Ca[ia][i]-B_Ca[ia][i])*(A_Ca[ia][i]-B_Ca[ia][i]);
		if(d>INTERFRMSDIS2)
			return rmsdcutoff*1.5+0.1;
		Lrmsd+=d;
		an++;
	}
	if(an==0)
		Lrmsd=rmsdcutoff*1.5+0.1;
	else
		Lrmsd=sqrt(Lrmsd/an);
	
	return Lrmsd;
}


void cluster_result(int RCaN,float RCa[][3],int CaN,float Ca[][3], int mod)
{
	int i,j;
	float rmsd;
	float (*A_Ca)[3],(*B_Ca)[3];
	
	clusterN=0;
	new_cluster(0,mod);

	snew(A_Ca,CaN);
	snew(B_Ca,CaN);
	for(i=1;i<RESULTKEEPN;i++){
		for(j=0;j<clusterN;j++){
			rmsd=cal_Lrmsd(RCaN,RCa,CaN,Ca,A_Ca,B_Ca,i,cluster[j].id[0],cluster[j].rmsdcutoff);
			if(rmsd<cluster[j].rmsdcutoff){
				add2cluster(i,j,rmsd);
				break;
			}
		}
		if(j>=clusterN){
			new_cluster(i,mod);
			if(clusterN>=max_clusterN)
				break;
		}
	}
	sfree(A_Ca);
	sfree(B_Ca);
}

void print_cluster(int mode,char *pA, char *pB, int t,char *fn,char *rfn, float cw, float ew, float sw, float ww, float iw, float nw, float ow)
{
	FILE *clusterf;
	int i, j,k;

	if((clusterf=fopen(fn,"w"))==NULL)
		printf("ERROR: Can't create the sdock result file %s!\n",fn);

	fprintf(clusterf,"%%SDOCK result for the docking work:");
	if(mode==NOWATMODE) fprintf(clusterf," (without explict water)\n");
	else fprintf(clusterf," (with explict water)\n");
	if(t==0){
		fprintf(clusterf,"   Protein A: static %s\n",pA);
		fprintf(clusterf,"   Protein B: mobile %s\n",pB);
	}
	else{
		fprintf(clusterf,"   Protein A: mobile %s\n",pA);
		fprintf(clusterf,"   Protein B: static %s\n",pB);
	}
	fprintf(clusterf,"%%rotation sampling file: %s\n",rfn);
	if(mode==WATMODE) fprintf(clusterf,"%%weights of score terms (collision, electrostatic, desolvation, water, induction, bbNHB, bbOHB):%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",cw,ew,sw,ww,iw,nw,ow);
	else fprintf(clusterf,"%%weights of score terms (collision, electrostatic, desolvation, water, induction, bbNHB, bbOHB):%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",cw,ew,sw,iw,nw,ow);
	fprintf(clusterf,"%%cluster file: %s\n",fn);
	fprintf(clusterf,"%%cluster number: %d\n",clusterN);
	fprintf(clusterf,"%%cluster radius: %5.3f\n",LRMSDCUTOFF);
	if(mode==WATMODE) fprintf(clusterf,"%%cluster no MIRMSD  |    rota     x       y      z    |     score:collision  vdw     ele     sol     wat     ind     hbN     hbO\n");
	else fprintf(clusterf,"%%cluster no MIRMSD  |    rota     x       y      z    |     score:collision  vdw     ele     sol     ind     hbN     hbO\n");

	for(j=0;j<clusterN;j++){
		for(k=0;k<cluster[j].N;k++){
			i=result_index[cluster[j].id[k]];
			if(mode==WATMODE) fprintf(clusterf,"%-6d%-5d%8.3f |%8d%8.3f%8.3f%8.3f |  %8.3f:%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",j+1,k+1,cluster[j].Lrmsd[k],dock_result[i].roti,dock_result[i].mx,dock_result[i].my,dock_result[i].mz,dock_result[i].score,dock_result[i].rep,dock_result[i].vdw,dock_result[i].ele,dock_result[i].sol,dock_result[i].wat,dock_result[i].ind,dock_result[i].hbN,dock_result[i].hbO);
			else fprintf(clusterf,"%-6d%-5d%8.3f |%8d%8.3f%8.3f%8.3f |  %8.3f:%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",j+1,k+1,cluster[j].Lrmsd[k],dock_result[i].roti,dock_result[i].mx,dock_result[i].my,dock_result[i].mz,dock_result[i].score,dock_result[i].rep,dock_result[i].vdw,dock_result[i].ele,dock_result[i].sol,dock_result[i].ind,dock_result[i].hbN,dock_result[i].hbO);
		}
	}	
	fclose(clusterf);
}
