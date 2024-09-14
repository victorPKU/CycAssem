/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * record.c : save the docking results and sort them by docking score.
 ***************************************************************************/

#include "record.h"
#include <stdio.h>
#include "mem.h"
#include <stdlib.h>
#include <string.h>


void save_dock(int roti, float trans_x,float trans_y,float trans_z,float score,float rep,float vdw,float ele,float sol, float wat, float ind, float hbN, float hbO)
{
	struct LIST *p,*q=NULL;
	int index;
	
	if(count>=RESULTKEEPN && score>=dock_result[order->i].score)
		return;
	if(order==NULL) {
		snew(order,1);
		order->n=NULL;
		p=order;
	}
	else{
		p=order;
		while(p!=NULL){
			if(score<dock_result[p->i].score) {q=p;p=p->n;}
			else if(p!=order){
				snew(q->n,1);
				q->n->n=p;
				p=q->n;
				break;
			}
			else{
				snew(order,1);
				order->n=p;
				p=order;
				break;
			}
		}
		if(p==NULL){
			snew(q->n,1);
			q->n->n=NULL;
			p=q->n;
		}
	}
	if(count<RESULTKEEPN){
		index=count;
		p->i=count;
		count++;
	}
	else{
		index=order->i;
		p->i=order->i;
		p=order->n;
		sfree(order);
		order=p;		
	}
	
	dock_result[index].roti=roti;
	dock_result[index].mx=trans_x;
	dock_result[index].my=trans_y;
	dock_result[index].mz=trans_z;
	dock_result[index].score=score;
	dock_result[index].rep=rep;
	dock_result[index].vdw=vdw;
	dock_result[index].ele=ele;
	dock_result[index].sol=sol;
	dock_result[index].wat=wat;
	dock_result[index].ind=ind;
	dock_result[index].hbN=hbN;
	dock_result[index].hbO=hbO;
}
void arrange_record()
{
	struct LIST *p,*q;
	int i=RESULTKEEPN;

	for(p=order;p!=NULL;){
		i--;
		result_index[i]=p->i;
		q=p->n;
		sfree(p);
		p=q;
	}
}
void init_record()
{
	snew(dock_result,RESULTKEEPN);
	count=0;
	order=NULL;
	snew(result_index,RESULTKEEPN);
}
void free_record()
{
	sfree(dock_result);
	sfree(result_index);
}
#define CLUSTERFILELINE  160
int read_record(char *fn, int *t,int cluster, int number)
{
	FILE *clusterf;
	int i,j;
	char line[CLUSTERFILELINE];
	int rn;
	
	if((clusterf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can't open the sdock result file %s!\n",fn);
		exit(0);
	}
	rn=0;dock_resultf=NULL;
	while(fgets(line, CLUSTERFILELINE,clusterf)){
		if(line[0]=='%'||line[0]==' ') continue;
		i=atoi(line);
		j=atoi(line+6);
		if((cluster==0||(cluster==i))&&((number==0||(number==j)))){
			rn++;
		}
	}
	snew(dock_resultf,rn);

	rewind(clusterf);
	(*t)=0;
	rn=0;
	while(fgets(line, CLUSTERFILELINE,clusterf)){
		if(line[0]=='%') continue;
		if(!strncmp(line+3,"Protein A",9)){
			if(!strncmp(line+14,"mobile",6))
				(*t)=1;
			continue;
		}
		if(!strncmp(line+3,"Protein B",9)) continue;
		i=atoi(line);
		j=atoi(line+6);

		if((cluster==0||(cluster==i))&&((number==0||(number==j)))){
			dock_resultf[rn].cluster=i;
			dock_resultf[rn].num=j;
			dock_resultf[rn].RMSD=atof(line+12);	
			dock_resultf[rn].roti=atoi(line+22);
			dock_resultf[rn].mx=atof(line+29);
			dock_resultf[rn].my=atof(line+37);
			dock_resultf[rn].mz=atof(line+45);
			dock_resultf[rn].score=atof(line+56);
			rn++;
		}
	}

	fclose(clusterf);

	return rn;
}

