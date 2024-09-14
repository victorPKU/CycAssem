#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXMODEL 25
int selectmod(char *fn, float dG, int hbn, float tot, int bestN, int *index)
{
	FILE *pf;
	float score[MAXMODEL];
	int HB[MAXMODEL];
	int modn=0,imod;
	char line[100];
	int i,j,k;
	int maxi;

	if((pf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",fn);
		exit(0);
	}
	while(fgets(line, 100,pf)){
		if(!strncmp(line,"TITLE",5)){
			score[modn]=0;
			HB[modn]=0;
			modn++;
		}
		else if(!strncmp(line,"REMARK",6)){
			if(!strncmp(line+8,"score",5)){
				score[modn-1]=atof(line+13);
			}
			else if(!strncmp(line+8,"backbone HB No",14)){
				HB[modn-1]=atoi(line+23);
			}
		}
	}
	fclose(pf);

	i=0;
	for(imod=0;imod<modn;imod++){
		if(score[imod]<=dG&&HB[imod]>=hbn&&(score[imod]-HB[imod]*1.5)<=tot){
			if(i<bestN){
				index[i]=imod;
				i++;
			}
			else{
				index[bestN]=imod; maxi=bestN;
				for(j=0;j<bestN;j++){
					if((score[maxi]-HB[maxi]*1.5)<(score[j]-HB[j]*1.5)) maxi=j;
				}
				if(maxi!=bestN) index[maxi]=imod;
			}
		}
	}
	return i;
}
void writetarget(FILE *f1, char *target)
{
	FILE *pf;
	char line[100];
	
	if((pf=fopen(target,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",target);
		exit(0);
	}
	while(fgets(line, 100,pf)){
		if(!strncmp(line,"ATOM",4)){
			line[56]='\n';
			fprintf(f1,"%s",line);
		}
	}
	fclose(pf);
}
int selected(int modn, int N, int *index)
{
	int i;
	for(i=0;i<N;i++) if(modn==index[i]) return 1;
	return 0;
}

void getname(char *cmplx, char *fn, char *dir, int modn)
{
	int i,j=0;
	char cyc[100];

	i=0;while(fn[i]!='\0'){
		cyc[j]=fn[i];
		if(fn[i]=='.') { cyc[j]='\0';break;} 
		j++;
		if(fn[i]=='/')j=0;
		i++;
	}
	sprintf(cmplx, "%s/%s_%d.pdb", dir, cyc,modn+1);
}
void buildcomplex(char *fn,int N,int *index,char *dir,char *target)
{
	FILE *pf,*cf;

	int modn=0;
	char line[100];
	char cmplx[100];
	int start=0;
	int atmstart=0;

	if((pf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",fn);
		exit(0);
	}
	while(fgets(line, 100,pf)){
		if(!strncmp(line,"TITLE",5)){
			if(selected(modn,N,index)==1){
				getname(cmplx,fn,dir,modn);
				if((cf=fopen(cmplx,"w"))==NULL){
					printf("ERROR: Can not open protein structure file %s\n",cf);
					exit(0);
				}
				start=1;
			}
			modn++;
		}
		else if(start){
			if(!strncmp(line,"ATOM",4)){
				if(atmstart==0){
					writetarget(cf,target);
					fprintf(cf,"TER\n");
					atmstart=1;
				}	
			}
			fprintf(cf,"%s",line);
			if(!strncmp(line,"END",3)){
				fclose(cf);atmstart=0;start=0;
			}
		}
	}
	fclose(pf);
	
}

int main(int argc, char *argv[])
{
	int index[MAXMODEL];
	float dG=atof(argv[2]);
	int hbn=atoi(argv[3]);
	float tot=atof(argv[4]);
	int bestN=atoi(argv[5]);
	int N;

	N=selectmod(argv[1], dG, hbn, tot, bestN, index);
	buildcomplex(argv[1],N,index,argv[6],argv[7]);
	printf("====> %2d complex model built for %s\n",N,argv[1]);
}
