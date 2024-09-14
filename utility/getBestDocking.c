#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BESTNUM 10000

int add2best(float score,float *bestscore,int *index, int c, int N)
{
	int i,j;
	int w;

	if(c<N){
		i=c;w=c;
	}
	else{
		if(score>bestscore[index[N-1]]) return -1;
		i=N;w=index[N-1];
	}
	while(i>=1&&score<bestscore[index[i-1]]){
		index[i]=index[i-1];
		i--;
	}
	index[i]=w;
	return w;
}
int main(int argc, char *argv[])
{
	FILE *f1,*f2;
	char line0[100];
	char dockf[100];
	char line[100];
	int i,j;
	float score;
	int id;
	float bestscore[BESTNUM];
	int index[BESTNUM];
	char bestpep[BESTNUM][100];
	char bestdock[BESTNUM][100];
	char pep[100];
	int bestno[BESTNUM];
	int bestN;
	int c=0;
	int l,m;

	if(argc>3){
		bestN=atoi(argv[3]);
		if(bestN<1||bestN>BESTNUM) bestN=BESTNUM;
	}
	else bestN=BESTNUM;

	for(i=0;i<bestN;i++){
		bestscore[i]=1.0;
		index[i]=i;
		bestpep[i][0]='\0';
		bestdock[i][0]='\0';
	}
	if((f1=fopen(argv[1],"r"))==NULL){
		printf("ERROR: Can not open file %s.\n", argv[1]);
		exit(0);
	}
	while(fgets(line0,100,f1)){
		l=strlen(line0);
		if(l<=1||l>=100) continue; 
		strncpy(dockf,line0,l-1);
		dockf[l-1]='\0';

		if((f2=fopen(dockf,"r"))==NULL){
			printf("ERROR: Can not open file %s.\n", dockf);
			exit(0);
		}
		while(fgets(line,100,f2)){
			if(line[0]=='%')continue;
			if(line[0]==' '){
				if(!strncmp("mobile",line+14,6)){
					m=strlen(line);
					strncpy(pep,line+21,m-22);
					pep[m-22]='\0';
				}
				continue;
			}
			i=atoi(line);
			j=atoi(line+6);
			if(j>1) continue;
			score=atof(line+56);
			id=add2best(score,bestscore,index,c,bestN);
			if(id!=-1){
				bestscore[id]=score;
				strncpy(bestdock[id],dockf,l-1);
				bestdock[id][l-1]='\0';
				strcpy(bestpep[id],pep);
				bestno[id]=i;
			} 
			c++;
		}
		fclose(f2);

	}
	fclose(f1);
	if((f1=fopen(argv[2],"w"))==NULL){
		printf("ERROR: Can not open file %s.\n", argv[2]);
		exit(0);
	}

	for(i=0;i<bestN;i++){
		fprintf(f1,"%s %s %2d %8.3f\n",bestpep[index[i]],bestdock[index[i]], bestno[index[i]], bestscore[index[i]] );
	}

  return EXIT_SUCCESS;

}
