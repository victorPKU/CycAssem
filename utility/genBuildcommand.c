#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINELEN 300

int main(int argc, char *argv[])
{
	FILE *in,*out;
	char line[LINELEN];
	char modelname[100];
	int No;
	int i=0,j,k;
	int n=1;

	if((out=fopen(argv[2],"w"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",argv[2]);
		exit(0);
	}
	if((in=fopen(argv[1],"r"))==NULL){
		printf("ERROR: Can not open protein structure file %s\n",argv[1]);
		exit(0);
	}

	while(fgets(line, LINELEN,in)){
		fprintf(out,"%s %s  ", argv[3],argv[4]); /*../SDOCK2.0-restrict/build ALK1box1.pdb*/
		i=0; while(line[i]!=' '){
			fprintf(out,"%c",line[i]);
			i++; 
		}
		fprintf(out,"  -o ");
		while(line[i]==' ') i++;
		k=0;
		while(line[i]!=' '){
			fprintf(out,"%c",line[i]);
			modelname[k]=line[i];
			k++;
			if(line[i]=='/') k=0;
			i++; 
		}
		modelname[k]='\0';
		fprintf(out,"  -n ");
		fprintf(out,"%s",modelname);
		fprintf(out,"  -d ");
		fprintf(out,"%s",argv[5]);/*dockmodel: docking model directory*/

		No=atoi(line+i);
		fprintf(out,"  -c ");
		fprintf(out,"%2d",No);

		fprintf(out,"  -m 1  -l 1  -r %s  -s %s\n",argv[6],argv[7]);
	}
	
	fclose(in);
	fclose(out);
}

