/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * build.c : build models according to the docking records.
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "build.h"
#include "rotplan.h"
#include "record.h"
#include "mem.h"

void get_model_name(char *mfnall, int i,int j,char *code,char *route)
{
	char cluster[5], num[4];

	cluster[0]='0'+(i)%10000/1000;
	cluster[1]='0'+(i)%1000/100;
	cluster[2]='0'+(i)%100/10;
	cluster[3]='0'+(i)%10;
	cluster[4]='\0';

	num[0]='0'+(j)%1000/100;
	num[1]='0'+(j)%100/10;
	num[2]='0'+(j)%10;
	num[3]='\0';

	sprintf(mfnall,"%s/SDOCK_%s_%s_%s.pdb",route,code,cluster,num);
}


void print_build_usage()
{
	printf(" Usage: preprocessed_protein_structure1 preprocessed_protein_structure2 \n");
	printf("    [-o  sdock_result_file_name (sdock_record)]\n");
	printf("    [-r rotation_sampling_file (so3layer.qua)] \n");
	printf("    [-s shift_vector (0,0,0)] \n");
	printf("    [-c cluster_number (1)] \n");
	printf("    [-m cluster_member_number (1)]\n");
	printf("    [-n model_name (AAAA)]\n");
	printf("    [-l ligand_only (0)]\n");
	printf("    [-d model_file_route (./)]\n");
	printf(" \n");
	printf(" Examples:\n");
	printf(" build 1AY7R.pdb 1AY7L.pdb\n");
	printf(" build 1AY7R.pdb 1AY7L.pdb -n 1ay7\n      #The first 4 letter of file name of the built complex model is 1ay7\n");
	printf(" build 1AY7R.pdb 1AY7L.pdb -r ../data/so3layer.qua\n");
	printf(" build 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_resultf -c 0 -m 1 -n 1AY7\n      #The first structure of all clusters\n");
	printf(" build 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_resultf -c 0 -m 0 -n 1AY7 -d model/\n      #build the structure of all results to the direction model/\n");
	printf(" build target.pdb frag.pdb -o frag_dock -c 1 -m 1 -n frag -l 1 -s 11.798,0.005,-0.002\n      #Only output liand and Shift the structure with vector [11.798,0.005,-0.002]\n");
}
int get_shift(float shift[3], char *arg)
{
	int i=0,x=0;
	int l;
	l=strlen(arg);

	shift[0]=atof(arg);
	while(i<l){
		if(arg[i]==','){
			x++;
			if(x<3){
				shift[x]=atof(arg+i+1);
			}
		}
		i++;
	}
	return x;
}

int main(int argc, char *argv[])
{
	int writeN;
	int t=0;
	int mi;
	char mfnall[150];
	float rotm[3][3];
	/*char mfn[100]="SDOCK_AAAA_****_***.pdb";*/

	int cluster;
	int number;	
	int arg_result;
	int arg_rotfile;
	int arg_name;
	int arg_route;

	STRUCTURE receptor,ligand;

	int arg_i;
	int canshift=0;
	float shift[3];
	int liandOnly=0;

	printf("                             SDOCK2.0\n");
	printf("  Global Protein-Protein Docking Using Stepwise Force-Field Potentials   \n");
	printf(" Copyright (c) 2024 Molecular design lab in Peking university, Beijing, China\n\n");
	printf(" build : build models according to the docking records.\n");

	if(argc==1){
		print_build_usage();
		return EXIT_SUCCESS;
	}

	if(argv[1][0]=='-'){
		printf("\nERROR: The first two parameters should be two structure files, and not %s!\n\n",argv[1]);
		print_build_usage();
		return EXIT_SUCCESS;
	}
	if(argc==2){
		printf("\nERROR: You only input one structure file and you should input two!\n\n");
		print_build_usage();
		return EXIT_SUCCESS;
	}

	if(argv[2][0]=='-'){
		printf("\nERROR: The first two parameters should be two structure files, and not %s!\n\n",argv[2]);
		print_build_usage();
		return EXIT_SUCCESS;
	}

	for(arg_i=3;arg_i<argc;arg_i+=2){
		if(argv[arg_i][0]!='-'||(argv[arg_i][1]!='o'&&argv[arg_i][1]!='r'&&argv[arg_i][1]!='c'&&argv[arg_i][1]!='s'&&argv[arg_i][1]!='l'&&argv[arg_i][1]!='m'&&argv[arg_i][1]!='n'&&argv[arg_i][1]!='d')){
			printf("\nERROR at '%s': should be one of these options: -o -r -c -s -l -m -n -d!\n\n",argv[arg_i]);
			print_build_usage();
			return EXIT_SUCCESS;
		}
		if(arg_i+1>=argc||(argv[arg_i][1]!='s'&&argv[arg_i+1][0]=='-')){
			if(arg_i+1<argc)
				printf("\nERROR at '%s %s': ",argv[arg_i],argv[arg_i+1]);
			else
				printf("\nERROR at '%s': ",argv[arg_i]);
			if(argv[arg_i][1]=='o')
				printf("please input the docking result file name after -o!\n\n");
			else if(argv[arg_i][1]=='r')
				printf("please input the rotation sampling file after -r!\n\n");
			else if(argv[arg_i][1]=='c')
				printf("please input a cluster number after -c!\n\n");
			else if(argv[arg_i][1]=='m')
				printf("please input a cluster member number after -m!\n\n");
			else if(argv[arg_i][1]=='n')
				printf("please input model name after -n!\n\n");
			else if(argv[arg_i][1]=='d')
				printf("please input the directory for the to-be-built models after -d!\n\n");
			else
				printf("please input the option value after %s!\n\n",argv[arg_i]);
			print_build_usage();
			return EXIT_SUCCESS;
		}
	}

	cluster=1;
	number=1;	
	arg_result=0;
	arg_rotfile=0;
	arg_name=0;
	arg_route=0;
	shift[0]=shift[1]=shift[2]=0.0;

	for(arg_i=3;arg_i<argc;arg_i+=2){
		if(argv[arg_i][1]=='o')
			arg_result=arg_i+1;
		else if(argv[arg_i][1]=='r'){
			arg_rotfile=arg_i+1;
		}
		else if(argv[arg_i][1]=='c'){
			cluster=atoi(argv[arg_i+1]);
		}
		else if(argv[arg_i][1]=='m'){
			number=atoi(argv[arg_i+1]);
		}
		else if(argv[arg_i][1]=='l'){
			liandOnly=atoi(argv[arg_i+1]);
		}
		else if(argv[arg_i][1]=='n'){
			arg_name=arg_i+1;
		}
		else if(argv[arg_i][1]=='d'){
			arg_route=arg_i+1;
		}
		else if(argv[arg_i][1]=='s'){
			canshift=get_shift(shift,argv[arg_i+1]);
			if(canshift<2){
				shift[0]=shift[1]=shift[2]=0.0;
				printf("\nWARNING wrong shift vector '%s'. No shift instead\n",argv[arg_i+1]);
			}
		}
	}

	if(arg_rotfile==0)
		ROTN=readrotplan("so3layer.qua");
	else
		ROTN=readrotplan(argv[arg_rotfile]);
	read_structure(&receptor,argv[1]);
	read_structure(&ligand,argv[2]);
	/*centerpro(&receptor);
	centerpro(&ligand);*/

	
	if(arg_result==0){
		writeN=read_record("sdock_record", &t,cluster,number);
	}
	else{
		writeN=read_record(argv[arg_result], &t,cluster,number);
	}
	if(writeN==0){
		if(arg_result==0)
			printf("WARNING: The cluster %d member %d is not found in the docking result file sdock_record!\n",cluster,number);
		else
			printf("WARNING: The cluster %d member %d is not found in the docking result file %s!\n",cluster,number,argv[arg_result]);
	}
	
	for(mi=0; mi<writeN; mi++){
		quat2matrix(ROTPLAN[dock_resultf[mi].roti], rotm);
		build_model(&receptor, &ligand, t,rotm, dock_resultf[mi].mx, dock_resultf[mi].my, dock_resultf[mi].mz,shift);
		if(arg_name==0&&arg_route==0)
			get_model_name(mfnall, dock_resultf[mi].cluster,dock_resultf[mi].num,"AAAA","./");
		else if(arg_name!=0&&arg_route==0)
			get_model_name(mfnall, dock_resultf[mi].cluster,dock_resultf[mi].num,argv[arg_name],"./");
		else if(arg_name==0&&arg_route!=0)
			get_model_name(mfnall, dock_resultf[mi].cluster,dock_resultf[mi].num,"AAAA",argv[arg_route]);
		else
			get_model_name(mfnall, dock_resultf[mi].cluster,dock_resultf[mi].num,argv[arg_name],argv[arg_route]);
		if(arg_result==0)
			write_complex(&receptor, &ligand,mfnall,dock_resultf[mi].cluster,dock_resultf[mi].num,dock_resultf[mi].score,argv[1],argv[2],t,"sdock_record",shift,liandOnly);
		else
			write_complex(&receptor, &ligand,mfnall,dock_resultf[mi].cluster,dock_resultf[mi].num,dock_resultf[mi].score,argv[1],argv[2],t,argv[arg_result],shift,liandOnly);
	}
	sfree(dock_resultf);
	sfree(receptor.a);
	sfree(ligand.a);
	return EXIT_SUCCESS;
}
