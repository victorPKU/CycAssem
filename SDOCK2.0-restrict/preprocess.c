/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * preprocess.c : preprocess the origin protein structure file (PDB format) to the SDOCK input file.
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "protein.h"
#include "GNM.h"
#include "mem.h"


void print_preprocess_usage()
{
	printf(" Usage: preprocess pdb_file \n");
	printf("    [-c chain_id (all chains in the pdb_file)]  \n");
	printf("    [-o output_file (preprocessed.pdb)]  \n");
	printf("    [-d dynamic_or_not (1)]  \n");
	printf("    [-u uncommon_amino_acid (no)] \n");
	printf("    [-s small_molecule (no)]  \n");
	printf("    [-a atom_parameter_file (ATM)]\n");
	printf("    [-m center_coordination_of_docking_box (0)\n\n");
	printf(" Format of uncommon_amino_acid or small_molecule: \n");
	printf("    chain_id1+name1:chain_id2+name2:chain_id3+name3:...  '@' means any chain, '_' means blank\n");
	printf(" \n");
	printf(" Examples:\n");
	printf(" preprocess 1SBB.pdb\n");
	printf(" preprocess 1SBB.pdb -u no -d 0    #do not consider dynamic using GNM model\n");
	printf(" preprocess 1FKM.pdb -c A -o 2G77L.pdb -u AMSE\n");
	printf(" preprocess 1EXB.pdb -o 1EXBR.pdb -s  @NDP -a myATM\n");
	printf(" preprocess 1ZMC.pdb -c AB  -o 2F5ZR.pdb -s  AFAD:BFAD -a ../data/ATM  # A and B chains\n");
	printf(" preprocess target.pdb -o targetbox.pdb -m -4.173,-29.160,-4.232\n     # put [-4.173,-29.160,-4.232] as the box center\n");
}
int get_het(char het[][5],char *arg)
{
	int i,j=0;
	int c=0;
	int strn=strlen(arg);
	
	if(!strcmp(arg,"no")) return 0;
	for(i=0;i<strn;i++){
		if(arg[i]!=':'){
			if(arg[i]=='_')
				het[c][j]=' ';
			else
				het[c][j]=arg[i];
			j++;
		}
		else{
			het[c][j]='\0';
			c++;
			j=0;
		}
		
	}
	het[c][j]='\0';
	return c+1;
}
int get_fixcenter(float fixcenter[3], char *arg)
{
	int i=0,x=0;
	int l;
	l=strlen(arg);

	fixcenter[0]=atof(arg);
	while(i<l){
		if(arg[i]==','){
			x++;
			if(x<3){
				fixcenter[x]=atof(arg+i+1);
			}
		}
		i++;
	}
	return x;
}

int main(int argc, char *argv[])
{
	PROTEIN pro;

	int arg_chain;
	int arg_outfile;
	int arg_atmfile;
	int arg_ucAA;
	int arg_smol;
	
	int ucAAN;
	char ucAA[20][5];
	int smolN;
	char smol[20][5];

	double *matrix=NULL;
	int arg_i;
	int dynamic=0;

	float fixcenter[3];
	int usefixcenter=0;	
	

	printf("                             SDOCK2.0\n");
	printf("  Global Protein-Protein Docking Using Stepwise Force-Field Potentials   \n");
	printf(" Copyright (c) 2024 Molecular design lab in Peking university, Beijing, China\n\n");
	printf(" preprocess: preprocess the origin protein structure file (PDB format) to the SDOCK input file\n");
	
	if(argc==1){
		print_preprocess_usage();
		return EXIT_SUCCESS;
	}

	if(argv[1][0]=='-'){
		printf("\nERROR: The first option should be a pdb file, and not %s!\n\n",argv[1]);
		print_preprocess_usage();
		return EXIT_SUCCESS;
	}

	for(arg_i=2;arg_i<argc;arg_i+=2){
		if(argv[arg_i][0]!='-'||(argv[arg_i][1]!='c'&&argv[arg_i][1]!='o'&&argv[arg_i][1]!='d'&&argv[arg_i][1]!='u'&&argv[arg_i][1]!='s'&&argv[arg_i][1]!='m'&&argv[arg_i][1]!='a')){
			printf("\nERROR at '%s': should be one of these options: -c -o -d -u -s -a -m!\n\n",argv[arg_i]);
			print_preprocess_usage();
			return EXIT_SUCCESS;
		}
		if(arg_i+1>=argc||(argv[arg_i][1]!='m'&&argv[arg_i+1][0]=='-')){
			printf("\nplease input the option value after %s!\n\n",argv[arg_i]);
			print_preprocess_usage();
			return EXIT_SUCCESS;
		}
	}	

	arg_chain=0;
	arg_outfile=0;
	arg_atmfile=0;
	arg_ucAA=0;
	arg_smol=0;
	ucAAN=0;
	smolN=0;
	dynamic=1;

	for(arg_i=2;arg_i<argc;arg_i+=2){
		if(argv[arg_i][1]=='c')
			arg_chain=arg_i+1;
		else if(argv[arg_i][1]=='o')
			arg_outfile=arg_i+1;
		else if(argv[arg_i][1]=='d')
			dynamic=atoi(argv[arg_i+1]);
		else if(argv[arg_i][1]=='u'){
			arg_ucAA=arg_i+1;
			ucAAN=get_het(ucAA,argv[arg_i+1]);
		}
		else if(argv[arg_i][1]=='s'){
			arg_smol=arg_i+1;
			smolN=get_het(smol,argv[arg_i+1]);
		}
		else if(argv[arg_i][1]=='a')
			arg_atmfile=arg_i+1;
		else if(argv[arg_i][1]=='m'){
			usefixcenter=get_fixcenter(fixcenter,argv[arg_i+1]);
		}
	}

	if(arg_chain==0)
		readprotein(&pro, argv[1], "@", ucAAN, ucAA, smolN, smol);
	else
		readprotein(&pro, argv[1], argv[arg_chain], ucAAN, ucAA, smolN, smol);
	if(arg_atmfile==0)
		get_atom_para(&pro, "ATM");
	else
		get_atom_para(&pro, argv[arg_atmfile]);

	if(usefixcenter==0)
		centerpro(&pro);
	else
		fixcenterpro(&pro,fixcenter);

	cal_surface(&pro,WATERRAD);
	if(dynamic) {
		genblocks(&pro);
	      snew(matrix,pro.blockN*pro.blockN);
		genmatrix(&pro, matrix);
		calfluctuation(matrix, &pro);
		sfree(matrix);
	}
	if(arg_outfile==0)
		write_structure(&pro, "preprocessed.pdb",argv,arg_chain,arg_atmfile,arg_ucAA,arg_smol);
	else
		write_structure(&pro, argv[arg_outfile],argv,arg_chain,arg_atmfile,arg_ucAA,arg_smol);

	free_pro(&pro);

	return EXIT_SUCCESS;
}
