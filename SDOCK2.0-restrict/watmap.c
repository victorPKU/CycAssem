/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * watmap.c : generate waters around the protein surface.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "geometry.h"
#include "mem.h"
#include "structure.h"
#include "watgrid.h"

void calcorefluc(STRUCTURE *pro)
{
	float averagefluc;
	int ia,an;

	averagefluc=0.0;an=0;
	for(ia=0;ia<pro->atmN;ia++){
		if(pro->a[ia].surf_core==COREATOM){
			if(pro->a[ia].fluc>0.01){
				averagefluc+=pro->a[ia].fluc;
				an++;
			}
		}
	}
	averagefluc/=an;

	pro->fluc0=averagefluc;
}
void print_watmap_usage()
{
	printf(" Usage: watmap preprocessed_protein_structure watermap_for_protein\n\n");
	printf(" Examples:\n");
	printf(" watmap 1AY7R.pdb 1AY7Rwat.pdb \n");
}

int main(int argc, char *argv[])
{
	STRUCTURE pro;
	int *mattergrid;
	int *pocketgrid;
	float *HBgrid;
	float *watenegrid;
	int enewatn;
	float (*enewatxyz)[4];

	printf("                             SDOCK2.0\n");
	printf("  Global Protein-Protein Docking Using Stepwise Force-Field Potentials   \n");
	printf(" Copyright (c) 2024 Molecular design lab in Peking university, Beijing, China\n\n");
	printf(" watmap: build the water map on the protein surface\n");

	if(argc==1){
		print_watmap_usage();
		return EXIT_SUCCESS;
	}

	if(argv[1][0]=='-'){
		printf("\nERROR: The first option should be a preprocessed protein structure file, and not %s!\n\n",argv[1]);
		print_watmap_usage();
		return EXIT_SUCCESS;
	}
	if(argc<=2 || argv[2][0]=='-'){
		printf("\nERROR: You only input the name for the output water map file!\n\n");
		print_watmap_usage();
		return EXIT_SUCCESS;
	}

	read_structure(&pro,argv[1]);
	fixpro(&pro);
	calcorefluc(&pro);

	decide_gridsize(&pro);
	snew(mattergrid, GRID3SIZE);
	init_mattergrid(mattergrid, GRID3SIZE);
	gen_mattergrid(&pro, mattergrid);
	printf("generate mattergrid\n");

	snew(HBgrid, GRID3SIZE);
	init_HBgrid(HBgrid, GRID3SIZE);
	gen_HBgrid(&pro, HBgrid);
	printf("generate HBgrid\n");

	snew(pocketgrid, GRID3SIZE);
	init_pocketgrid(mattergrid,pocketgrid, GRID3SIZE);
	gen_pocketgrid(&pro, pocketgrid,mattergrid);
	printf("generate pocketgrid\n");

	snew(watenegrid, GRID3SIZE);
	init_watenegrid(watenegrid, GRID3SIZE);
	gen_watenegrid(watenegrid,mattergrid, pocketgrid,HBgrid);
	printf("generate watenegrid\n");

	snew(enewatxyz, pro.surfatN*10);
	enewatn=enegrid2watpositin(mattergrid, pocketgrid,watenegrid,HBgrid,enewatxyz,0,1);
	enewatn+=enegrid2watpositin(mattergrid, pocketgrid,watenegrid,HBgrid,enewatxyz,enewatn,-1);
	sfree(watenegrid);
	printwatxyz(argv[1],argv[2],enewatn,enewatxyz);
	sfree(enewatxyz);
	printf("print pocket water %s!\n", argv[2]);

	sfree(mattergrid);
	sfree(HBgrid);
	sfree(pocketgrid);

	sfree(pro.a);

  return EXIT_SUCCESS;
}
