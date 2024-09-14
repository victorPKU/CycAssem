/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * mem.c : allocate and free memory.
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

void *scalloc(char *name, unsigned nelem,unsigned elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL){
		  printf("ERROR:  Memery Alloc Error for %s\n", name);
		 exit(0); 
	  }
     }
  return p;
}
void *srealloc(char *name, void *ptr,unsigned size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if (ptr==NULL) 
	p=malloc((size_t)size); 
      else 
	p=realloc(ptr,(size_t)size);
      if (p==NULL){ 
		 printf("ERROR:  Memery Alloc Error for %s\n", name);
		 exit(0); 
	  }
    }
  return p;
}

void sfree(void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}
