/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * mem.h : allocate and free memory.
 ***************************************************************************/
#ifndef _mem_h
#define _mem_h


#define snew(ptr,nelem) (ptr)=scalloc(#ptr,(nelem),sizeof(*(ptr)))
#define srenew(ptr,nelem) (ptr)=srealloc(#ptr,(ptr),(nelem)*sizeof(*(ptr)))

void *scalloc(char *name, unsigned nelem,unsigned elsize); 
void *srealloc(char *name, void *ptr,unsigned size);
void  sfree( void *ptr);


#endif

