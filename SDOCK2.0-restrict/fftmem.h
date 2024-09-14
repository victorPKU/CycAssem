/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * fftmem.h : allocate and free memory.
 ***************************************************************************/

#define fftnew(ptr,nelem) (ptr)=fftalloc((nelem)*sizeof(*(ptr)))
void *fftalloc(unsigned size);
void  fftfree( void *ptr);
