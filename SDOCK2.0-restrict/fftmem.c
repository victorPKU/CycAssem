/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * fftmem.c : allocate and free memory.
 ***************************************************************************/

#include "fftmem.h"
#include <fftw3.h>

void *fftalloc(unsigned size)
{
	void *p;
	p=fftwf_malloc(size);

	return p;
}
void  fftfree( void *ptr)
{
	fftwf_free(ptr);
}
