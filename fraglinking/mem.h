#ifndef _mem_h
#define _mem_h


#define snew(ptr,nelem) (ptr)=scalloc(#ptr,(nelem),sizeof(*(ptr)))
#define srenew(ptr,nelem) (ptr)=srealloc(#ptr,(ptr),(nelem)*sizeof(*(ptr)))

void *scalloc(char *name, unsigned nelem,unsigned elsize); 
void *srealloc(char *name, void *ptr,unsigned size);
void  sfree( void *ptr);


#endif

