/* call to the C function realloc is necessary as memory is adjusted at
   run time 			S. Sridharan. Sep 95 */
#include <stdlib.h>
#include <stdio.h> 

#ifdef IRIX
#define memalloc memalloc_
#elif defined (LINUX)
#define memalloc memalloc_
#elif defined (CRAY)
#define memalloc MEMALLOC
#elif defined (PC)
#define memalloc (_stdcall MEMALLOC)
#endif 
  
void *memalloc( void **ptr,int *entry_size, int *new_size )
{       void *newptr;   
/*      printf("\n Entry_size: %d, new_size: %d\n",*entry_size, *new_size);*/
	if(!*new_size) 
    {
       if(*ptr) free(*ptr);
       return NULL;
	}
	if(!*ptr) { 
                  newptr=calloc((size_t)(*new_size),(size_t)(*entry_size));  
                 /*printf("\n *ptr: %p, newptr: %p, newptr int : %ld\n",*ptr,newptr,newptr);*/
                  }
	else      { newptr=realloc(*ptr, (size_t)((*new_size)*(*entry_size)));
                    /*printf("\n realloc: newptr: %p, newptr int : %d\n",*ptr,newptr);*/ }

	if (newptr == 0 && *new_size != 0) 
        {
          perror("memalloc");
	  exit(EXIT_FAILURE);
	}

	return newptr;
}
