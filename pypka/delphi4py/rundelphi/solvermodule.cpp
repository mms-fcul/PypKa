#include <stdlib.h>
#include <stdio.h>
		
#if  (defined _WIN32) || (defined __WIN32__) || (defined __TOS_WIN__) || (defined __WINDOWS__) 
		#define WIN_OS
#elif (defined __linux) || (defined __unix) || (defined macintosh) || (defined Macintosh) || (defined __APPLE__ && defined __MACH__) || (defined __ANDROID__)
		#define LINUX_OS
#endif

#ifdef LINUX_OS
    #include <dlfcn.h>
#elif defined WIN_OS
    #include <windows.h>		
#endif

#  ifdef __cplusplus
	extern "C" {
#  endif
	
/* function pointer type*/
typedef void (*fptr)(void**, void**, void**, void**, void**, void**, int*,
	int*, int*, int*, double*, void**, void**, void**, void**, void**, void**, void**,
	void**, void**, void**, void**, void**, int*, int*, double*, int*, double*,
	int*, double*, void**,double*,void**,double*,void**);
			
/** This is the DelPhi Fortran-C binding. In input there is the needed information to link
with an external module for surface computation. In this version it is assumed that all memory has
already been allocated*/

void solvermodule_(char* module1,
	void** idpos_,void**  db_,void** sf1_,void** sf2_, void** iqpos_, void** qval_,int* icount2a_,
	int* icount2b_,int* icount1a_,int* icount1b_, double* spec_, void** phimap_, void** phimap1_,void** phimap2_,
	void** phimap3_, void** ibndx_,void** ibndy_,void** ibndz_,void** bndx1_, void** bndx2_, void**  bndx3_,
	void** bndx4_,void**  gval_, int* ibctyp, int* isolvarch, double* rionst, int* igrid_, double* scale_,
	int* nlit, double* res2, void** iper, double* qplus_, void** cqplus_, double* qmin_, void** cqmin_)
{
    char module[1024];
    
    printf("\n<<INFO>> Loading module %s....",module1);
    
#ifdef LINUX_OS	

    void* tempPoint;
    char* error;
    void *handle;	 	  

    sprintf(module,"%s.so",module1);	
    handle = dlopen (module, RTLD_LAZY);
    if (!handle) 
    {	    
		fputs (dlerror(), stderr);
		printf("\n\n<<ERROR>> Cannot load Solver module named %s\n\n",module);
		exit(1);
    }

    tempPoint = dlsym(handle, "delphiSolverModule");
    if ((error = dlerror()) != NULL)  
    {
		fputs(error, stderr);
		printf("\n\n<<ERROR>> Cannot load function delphiSolverModule in plug-in named %s\n",module);
		exit(1);
    }
    
    printf("ok!\n");
	    
    fptr libPbsolver4Delphi = (fptr)tempPoint;
				
    // propagates the called function
    (*libPbsolver4Delphi)(idpos_,db_,sf1_,sf2_,iqpos_,qval_,icount2a_,icount2b_,icount1a_,icount1b_,spec_,phimap_,
	phimap1_,phimap2_,phimap3_,ibndx_,ibndy_,ibndz_,bndx1_,bndx2_,bndx3_,bndx4_,gval_,ibctyp,isolvarch,
	rionst,igrid_,scale_,nlit,res2,iper,qplus_,cqplus_,qmin_,cqmin_);
	
    dlclose(handle);

#elif defined WIN_OS
    
    sprintf(module,"%s.dll",module1);	
    HMODULE delphi_solv = LoadLibraryEx(TEXT(module),NULL,LOAD_WITH_ALTERED_SEARCH_PATH);
    if( delphi_solv == NULL) 
    {
		printf("\n\n<<ERROR>> Cannot load Solver module named %s\n\n",module);
		exit(1);
    } 
    else 
    {
		fptr libPbsolver4Delphi = (fptr)GetProcAddress(delphi_solv,"libPbsolver4Delphi");
		if(libPbsolver4Delphi == NULL) 
		{			
		    printf("\n\n<<ERROR>> Cannot load function libPbsolver4Delphi in plug-in named %s\n",module);
		    exit(1);
		} 
		else 
		{								
		    (*delphiSolverModule)(idpos_,db_,sf1_,sf2_,iqpos_,qval_,icount2a_,icount2b_,icount1a_,icount1b_,spec_,phimap_,
		    phimap1_,phimap2_,phimap3_,ibndx_,ibndy_,ibndz_,bndx1_,bndx2_,bndx3_,bndx4_,gval_,ibctyp,isolvarch,
		    rionst,igrid_,scale_,nlit,res2,iper,qplus_,cqplus_,qmin_,cqmin_);
		}
    }
    
#else

    printf("\n<<ERROR>> Cannot load %s, no OS loader defined",module1);
    exit(1);
    
#endif
    
    return;	  
}

#  ifdef __cplusplus
}
#  endif
