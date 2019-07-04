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
typedef void (*fptr)(void**,void**,void**,void**,void**,void**,int*,int*,int*,int*,double*,void**,
			void**,void**,void**,void**,void**,void**,void**,void**,void**,void**,
		        void**, int*, int*,double*, int*, double*, int*, double*,void**, bool*, int*,int*,int*,int*,int*, 
		        double*,double*, void**, void**, double*, double*, double*, double*);
			
/** This is the DelPhi Fortran-C binding. In input there is the needed information to link
with an external module for surface computation. In this version it is assumed that all memory has
already been allocated*/
void nlsolvermodule_(char* module1,
		    void** idpos, void** db, void** sf1, void** sf2, void** iqpos, void** qval,
		    int* icount2a, int* icount2b, int* icount1a, int* icount1b, double* spec, void** phimap, void** phimap1,
		    void** phimap2, void** phimap3, void** ibndx, void** ibndy, void** ibndz, void** bndx1,
		    void** bndx2, void** bndx3, void** bndx4, void** gval, int* ibctyp, int* isolvarch, 
		    double* rionst, int* igrid, double* scale, int* nlit, double* res2, void** inout,
		    bool* imanual, int* nnit, int* ival1a, int* ival1b, int* ival2a, int* ival2b, 
		    double* conc1, double* conc2, void** qmap1, void** qmap2, double* qfact, double* epsout, 
		    double* deblen, double* relpar)
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

    tempPoint = dlsym(handle, "NlDelphiSolverModule");
    if ((error = dlerror()) != NULL)  
    {
	fputs(error, stderr);
	printf("\n\n<<ERROR>> Cannot load function NlDelphiSolverModule in plug-in named %s\n",module);
	exit(1);
    }
    printf("ok!\n");
	    
    fptr NlDelphiSolverModule = (fptr)tempPoint;
				
    // propagates the called function
    (*NlDelphiSolverModule)(idpos,db,sf1,sf2,iqpos,qval,icount2a,icount2b,icount1a,icount1b,spec,phimap,
	    phimap1,phimap2,phimap3,ibndx,ibndy,ibndz, bndx1,bndx2,bndx3,bndx4,gval,ibctyp,isolvarch,rionst,
	    igrid,scale,nlit,res2,inout,imanual,nnit,ival1a,ival1b,ival2a,ival2b,conc1,conc2,qmap1,qmap2,qfact,epsout,deblen,relpar);		
		
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
	fptr NlDelphiSolverModule = (fptr)GetProcAddress(delphi_solv,"NlDelphiSolverModule");
	if(NlDelphiSolverModule == NULL) 
	{			
	    printf("\n\n<<ERROR>> Cannot load function NlDelphiSolverModule in plug-in named %s\n",module);
	    exit(1);
	} 
	else 
	{								
	    (*NlDelphiSolverModule)(idpos,db,sf1,sf2,iqpos,qval,icount2a,icount2b,icount1a,icount1b,spec,phimap,
	    phimap1,phimap2,phimap3,ibndx,ibndy,ibndz, bndx1,bndx2,bndx3,bndx4,gval,ibctyp,isolvarch,rionst,
	    igrid,scale,nlit,res2,inout,imanual,nnit,ival1a,ival1b,ival2a,ival2b,conc1,conc2,qmap1,qmap2,qfact,epsout,deblen,relpar);		
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
