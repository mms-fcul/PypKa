	#include <stdlib.h>
	#include <stdio.h>
		
	#if  (defined _WIN32) || (defined __WIN32__) || (defined __TOS_WIN__) || (defined __WINDOWS__) 
			#define WIN_OS
	#elif (defined __linux) || (defined __unix) ||  (defined __ANDROID__)
			#define LINUX_OS
	#elif (defined macintosh) || (defined Macintosh) || (defined __APPLE__ && defined __MACH__)
			#define MAC_OS
	#endif

	#if (defined LINUX_OS) || (defined MAC_OS)
        #include <dlfcn.h>
	#elif defined WIN_OS
		#include <windows.h>
	#endif

	#if (defined MSDIG) || (defined _WIN32 && defined __INTEL_COMPILER)
		#define surfacemodule SURFACEMODULE
	#else
		#define surfacemodule surfacemodule_
	#endif
 

	#  ifdef __cplusplus
		extern "C" {
	#  endif
		
	/* function pointer type*/
    typedef void (*fptr)(double,double,double,
                         double,double,double,
                         double,double,double,double,double,
                         int*,int,double,double*,
                         double*,bool*,int*,int,int*,int,
                         double*,double*,int*,int,double,double,char*,int*);
	  
					  
	/** This is the DelPhi Fortran-C binding. In input there is the needed information to link
	with an external module for surface computation*/
	void surfacemodule(		char* module1,double* xmin,double* ymin,double* zmin,
							double* xmax,double* ymax,double* zmax,
							double* c1,double* c2,double* c3,	double* rmax,double* perf,
							int** i_epsmap,int* igrid,double* scale,double** i_scspos,
							double** i_scsnor,bool** i_idebmap,int** i_ibgp,int* inside,int* ibnum,int*
							natom, double** xn1,double** rad,int** d,int* maxbgp,double* probe,double*
							exrad,char** atinf,int** atsurf)
	{
	  
	  	  
	  char module[1024];
	  
	  printf("\n<<INFO>> Loading module %s....",module1);
	  
  	#if (defined LINUX_OS) || (defined MAC_OS)
	
	  void* tempPoint;
	  char* error;
	  void *handle;	 	  
	
	  #ifdef LINUX_OS
	    sprintf(module,"%s.so",module1);
	  #else
	    sprintf(module,"%s.dylib",module1);
	  #endif
		
	  handle = dlopen (module, RTLD_LAZY);
	  if (!handle) 
	  {	    
            fputs (dlerror(), stderr);
	    printf("\n\n<<ERROR>> Cannot load Surface module named %s\n\n",module);
            exit(1);
	  }

	  tempPoint = dlsym(handle, "epsmakemodule");
	  if ((error = dlerror()) != NULL)  
	  {
            fputs(error, stderr);
	    printf("\n\n<<ERROR>> Cannot load function epsmakemodule in plug-in named %s\n",module);
            exit(1);
	  }
	  printf("ok!\n");
	  	  
	  fptr epsmakemodule = (fptr)tempPoint;
	  
	  // propagates the called function
	  (*epsmakemodule)(*xmin,*ymin,*zmin,*xmax,*ymax,*zmax,*c1,*c2,*c3,*rmax,*perf,*i_epsmap,
			   *igrid,*scale,*i_scspos,*i_scsnor,*i_idebmap,*i_ibgp,
			   *inside,ibnum,*natom,*xn1,*rad,*d,*maxbgp,*probe,*exrad,*atinf,*atsurf);
        		
	  dlclose(handle);		

	  #elif defined WIN_OS
	  
		sprintf(module,"%s.dll",module1);	
		HMODULE delphi_surf = LoadLibraryEx(TEXT(module),NULL,LOAD_WITH_ALTERED_SEARCH_PATH);
		if( delphi_surf == NULL) 
		{
			printf("\n\n<<ERROR>> Cannot load Surface module named %s\n\n",module);
            exit(1);
		} 
		else 
		{
			fptr epsmakemodule = (fptr)GetProcAddress(delphi_surf,"epsmakemodule");
			if(epsmakemodule == NULL) 
			{			
				printf("\n\n<<ERROR>> Cannot load function epsmakemodule in plug-in named %s\n",module);
				exit(1);
			} 
			else 
			{								
				(*epsmakemodule)(*xmin,*ymin,*zmin,*xmax,*ymax,*zmax,*c1,*c2,*c3,*rmax,*perf,*i_epsmap,
			   *igrid,*scale,*i_scspos,*i_scsnor,*i_idebmap,*i_ibgp,
			   *inside,ibnum,*natom,*xn1,*rad,*d,*maxbgp,*probe,*exrad,*atinf,*atsurf);
			// ....it must not be freed because memory allocated inside the dll gets invalidated because
			// in Windows the dell memory has its own heap space...
			//	FreeLibrary(delphi_surf);				
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
