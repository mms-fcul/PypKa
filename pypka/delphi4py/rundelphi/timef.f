        subroutine datime(day)
   
        character*24 day

#ifdef PC
cDEC$ IF DEFINED (PC)
        character*24 fdate
        external fdate

        day=fdate()
cDEC$ ELSE
#elif (AIX)
        external fdate_
        call fdate_(day)
#else
        character*24 fdate
#ifndef GNU
        external fdate
#endif
        day=fdate()
#endif
CDEC$ END IF

        end

        subroutine ddtime(tarray)
  
#ifdef PC
cDEC$ IF DEFINED (PC)
cDEC$ ELSE
#elif (AIX)
#define dtime dtime_
#endif
CDEC$ END IF

        real*4 tarray(2)
        real*4 dtime
#ifndef GNU
        external dtime
#endif

        sum=dtime(tarray)

        end

        subroutine eetime(tarray)

#ifdef PC
cDEC$ IF DEFINED (PC)
cDEC$ ELSE
#elif (AIX)
#define etime etime_
#endif
CDEC$ END IF
 
        real*4 tarray(2)
        real*4 etime
#ifndef GNU
        external etime
#endif

        sum=etime(tarray)

        end

