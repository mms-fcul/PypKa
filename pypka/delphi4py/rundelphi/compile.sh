### USER VARIABLES ###

precision=$1

recompile=$2

fortran_compiler=/usr/bin/gfortran
c_compiler=/usr/bin/gcc
cplus_compiler=/usr/bin/g++


### MAIN ###
if [[ ! $1 || ( $1 != "single" && $1 != "double" ) ]]
then
    echo Please define the precision argument as \"single\" or \"double\"
    exit
fi

if [[ ! $2 || ( $2 != "true" && $2 != "false" ) ]]
then
    echo Please define the recompile argument as \"true\" or \"false\"
    exit
fi

if [ $precision == "single" ]
then
    flags="-DGNU -fcray-pointer -ffixed-line-length-none -O3 -cpp"
elif [ $precision == "double" ]
then
    flags="-DGNU -fcray-pointer -ffixed-line-length-none -O3 -cpp -DDP -fdefault-real-8"
    flags="-g -O3 -fcray-pointer -ffixed-line-length-none -DGNU -DIFC -DDP -cpp -fdefault-real-8 -lstdc++"
fi

if [ $recompile == true ]
then    
    # remove all old shared objects
    rm -f rundelphi.so

    # f2py is not able to compile wrtsit4.f for some unknown reason
    # so wrtsit4 will also be recompiled
    gfortran $flags -I. -c  wrtsit4.f -shared -fPIC -o wrtsit4.so
else
    # remove all old shared objects
    rm -f *.so

    # c and c++ files have to be precompiled
    $c_compiler -O3  -DLINUX -c creapdb.c -shared -o creapdb.so -fPIC
    $c_compiler -O3  -DLINUX -c memalloc.c -shared -o memalloc.so -fPIC
    $cplus_compiler -O3  -DLINUX -c surfacemodule.cpp -shared -o surfacemodule.so -fPIC
    $cplus_compiler -O3  -DLINUX -c solvermodule.cpp -shared -o solvermodule.so -fPIC
    $cplus_compiler -O3  -DLINUX -c nlsolvermodule.cpp -shared -o nlsolvermodule.so -fPIC
    
    for i in `ls *.f`
    do
	gfortran $flags -I. -c  $i -shared -fPIC -o ${i:0:-2}.so
    done
fi

# f2py command
# only compiler tested was gfortran
f2py2.7 -c -m rundelphi creapdb.so memalloc.so surfacemodule.so solvermodule.so nlsolvermodule.so wrtsit4.so delphi.f elb.f namlen3.f qinttot.f cputime.f timef.f up.f wrt4.f rent4.f rdhrad.f ichash4.f rdhcrg.f cent4.f setrc4d.f getatm2.f radass4.f crgass4.f	chkcrg.f form2.f rfind4.f cfind4.f extrmobjects.f extrm.f off4.f wrtprm.f grdatm.f crgarr.f distrTOpoint.f scaler4.f epsmak.f setout.f vwtms2.f distobj.f sas.f cube.f indver.f watput.f scale.f msrf.f ex.f mkvtl.f fxhl.f wrtsurf.f objvertices.f wrtatd.f dbsfd.f mkdbsf.f wrteps4b.f setfcrg.f setcrg.f setbc4.f relfac4b.f wrtgcrg.f phintp4.f nitit.f conplt.f itit4j.f omalt.f encalc4.f nlener.f anagrd4.f react2.f clb.f rdiv.f watpte.f rforce.f ts.f phicon.f expand4.f wrtphi.f --f77exec="$fortran_compiler" --f77flags="$flags"

ln -s ../readFiles/readFiles.so .

