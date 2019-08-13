### USER VARIABLES ###

precision=$1
compiler=/usr/bin/gfortran

### MAIN ###
if [[ ! $1 || ( $1 != "single" && $1 != "double" ) ]]
then
    echo Please define the precision argument as \"single\" or \"double\"
    exit
fi

# remove old python module file
rm -f *.so

# c and c++ files have to be precompiled
gcc -O3  -DLINUX -c memalloc.c -shared -o memalloc.so -fPIC


if [ $precision == single ]
then
    flags="-DGNU -fcray-pointer -ffixed-line-length-none -O3 -cpp"
elif  [ $precision == double ]
then
    flags="-DGNU -fcray-pointer -ffixed-line-length-none -O3 -cpp -DDP -fdefault-real-8 "
fi


# f2py command
f2py3.6 -c -m readFiles memalloc.so delphi1.f qinttot.f namlen3.f up.f rdhrad.f elb.f rent4.f rdhcrg.f cent4.f setrc4d.f getatm2.f form2.f radass4.f rfind4.f crgass4.f cfind4.f chkcrg.f ichash4.f extrm.f off4.f --f77exec="$compiler" --f77flags="$flags"

mv readFiles*.so readFiles.so
