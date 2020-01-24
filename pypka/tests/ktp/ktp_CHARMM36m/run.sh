#! /bin/bash -e

for i in {000..100}
do
    echo $i
    sed "s/XXX/$i/" parameters.dat > parameters_$i.dat   
    python3 ../../../pypka.py parameters_$i.dat
    exit
done
