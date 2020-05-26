#! /bin/bash -e

for i in {000..100}
do
    echo $i
    sed "s/XXX/$i/" parameters_G54a7.dat > parameters_${i}_G.dat
    python3 ../../../pypka.py parameters_${i}_G.dat
    exit
done
