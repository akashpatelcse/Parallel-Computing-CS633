#!/bin/bash
PATH=$HOME/mpich-install/bin:$PATH;        
export PATH;
which mpirun

PATH=$HOME/Sem2/mpich-install/bin:$PATH;
export PATH;
which mpirun
mpicc src.c                     #Compilation

for x in 0 1 2 3 4 
do
    bash NodeAllocatorLite.sh
    for mode in 1 2 3           # 1 = Normal Mode, 2 = Vector Mode, 3 = Pack Mode
    do
        for P in 16 36 49 64
        do
            for N in 16 32 64 128 512 1024
            do  
                mpirun -n $P -f hosts ./a.out $N $mode
            done
        done
    done
done

python3 PlotGraphs.py           #Plot Graph
echo "Plot Created"

