clear
mpicc src.c

for execution in 1 to 10
do
    bash NodeAllocator.sh                                       #Will check avilable nodes.
    for P in 4 16                                           #Number of Nodes
    do                          
        for ppn in 1 8                                      #(number of processes/cores per node)
        do
            python3 CreateHostFile.py $ppn                   #This will Create host file according to number of processes per nodes
            for D in 16 256 2048                             #(doubles)
            do
                for operation in 1 2 3 4                         # 1 for Broadcast, 2 for Reduce, 3 for Gather, 4 for AlltoAllv
                do
                    mpirun -np $((P * ppn)) -f hostfile ./a.out $D $operation
                done
            done
        done
    done
done

echo "Creating Plots"
python3 plot.py
echo "Plots Created"

