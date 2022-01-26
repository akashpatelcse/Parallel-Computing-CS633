clear
mpicc -o code src.c

for execution in 1 2 3 4 5
do
    bash NodeAllocator.sh                                       #Will check avilable nodes.
    for P in 1 2                                           #Number of Nodes
    do                          
        for ppn in 1 2 4                            #(number of processes/cores per node)
        do
            echo -n "$P $ppn "  | tee -a Temp_output.txt
            python3 CreateHostFile.py $ppn                   #This will Create host file according to number of processes per nodes
            mpirun -np $((P * ppn)) -f hostfile ./code tdata.csv
        done
    done
done

echo "Creating Plots"
python plot.py
echo "Plots Created"

