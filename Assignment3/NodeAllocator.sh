#!/bin/bash
echo "Node Allocator starting"
echo '' > TemporaryFile.txt
for i in 11 25 41 54 70 7 9 13
do 
echo 'csews'$i >> TemporaryFile.txt
echo $i
ssh csews$i uptime | awk '{print $10}' >> TemporaryFile.txt
done

python3 CreateHostFile.py 1
echo "Host File Created"


