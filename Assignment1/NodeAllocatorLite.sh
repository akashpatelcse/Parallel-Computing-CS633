#!/bin/bash
echo "Node Allocator starting"
echo '' > TemporaryFile.txt
for i in {1..40}
do 
echo 'csews'$i >> TemporaryFile.txt
echo $i
ssh csews$i uptime | awk '{print $6}' >> TemporaryFile.txt
done

python3 CreateHostFile.py
echo "Host File Created"
