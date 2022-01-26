#!/bin/bash
echo "Node Allocator starting"
echo '' > TemporaryFile.txt
for i in 2 3 5 6 7 8 9 10 11 12 14 15 16 31 13 17 18 19 20 21 22 23 24 25 26 27 28 29 30 32 33 34 35 36 37 38 39 40 41 42 43 44 46 45 47 48 49 50 51 52 53 54 56 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 
#for i in 2 3 4 5 6 7 8 9 10 11 12 14 15 16 31 13 17 18 19 20 21 22 23 24 25 26 27 28 29 30 32 33 34 35 36 37 38 39 40 41 42 43 44 46 45 47 48 49 50 51 52 53 54 56 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 
do 
echo 'csews'$i >> TemporaryFile.txt
echo $i
ssh csews$i uptime | awk '{print $10}' >> TemporaryFile.txt
#ssh csews$i uptime
done

python3 CreateHostFile.py 1
echo "Host File Created"


#6

