import sys

print (sys.argv[1])
process_per_node = int(sys.argv[1])

f = open("TemporaryFile.txt", "r")
arr = []; dist = {} ; temp = -1
for i in f:
    arr.append(str(i.replace("\n", "")).replace(",",""))
for i in range(len(arr)-1, 0, -1):
    ttt = arr[i]
    if ttt.replace(".", "").isnumeric():
        temp = float(arr[i])
        continue
    elif temp != -1:
        dist[arr[i]] = temp
        temp = -1

sorted_dict = {}
sorted_keys = sorted(dist, key=dist.get)  # [1, 3, 2]
f.close()

for w in sorted_keys:
    sorted_dict[w] = dist[w]

f = open("hostfile", "w")
for i in sorted_keys:
    f.write(i + ":" + str(process_per_node) + "\n")



