f = open("TemporaryFile.txt", "r")
arr = []; dist = {} ; temp = -1
for i in f:
    arr.append(str(i.replace("\n", "")))
for i in range(len(arr)-1, 0, -1):
    if arr[i].isnumeric():
        temp = int(arr[i])
        continue
    elif temp != -1:
        dist[arr[i]] = temp
        temp = -1

sorted_dict = {}
sorted_keys = sorted(dist, key=dist.get)  # [1, 3, 2]
f.close()

for w in sorted_keys:
    sorted_dict[w] = dist[w]

f = open("hosts", "w")
for i in sorted_keys:
    f.write(i + ":" + str(8) + "\n")


