import matplotlib.pyplot as plt


filename = "../visfiles/compress_M1000_N256.txt"
vz = 0.0001
dt = 0.005
area = 5184

displ = []
fzs   = []
with open(filename, 'r') as file:
    for line in file:
        line = file.readline()
        line = line.strip()
        print(line)
        arr = line.split()
        if len(arr) != 11 or not arr[0].isdigit(): continue
        d =  float(arr[0]) * vz * dt
        displ.append(d)
        fzs.append(float(arr[10])/area)

plt.plot(displ, fzs)
plt.show()
