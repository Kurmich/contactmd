import matplotlib.pyplot as plt


displ = []
fzs   = []
with open('conethermodata.txt', 'r') as file:
    for line in file:
        line = file.readline()
        line = line.strip()
        print(line)
        arr = line.split()
        if len(arr) != 11: continue
        displ.append(-float(arr[1]))
        fzs.append(float(arr[10]))

plt.plot(displ, fzs)
plt.show()
