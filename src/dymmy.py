class A:
    def __init__(self, a):
        self.a = a

l = []
v = A(5)
l.append(v)
d = set()
d.add((v,v))

for p in l: print(p.a)
for p in d: print(p[0].a)
v.a = 7
d.add((v,v))
for p in l: print(p.a)
for p in d: print(p[0].a)
