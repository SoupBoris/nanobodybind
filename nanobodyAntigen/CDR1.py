#Starting sequence for mutations
CDR1 = "GNIFYYYYM"
lenseq1 = len(CDR1)

#A vector of the length of the sequence I wish to create all mutations of
mutsinCDR1 = [None] * 5000

cdr1mutations1 = ("N", "S", "T", "Y")
cdr1mutations2 = ("F", "S")
cdr1mutations3 = ("Y", "G", "S", "D", "T")

for t in range(1250):
    for i in range(4):
        mutsinCDR1[i] = CDR1[0] + cdr1mutations1[i] + CDR1[2:]
        mutsinCDR1[i+t*4] = mutsinCDR1[i]

for j in range(5000):
    if j<2500:
        mutsinCDR1[j] = mutsinCDR1[j][0:3] + cdr1mutations2[0] + mutsinCDR1[j][4:]
    elif j>2499:
        mutsinCDR1[j] = mutsinCDR1[j][0:3] + cdr1mutations2[1] + mutsinCDR1[j][4:]

dummymat=[]
for m in range(250):
    for n in range(5):
        for l in range(4):
            mutsinCDR1[l+n*4+(m*20)] = mutsinCDR1[l+n*4+(m*20)][0:4] + cdr1mutations3[n] + mutsinCDR1[l + n * 4 + (m * 20)][5:]
            dummymat.append(mutsinCDR1[l + (n * 4) + (m * 20)])

dummy2=[]
for p in range(50):
    for q in range(5):
        for r in range(20):
            dummymat[r+(q*20)+(p*100)] = dummymat[r+(q*20)+(p*100)][0:5] + cdr1mutations3[q] + dummymat[r + (q * 20) + (p * 100)][6:]
            dummy2.append(dummymat[r+(q*20)+(p*100)])

dummy3=[]
for g in range(10):
    for h in range(5):
        for k in range(100):
            dummy2[k+(h*100)+(g*500)] = dummy2[k+(h*100)+(g*500)][0:6] + cdr1mutations3[h] + dummy2[k + (h * 100) + (g * 500)][7:]
            dummy3.append(dummy2[k+(h*100)+(g*500)])

finalmat=[]
for d in range(2):
    for w in range(5):
        for z in range(500):
            dummy3[z+(w*500)+(d*2500)] = dummy3[z+(w*500)+(d*2500)][0:7] + cdr1mutations3[w] + dummy3[z + (w * 500) + (d * 2500)][8:]
            finalmat.append(dummy3[z+(w*500)+(d*2500)])

print(finalmat)