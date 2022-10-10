#Starting sequence for mutations
CDR3 = "AAYYYYYYYFYY"
lenseq1 = len(CDR3)

#A vector of the length of the sequence I wish to create all mutations of
mutsinCDR3 = [None] * 3125000

cdr3muts1 = ("A", "V")
cdr3muts2 = ("Y", "G", "S", "D", "T")
cdr3muts3 = ("F", "H", "L", "Y")


for j2 in range(625000):
    for j3 in range(5):
        mutsinCDR3[j3+j2*5] = CDR3[0:2] + cdr3muts2[j3] + CDR3[3:]


for j1 in range(3125000):
    if j1<1562500:
        mutsinCDR3[j1] = mutsinCDR3[j1][0] + cdr3muts1[0] + mutsinCDR3[j1][2:]
    elif j1>1562499:
        mutsinCDR3[j1] = mutsinCDR3[j1][0] + cdr3muts1[1] + mutsinCDR3[j1][2:]


cdr3dummy2=[]
for j4 in range(125000):
    for j5 in range(5):
        for j6 in range(5):
            mutsinCDR3[j6+(j5*5)+(j4*25)] = mutsinCDR3[j6+(j5*5)+(j4*25)][0:3] + cdr3muts2[j5] + \
                                            mutsinCDR3[j6+(j5*5)+(j4*25)][4:]
            cdr3dummy2.append(mutsinCDR3[j6+(j5*5)+(j4*25)])


cdr3dummy3=[]
for j7 in range(25000):
    for j8 in range(5):
        for j9 in range(25):
            cdr3dummy2[j9+(j8*25)+(j7*125)] = cdr3dummy2[j9+(j8*25)+(j7*125)][0:4] + cdr3muts2[j8] + \
                                                      cdr3dummy2[j9+(j8*25)+(j7*125)][5:]
            cdr3dummy3.append(cdr3dummy2[j9 + (j8 * 24) + (j7 * 120)])


cdr3dummy4=[]
for j10 in range(5000):
    for j11 in range(5):
        for j12 in range(125):
            cdr3dummy3[j12+(j11*125)+(j10*625)] = cdr3dummy3[j12+(j11*125)+(j10*625)][0:5] + cdr3muts2[j11] \
                                                          + cdr3dummy3[j12 + (j11 * 120) + (j10 * 480)][6:]
            cdr3dummy4.append(mutsinCDR3[j12 + (j11 * 120) + (j10 * 480)])


cdr3dummy5=[]
for j13 in range(1000):
    for j14 in range(5):
        for j15 in range(625):
            cdr3dummy4[j15 + (j14*625) + (j13 * 3125)] = cdr3dummy4[j15 + (j14 * 625) + (j13 * 3125)][0:6] + cdr3muts2[j14] \
                                                         + cdr3dummy4[j15 + (j14*625) + (j13 * 3125)][7:]
            cdr3dummy5.append(cdr3dummy4[j15 + (j14*625) + (j13 * 3125)])


cdr3dummy6=[]
for j16 in range(200):
    for j17 in range(5):
        for j18 in range(3125):
            cdr3dummy5[j18+(j17*3125)+(j16*15625)] = cdr3dummy5[j18+(j17*3125)+(j16*15625)][0:7] + cdr3muts2[j17] \
                                                         + cdr3dummy5[j18+(j17*3125)+(j16*15625)][8:]
            cdr3dummy6.append(cdr3dummy5[j18+(j17*3125)+(j16*15625)])


cdr3dummy7=[]
for j19 in range(40):
    for j20 in range(5):
        for j21 in range(15625):
            cdr3dummy6[j21+(j20*15625)+(j19*78125)] = cdr3dummy6[j21+(j20*15625)+(j19*78125)][0:8] + cdr3muts2[j20] \
                                                         + cdr3dummy6[j21+(j20*15625)+(j19*78125)][9:]
            cdr3dummy7.append(cdr3dummy6[j21+(j20*15625)+(j19*78125)])


cdr3dummy8=[]
for j22 in range(10):
    for j23 in range(4):
        for j24 in range(78125):
            cdr3dummy7[j24+(j23*78125)+(j22*312500)] = cdr3dummy7[j24+(j23*78125)+(j22*312500)][0:9] + cdr3muts3[j23] \
                                                             + cdr3dummy7[j24+(j23*78125)+(j22*312500)][10:]
            cdr3dummy8.append(cdr3dummy7[j24+(j23*78125)+(j22*312500)])

finalmat3=[]
for j25 in range(2):
    for j26 in range(5):
        for j27 in range(312500):
            cdr3dummy8[j27+(j26*312500)+(j25*1562500)] = cdr3dummy8[j27+(j26*312500)+(j25*1562500)][0:10] + cdr3muts2[j26] \
                                                             + cdr3dummy8[j27+(j26*312500)+(j25*1562500)][11:]
            finalmat3.append(cdr3dummy8[j27+(j26*312500)+(j25*1562500)])

print(finalmat3)