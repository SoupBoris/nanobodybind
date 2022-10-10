#Starting sequence for mutations
CDR2 = "EFVAAIAYGAITNY"
lenseq1 = len(CDR2)

#A vector of the length of the sequence I wish to create all mutations of
mutsinCDR2 = [None] * 7680

cdr2mutations1 = ("F", "L")
cdr2mutations2 = ("A", "G", "S", "T")
cdr2mutations3 = ("A", "D", "G", "N", "S", "T")
cdr2mutations4 = ("Y", "G", "S", "D", "T")
cdr2mutations5 = ("A", "G", "S", "T")
cdr2mutations6 = ("I", "N", "S", "T")
cdr2mutations7 = ("N", "Y")


for it2 in range(1920):
    for it3 in range(4):
        mutsinCDR2[it3+it2*4] = CDR2[0:4] + cdr2mutations2[it3] + CDR2[5:]


for it1 in range(7680):
    if it1<3840:
        mutsinCDR2[it1] = mutsinCDR2[it1][0] + cdr2mutations1[0] + mutsinCDR2[it1][2:]
    elif it1>3839:
        mutsinCDR2[it1] = mutsinCDR2[it1][0] + cdr2mutations1[1] + mutsinCDR2[it1][2:]


cdr2dummy1=[]
for it4 in range(320):
    for it5 in range(6):
        for it6 in range(4):
            mutsinCDR2[it6+(it5*4)+(it4*24)] = mutsinCDR2[it6+(it5*4)+(it4*24)][0:6] + cdr2mutations3[it5] + \
                                               mutsinCDR2[it6+(it5*4)+(it4*24)][7:]
            cdr2dummy1.append(mutsinCDR2[it6+(it5*4)+(it4*24)])


cdr2dummy2=[]
for it7 in range(64):
    for it8 in range(5):
        for it9 in range(24):
            cdr2dummy1[it9+(it8*24)+(it7*120)] = cdr2dummy1[it9+(it8*24)+(it7*120)][0:7] + cdr2mutations4[it8] +\
                                                 cdr2dummy1[it9+(it8*24)+(it7*120)][8:]
            cdr2dummy2.append(cdr2dummy1[it9+(it8*24)+(it7*120)])


cdr2dummy3=[]
for it10 in range(16):
    for it11 in range(4):
        for it12 in range(120):
            cdr2dummy2[it12+(it11*120)+(it10*480)] = cdr2dummy2[it12+(it11*120)+(it10*480)][0:9] + cdr2mutations5[it11]\
                                                     + cdr2dummy2[it12+(it11*120)+(it10*480)][10:]
            cdr2dummy3.append(cdr2dummy2[it12+(it11*120)+(it10*480)])


cdr2dummy4=[]
for it13 in range(4):
    for it14 in range(4):
        for it15 in range(480):
            cdr2dummy3[it15+(it14*480)+(it13*1920)] = cdr2dummy3[it15+(it14*480)+(it13*1920)][0:10] + \
                                                      cdr2mutations6[it14] \
                                                      + cdr2dummy3[it15+(it14*480)+(it13*1920)][11:]
            cdr2dummy4.append(cdr2dummy3[it15+(it14*480)+(it13*1920)])


finalmat2=[]
for it16 in range(2):
    for it17 in range(2):
        for it18 in range(1920):
            cdr2dummy4[it18+(it17*1920)+(it16*3840)] = cdr2dummy4[it18+(it17*1920)+(it16*3840)][0:12] + \
                                                       cdr2mutations7[it17] \
                                                       + cdr2dummy4[it18+(it17*1920)+(it16*3840)][13:]
            finalmat2.append(cdr2dummy4[it18+(it17*1920)+(it16*3840)])

print(finalmat2)