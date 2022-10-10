from CDR1 import finalmat
from CDR2 import finalmat2
from CDR3 import finalmat3

CombinationsofCDRs = []

for ay in range(5000):
    for aydios in range(7680):
        for aydiosmio in range(3125000):
            CombinationsofCDRs.append(f'{finalmat[ay]}{finalmat2[aydios]}{finalmat3[aydiosmio]}')
            print(aydiosmio)
    print(ay)

print(CombinationsofCDRs)