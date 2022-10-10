from CDR1 import finalmat
from CDR2 import finalmat2
from CDR3 import finalmat3

CombinationsofCDRs = [None] * 120000000000000

for ay in range(5000):
    for aydios in range(7680):
        for aydiosmio in range(3125000):
            CombinationsofCDRs[aydiosmio+(aydios*50000)+(aydiosmio*24000000000)]=finalmat[ay]+finalmat2[aydios]+finalmat3[aydiosmio]

print(CombinationsofCDRs)