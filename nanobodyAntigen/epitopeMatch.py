# load nanobodies-epitopes db into a pandas df
# compare string distance to input value creating news column 'distance'

# collect rows containing top 3 similar 'distance' entries (min)

# return pandas df containing N rows of nanobody info

from pandas import DataFrame
from pandas import read_csv


class epitopeMatching():

    # epitope argument, a string of the
    # target epitope
    epitope = ''

    # database of matching epitopes to
    # nanobodies that bind together
    epitopeDatabase = DataFrame()


    # --- Constructor funct ---
    # Creates epitope variable as well as
    # loading up database of matching
    # epitopes to nanobodies
    def  __init__(self, epitopeInput):
        self.epitope = epitopeInput
        self.epitopeDatabase = read_csv('nanobodyAntigen/databases/realCDRS.csv',encoding='cp1252')
        return

    # --- Set input funct ---
    # Utility function to set input manually
    def setInputEpitope(self, epitopeInput):
       self.epitope = epitopeInput
       return

    # --- String distance calculator ---
    #
    #
    def stringDistanceCalc(self):
        
    # use string Distance with blosum62
    def string_Distance(string1, string2):
        string1 = string1.strip()
        string2 = string2.strip()
        n = len(string1)
        m = len(string2)
        if (n == 0):
            return m
        if (m == 0):
            return n

        dmatrix = [[0 for x in range(m)] for y in range(n)]

        for i in range(n):
            dmatrix[i][0] = i
        for j in range(m):
         dmatrix[0][j] = j

        for i in range(1, n):
            si = string1[i - 1]
            for j in range(1, m):
                sj = string2[j - 1]
                if (si.upper() == sj.upper()):
                    cost = 0
                else:
                    blosum = MatrixInfo.blosum62
                    pair = (si, sj)
                    myscore = score_match(pair, blosum)

                    pairi = (si, si)
                    myscorei = score_match(pairi, blosum)

                    pairj = (sj, sj)
                    myscorej = score_match(pairj, blosum)

                    testcost = (myscorei + myscorej) / 2 - myscore
                    cost = testcost / 2

                dmatrix[i][j] = min(min(dmatrix[i - 1][j] + 3, dmatrix[i][j - 1] + 3), dmatrix[i - 1][j - 1] + cost)
        return dmatrix[n - 1][m - 1]
            

    # --- Retreiver of nanobodies ---
    # Retreives nanobodies of interest based
    # on their binding properties to the antigens
    # similar to the input antigen as a NEW df
    def retreiveNanobodies(self, matchingEpitopes): 
        # ???
        val = 5
        return


a = epitopeMatching('ABC')
print(a.epitopeDatabase)
