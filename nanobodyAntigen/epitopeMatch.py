# load nanobodies-epitopes db into a pandas df
# compare string distance to input value creating news column 'distance'

# collect rows containing top 3 similar 'distance' entries (min)

# return pandas df containing N rows of nanobody info

from pandas import DataFrame
from pandas import read_csv
from Bio.SubsMat import MatrixInfo


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

        self.epitopeDatabase.drop(self.epitopeDatabase[self.epitopeDatabase['Epitope'] == 'Epitop'].index, inplace = True)
        self.epitopeDatabase.dropna(axis=0, inplace=True)
        return

    # --- Set input funct ---
    # Utility function to set input manually
    def setInputEpitope(self, epitopeInput):
       self.epitope = epitopeInput
       return

    def score_match(self, pair, matrix):
        if pair not in matrix:
            return matrix[(tuple(reversed(pair)))]
        else:
            return matrix[pair]

    # --- String distance calculator ---
    # Calculate string distance to the input epitope
    # for all entries of epitopes withing the database
    # Using blosum62 matrix algorithm

    # Note: good idea to optimize process by only calculating
    # string distance of non repeating epitopes
    def stringDistanceCalc(self, string1, string2):
        string1 = string1.strip()
        string1 = string1.replace(',','')
        string1 = string1.replace(' ','')
        string2 = string2.strip()
        string2 = string2.replace(',','')
        string2 = string2.replace(' ','')

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
                    myscore = self.score_match(pair, blosum)

                    pairi = (si, si)
                    myscorei = self.score_match(pairi, blosum)

                    pairj = (sj, sj)
                    myscorej = self.score_match(pairj, blosum)

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
        return

    # --- Find matching epitopes ---
    # Finds matching eptiopes, chooses top 3
    # and returns a list of all nanobodies associated
    # to these 3 eptiopes
    def findEpitopeMatch(self):
        self.epitopeDatabase['StringDistance'] = self.epitopeDatabase['Epitope'].apply(self.stringDistanceCalc, args = [self.epitope])
        self.epitopeDatabase.sort_values(axis=0, by='StringDistance', ascending=True, inplace=True)
        filteredEpitopes = self.epitopeDatabase.drop_duplicates(subset='Epitope', ignore_index=True)
        filteredEpitopes = filteredEpitopes.iloc[0:5,20].values
        # If you check the data on epitopes with least string-dist, 
        #these are the smallest ones, perhaps we should reconsider using a certain range of string distance rather 
        # #than a fixed amount of most similar entries

        filteredNanobodies = self.epitopeDatabase[self.epitopeDatabase['Epitope'].isin(filteredEpitopes)]
        print(filteredNanobodies)
        return

a = epitopeMatching('D,Y,G,NFNTQATNRNTDGSTDY,CNDGRTPGSR,VSDGNGM,CKGTD')

a.findEpitopeMatch()
#print(a.epitopeDatabase)
