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
    # Calculate string distance to the input epitope
    # for all entries of epitopes withing the database

    # Note: good idea to optimize process by only calculating
    # string distance of non repeating epitopes
    def stringDistanceCalc(self):
        # ???
        return

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
    def findEpitopeMatch():
        return

a = epitopeMatching('ABC')
print(a.epitopeDatabase)