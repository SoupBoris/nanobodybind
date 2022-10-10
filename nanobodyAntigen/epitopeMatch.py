# load nanobodies-epitopes db into a pandas df
# compare string distance to input value creating news column 'distance'

# collect rows containing top 3 similar 'distance' entries (min)

# return pandas df containing N rows of nanobody info

import this
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

    # --- Set input funct ---
    # Utility function to set input manually
    def setInputEpitope(epitopeInput):
        epitope = epitopeInput

    

