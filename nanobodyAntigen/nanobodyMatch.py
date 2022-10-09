# open big DB of nanobodies

# convert N rows pandas df of nanobodies to numpy array
# for each row in input_array compute string distance for H_CDR1-3 and do a weighted sum
# retreive top 10 minimum distance values from the existing nanobodies

# we now have a N*10 rows array of nanobodies
# remove duplicate results
# potentiallyBindNano -> an array with all the suspect nanobodies that could bind

# load up trained K-N and RF models

# compare CDRsm, pL and Hydrophobia using K-N and RF from potentiallyBindNano
# compute binding probabilty for all datapoints

# return consolidate dataframe of nanobodies with possibility to bind to antigen