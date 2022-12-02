
doTimestamps = True
DO_T1_DIAGN = True
VERBOSITY = 5
GENERAL_BASES = 'D T Q 5'
GENERAL_TIME_LIMIT = '100'
GENERAL_RAM_MEMORY = '64'
GENERAL_CORES = '8'

MEMORY_TO_CLUSTER = {64: 'ulysses', 100: 'telemachus'} # keys represent RAM in gb. Any job that requires less than 100GB will be run on ulysses. More than 100 - telemachus.
MEMORY_TO_CORES = {64: 16, 100: 16} # how many cores to allocate with a given memory