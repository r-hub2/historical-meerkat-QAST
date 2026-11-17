# Variances (jackknife/bootstrap) in the main process during the tests:
# covr only records the coverage of the main process, and the refits run
# in the fork children when there is more than one worker. The fork path
# is still tested explicitly in test-jackknife-se.R (cores = 2).
options(QAST.cores = 1L)
