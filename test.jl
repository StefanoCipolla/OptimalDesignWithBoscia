using ODWB
seed = 1
m = 100
n = 10
criterion = "D"
time_limit = 3600 # 1 hour
corr = true # this create correlated data, if false, it will create independent data.

# Boscia
ODWB.solve_opt(seed, m, n, time_limit, criterion,corr)

# SCIP
#if criterion in ["A", "D"]
#    error("SCIP OA does not work with the optimal problems!")
#end
#ODWB.solve_opt_scip(seed, m, n, time_limit, criterion, corr)

# Pajarito
#ODWB.solve_opt_pajarito(seed, m, n, time_limit, criterion, corr)

# Custom
ODWB.solve_opt_custom(seed, m, n, time_limit, criterion, corr)
