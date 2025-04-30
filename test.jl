using ODWB
seed = 1
m = 50
n = 20
criterion = "D"
time_limit = 36000 # 1 hour
corr     = true # this create correlated data, if false, it will create independent data.
long_run = false

A, C, N, ub, C_hat = ODWB.build_data_ext(seed, m, n, criterion in ["AF", "DF", "GTIF"], corr, scaling_C=long_run)

N = float(n)
@show ub
#solution_custom = ODWB.solve_opt_custom_ext(seed, m, n, 60, criterion, corr;
#                            A_external=A,
#                            C_external=C,
#                            N_external=N,
#                            ub_external=ub,
#                            C_hat_external=C_hat)
solution_boscia = ODWB.solve_opt_ext(seed, m, n, time_limit, criterion,corr;
                                 A_external=A,
                                 C_external=C,
                                 N_external=N,
                                 ub_external=ub,
                                 C_hat_external=C_hat)                       
file = ODWB.matopen("BosciaGenerated.mat", "w")
ODWB.write(file, "X", A)
ODWB.write(file, "N", N)
ODWB.close(file)