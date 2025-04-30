
"""
build_data

seed - for the Random functions.
m    - number of experiments.
fusion - boolean deiciding whether we build the fusion or standard problem.
corr - boolean deciding whether we build the independent or correlated data.   
"""
function build_data_ext(seed, m, n, fusion, corr; scaling_C=false)
    @show scaling_C
    # set up
    Random.seed!(seed)
    if corr 
        B = rand(m,n)
        B = B'*B
        @assert isposdef(B)
        D = MvNormal(randn(n),B)
        
        A = rand(D, m)'
        @assert rank(A) == n 
    else 
        A = rand(m,n)
        @assert rank(A) == n # check that A has the desired rank!
    end 
    C_hat = rand(2n, n)
    C = scaling_C ? 1/2n*transpose(C_hat)*C_hat : transpose(C_hat)*C_hat

    @assert rank(C) == n
    
    if fusion
        N = rand(floor(m/20):floor(m/3))
        ub = rand(1.0:m/10, m)
    else
        N = floor(1.5*n)
        #u = floor(N/3)
        #ub = rand(1.0:u, m)
        ub = N.* ones(m)
        #ub  = fill(Inf,m)
    end
        
    return A, C, N, ub, C_hat
end