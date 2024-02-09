# Test Fusion problem with Co-BnB
using optDesign
using Random

# Fusion Co-BNB

for criterion in ["AF","DF"]
    for m in [50] # 60,80
        for k in [10,4]
            n = Int(floor(m/k))
            for seed in 1:2
                @show m ,n, seed
                optDesign.solve_opt_custom(seed, m, n, 300, criterion, false, write=false,verbose=true)
                optDesign.solve_opt(seed,m,n,300,criterion,false,write=false,verbose=false)
            end
        end
    end
end


# Optimal with SCIP
#=
for criterion in ["A","D"]
    for m in [50,60]
        for k in [4,10]
            n = Int(floor(m/k))
            for seed in 1:2
                @show m ,n, seed
                optDesign.solve_opt_scip(seed, m, n, 1800, criterion, false, write=false, verbose=true)
            end
        end
    end
end 
=#