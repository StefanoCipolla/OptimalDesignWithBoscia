## Script for running sbatch

### For the solver comparison

# Instance settings
num_experiments = [50,60,80,100,120]
#num_experiments = [120]
#criteria = ["D", "A", "DF", "AF"]
criteria = ["DF","AF"]
data_types = ["CORR", "IND"]
#data_types = ["CORR"]
#solvers = ["Boscia", "Pajarito", "SCIP", "Custom"]
solvers = ["Custom"]

# for exclusive runs set in the sbatch file: #SBATCH --exclusive

# create instances
for criterion in criteria
    for data in data_types
        for solver in solvers
            for m in num_experiments
                # 
                run(`sbatch -A optimi -J Co-Fusion experiment.sbatch $criterion $solver $data $m`) # CB
            end
        end
    end
end



### For the Boscia settings
#=
# Instance settings
#num_experiments = [50,80,120]
num_experiments = [120]
#num_experiments = [300, 400]
#criteria = ["D", "A", "DF", "AF"]
criteria=["A"]
#data_types = ["CORR", "IND"]
data_types = ["CORR"]
options=["shadow_set"]
#options = ["default", "MIP_SCIP", "heuristics"]
#options = ["strong_branching", "tightening", "shadow_set", "tighten_lazification"]
#options=["long_run"]

# for exclusive runs set in the sbatch file: #SBATCH --exclusive

# create instances
for criterion in criteria
    for data in data_types
        for option in options
            for m in num_experiments
                run(`sbatch -A optimi -J Bo-Option settings.sbatch $criterion $data $m $option`) # CB
            end
        end
    end
end
=#


# to cancel jobs
## scancel --user dhendryc
