## Script for running sbatch

### For the solver comparison

# Instance settings
num_experiments = [50,60,80,100,120]
criteria = ["D", "A", "DF", "AF"]
data_types = ["CORR", "IND"]
solvers = ["Boscia", "Pajarito", "SCIP", "Custom"]

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
