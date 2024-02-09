# Build tables (csv files) comparing results from the different solvers
using CSV
using DataFrames

function build_non_grouped_csv(corr = true)
    modes = ["boscia", "FrankWolfe", "scip", "hypatia", "pajarito"]
    type = corr ? "correlated" : "independent"

    # prelimenary
    #dimensions = [20,30,50,60,80,100,120,150,180,200]
    #fracs = vcat(fill(4,5), fill(6,5), fill(8,5), fill(10,5))
    # exclusive
    dimensions = [50,60,80,100,120]
    fracs = vcat(fill(4,5), fill(10,5))

    #criteria = corr ? ["DF", "AF", "D"] : ["A", "AF", "D", "DF"] #  "DF" independent -> missing SCIP entries, "A" correlated missing Boscia entries
    criteria = ["A","AF","D","DF"]
    #criteria = ["A"] 
    for criterion in criteria 
        @show criterion
        println("\n")
        df = DataFrame()

        df_condi = if criterion in ["AF","DF"]
            DataFrame(CSV.File(joinpath(@__DIR__, "csv/fusion_" * type * "_data.csv")))
        else
            DataFrame(CSV.File(joinpath(@__DIR__, "csv/opt_" * type * "_data.csv")))
        end

        timeBoscia = []
        solutionBoscia = []
        terminationBoscia = []
        dualGapBoscia = []
        relDualGapBoscia = []
        lbBoscia = []
        numberNodesBoscia = []

        timeCustomBB = []
        solutionCustomBB = []
        terminationCustomBB = []
        dualGapCustomBB = []
        relDualGapCustomBB = []
        numberNodesCustomBB = []

        timeScip = []
        solutionScip = []
        terminationScip = []
        dualGapScip = []
        relDualGapScip = []
        numberCutsScip = []

        timePajarito = []
        solutionPajarito = []
        terminationPajarito = []
        dualGapPajarito = []
        relDualGapPajarito = []
        numberCutsPajarito = []
        numberIterPajarito = []

        for m in dimensions
            @show m
            df_dim = DataFrame()

            # load boscia
            println("Load boscia")
            df_boscia = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Boscia/boscia_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv")))

            df_boscia.termination .= replace.(df_boscia.termination, "Time limit reached" => "TIME_LIMIT")
            for row in eachrow(df_boscia)
                if row.time > 3600
                    row.termination = "TIME_LIMIT" 
                end
            end
            println("Termination boscia")
            termination_boscia = [row == "OPTIMAL" || row == "tree.lb>primal-dual_gap" || row == "primal>=tree.incumbent" ? 1 : 0 for row in df_boscia[!,:termination]]

            println("Set up data")
            df_dim[!,:seed] = df_boscia[!,:seed]
            df_dim[!,:numberOfExperiments] = df_boscia[!,:numberOfExperiments]
            df_dim[!,:numberOfParameters] = df_boscia[!, :numberOfParameters]
            df_dim[!,:numberOfAllowedEx] = df_boscia[!,:N]
            df_dim[!,:frac] = fracs
            

            println("Reading boscia data")
            df_dim[!,:timeBoscia] = df_boscia[!,:time]
            timeBoscia = vcat(timeBoscia, df_boscia[!,:time])
            df_dim[!,:solutionBoscia] = df_boscia[!,:scaled_solution]
            solutionBoscia = vcat(solutionBoscia, df_boscia[!,:scaled_solution])
            df_dim[!,:terminationBoscia] = termination_boscia
            terminationBoscia=vcat(terminationBoscia, termination_boscia)
            lowerBounds = m * (df_boscia[!, :solution] - df_boscia[!, :dual_gap])
            df_dim[!,:lbBoscia] = lowerBounds
            lbBoscia = vcat(lbBoscia, lowerBounds)
            df_dim[!, :dualGapBoscia] = m * df_boscia[!, :dual_gap]
            dualGapBoscia = vcat(dualGapBoscia, df_dim[!,:dualGapBoscia])
            df_dim[!,:numberNodesBoscia] = df_boscia[!,:num_nodes]
            numberNodesBoscia = vcat(numberNodesBoscia, df_boscia[!,:num_nodes])

            @show size(df_dim)
        
            # load scip 
            # only works with the fusion problems
            if criterion == "DF" || criterion == "AF"
                df_scip = DataFrame(CSV.File(joinpath(@__DIR__, "csv/SCIP/scip_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv")))
                df_scip.termination .= replace.(df_scip.termination, "Time limit reached" => "TIME_LIMIT")
                for row in eachrow(df_scip)
                    if row.time > 3600
                        row.termination = "TIME_LIMIT" 
                    end
                end

                termination_scip = [row == "OPTIMAL" ? 1 : 0 for row in df_scip[!,:termination]]

                println("Load SCIP")
                df_scip[!,:timeScip] = df_scip[!,:time]
                timeScip = vcat(timeScip, df_scip[!,:time])
                df_scip[!,:terminationScip] = termination_scip
                terminationScip = vcat(terminationScip, termination_scip)
                df_scip[!,:solutionScip] = df_scip[!,:solution]*m
                solutionScip = vcat(solutionScip, df_scip[!,:solution]*m)
                df_scip[!,:dualGapScip] = df_scip[!, :solution]*m - lowerBounds
                dualGapScip = vcat(dualGapScip, df_scip[!,:dualGapScip])
                df_scip[!,:numberCutsScip] = df_scip[!,:calls]
                numberCutsScip = vcat(numberCutsScip, df_scip[!,:calls])

                df_scip[!,:frac] = fracs
                df_scip[!,:numberOfAllowedEx] = df_scip[!,:N]
                df_scip = select(df_scip, [:terminationScip, :timeScip, :solutionScip, :dualGapScip, :numberCutsScip, :seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])
                df_dim = outerjoin(df_dim, df_scip, on = [:seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])

                @show size(df_dim)
         #=   else
                df_bb = DataFrame(CSV.File(joinpath(@__DIR__, "csv/CustomBB/customBB_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv")))
                df_bb.termination .= replace.(df_bb.termination, "Time limit reached" => "TIME_LIMIT")
                for row in eachrow(df_bb)
                    if row.time > 3600
                        row.termination = "TIME_LIMIT" 
                    end
                end

                termination_bb = [row == "optimal" ? 1 : 0 for row in df_bb[!,:termination]]

                println("Load Custom BB")
                df_bb[!,:timeCustomBB] = df_bb[!,:time]
                timeCustomBB = vcat(timeCustomBB, df_bb[!,:time])
                df_bb[!,:terminationCustomBB] = termination_bb
                terminationCustomBB = vcat(terminationCustomBB, termination_bb)
                df_bb[!,:solutionCustomBB] = df_bb[!,:solution_scaled]
                solutionCustomBB = vcat(solutionCustomBB, df_bb[!,:solution_scaled])
                df_bb[!,:dualGapCustomBB] = df_bb[!, :solution_scaled] - lowerBounds
                dualGapCustomBB = vcat(dualGapCustomBB, df_bb[!,:dualGapCustomBB])
                df_bb[!,:numberNodesCustomBB] = df_bb[!,:number_nodes]
                numberNodesCustomBB = vcat(numberNodesCustomBB, df_bb[!,:number_nodes])
        
                df_bb[!,:frac] = fracs
                df_bb[!,:numberOfAllowedEx] = df_bb[!,:N]
                df_bb = select(df_bb, [:terminationCustomBB, :timeCustomBB, :solutionCustomBB, :dualGapCustomBB, :numberNodesCustomBB, :seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])
                @show size(df_bb)
                df_dim = outerjoin(df_dim, df_bb, on = [:seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])

                @show size(df_dim)=#
            end

            df_bb = DataFrame(CSV.File(joinpath(@__DIR__, "csv/CustomBB/customBB_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv")))
                df_bb.termination .= replace.(df_bb.termination, "Time limit reached" => "TIME_LIMIT")
                for row in eachrow(df_bb)
                    if row.time > 3600
                        row.termination = "TIME_LIMIT" 
                    end
                end

                termination_bb = [row == "optimal" ? 1 : 0 for row in df_bb[!,:termination]]

                println("Load Custom BB")
                df_bb[!,:timeCustomBB] = df_bb[!,:time]
                timeCustomBB = vcat(timeCustomBB, df_bb[!,:time])
                df_bb[!,:terminationCustomBB] = termination_bb
                terminationCustomBB = vcat(terminationCustomBB, termination_bb)
                df_bb[!,:solutionCustomBB] = df_bb[!,:solution_scaled]
                solutionCustomBB = vcat(solutionCustomBB, df_bb[!,:solution_scaled])
                df_bb[!,:dualGapCustomBB] = df_bb[!, :solution_scaled] - lowerBounds
                dualGapCustomBB = vcat(dualGapCustomBB, df_bb[!,:dualGapCustomBB])
                df_bb[!,:numberNodesCustomBB] = df_bb[!,:number_nodes]
                numberNodesCustomBB = vcat(numberNodesCustomBB, df_bb[!,:number_nodes])
        
                df_bb[!,:frac] = fracs
                df_bb[!,:numberOfAllowedEx] = df_bb[!,:N]
                df_bb = select(df_bb, [:terminationCustomBB, :timeCustomBB, :solutionCustomBB, :dualGapCustomBB, :numberNodesCustomBB, :seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])
                @show size(df_bb)
                df_dim = outerjoin(df_dim, df_bb, on = [:seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])

                @show size(df_dim)

        

            # load pajarito
            println("Load Pajarito")
            df_paj = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Pajarito/pajarito_" * criterion * "_" * string(m) * "_" * type * "_optimality.csv")))

            for row in eachrow(df_paj)
                if row.time > 3600
                    row.termination = "TIME_LIMIT" 
                end
            end

            termination_paj = [row == "OPTIMAL" ? 1 : 0 for row in df_paj[!,:termination]]

            time_paj = df_paj[!,:time]
            for i in eachindex(time_paj)
                if df_paj[i,:termination] == "OTHER_ERROR"
                    time_paj[i] = 3600.00
                end
            end

            df_paj[!,:timePajarito] = time_paj
            timePajarito = vcat(timePajarito, time_paj)
            df_paj[!,:solutionPajarito] = df_paj[!,:solution]
            solutionPajarito = vcat(solutionPajarito, df_paj[!,:solution])
            df_paj[!,:terminationPajarito] = termination_paj
            terminationPajarito = vcat(terminationPajarito, termination_paj)
            df_paj[!,:dualGapPajarito] = df_paj[!, :solution] - lowerBounds
            dualGapPajarito = vcat(dualGapPajarito, df_paj[!,:dualGapPajarito])
            df_paj[!,:numberCutsPajarito] = df_paj[!,:numberCuts]
            numberCutsPajarito = vcat(numberCutsPajarito, df_paj[!,:numberCuts])
            df_paj[!,:numberIterPajarito] = df_paj[!,:numberIterations]
            numberIterPajarito = vcat(numberIterPajarito, df_paj[!,:numberIterations])

            df_paj[!,:frac] = fracs
            df_paj[!,:numberOfAllowedEx] = df_paj[!,:N]
            df_paj = select(df_paj, [:terminationPajarito, :timePajarito, :solutionPajarito, :dualGapPajarito, :numberCutsPajarito, :numberIterPajarito, :seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])
            df_dim = outerjoin(df_dim, df_paj, on = [:seed, :numberOfExperiments, :numberOfParameters, :numberOfAllowedEx, :frac])

            @show size(df_dim)
        
            # compute relative gap
            rel_gap_scip = []
            rel_gap_boscia = []
            rel_gap_pajarito = []
            rel_gap_custombb = []
            #rel_gap_hypatia_cont = []
            #rel_gap_hypatia_limit = []
            #rel_gap_frankwolfe_cont = []
            #rel_gap_frankwolfe_limit = []
            for row in eachrow(df_dim)
                if criterion == "AF" || criterion == "DF"
                    if min(abs(row.solutionScip), abs(row.lbBoscia)) == 0
                        push!(rel_gap_scip, row.solutionScip - row.lbBoscia)
                    elseif sign(row.lbBoscia) != sign(row.solutionScip)
                        push!(rel_gap_scip, Inf)
                    else
                        push!(rel_gap_scip, (row.solutionScip - row.lbBoscia)/min(abs(row.solutionScip), abs(row.lbBoscia)))
                    end
             #   else
             #       if min(abs(row.solutionCustomBB), abs(row.lbBoscia)) == 0
             #           push!(rel_gap_custombb, row.solutionCustomBB - row.lbBoscia)
             #       elseif sign(row.lbBoscia) != sign(row.solutionCustomBB)
             #           push!(rel_gap_custombb, Inf)
             #       else
             #           push!(rel_gap_custombb, (row.solutionCustomBB - row.lbBoscia)/min(abs(row.solutionCustomBB), abs(row.lbBoscia)))
             #       end
                end

                if min(abs(row.solutionCustomBB), abs(row.lbBoscia)) == 0
                    push!(rel_gap_custombb, row.solutionCustomBB - row.lbBoscia)
                elseif sign(row.lbBoscia) != sign(row.solutionCustomBB)
                    push!(rel_gap_custombb, Inf)
                else
                    push!(rel_gap_custombb, (row.solutionCustomBB - row.lbBoscia)/min(abs(row.solutionCustomBB), abs(row.lbBoscia)))
                end

                if min(abs(row.solutionBoscia), abs(row.lbBoscia)) == 0
                    push!(rel_gap_boscia, row.solutionBoscia - row.lbBoscia)
                elseif sign(row.lbBoscia) != sign(row.solutionBoscia) || row.solutionBoscia == Inf
                    push!(rel_gap_boscia, Inf)
                else
                    push!(rel_gap_boscia, (row.solutionBoscia - row.lbBoscia)/min(abs(row.solutionBoscia), abs(row.lbBoscia)))
                end
                if min(abs(row.solutionPajarito), abs(row.lbBoscia)) == 0
                    push!(rel_gap_pajarito, row.solutionPajarito - row.lbBoscia)
                elseif sign(row.lbBoscia) != sign(row.solutionPajarito) || row.solutionPajarito == Inf
                    push!(rel_gap_pajarito, Inf)
                else
                    push!(rel_gap_pajarito, (row.solutionPajarito - row.lbBoscia)/min(abs(row.solutionPajarito), abs(row.lbBoscia)))
                end
          end
          if criterion == "AF" || criterion == "DF"
            df_dim[!, :relGapScip] = round.(rel_gap_scip,digits=3)
            relDualGapScip = vcat(relDualGapScip, rel_gap_scip)
       #   else
       #     df_dim[!,:relGapCustomBB] = round.(rel_gap_custombb, digits=3)
       #     relDualGapCustomBB = vcat(relDualGapCustomBB, rel_gap_custombb)
          end
          df_dim[!,:relGapCustomBB] = round.(rel_gap_custombb, digits=3)
            relDualGapCustomBB = vcat(relDualGapCustomBB, rel_gap_custombb)
          df_dim[!, :relGapBoscia] = round.(rel_gap_boscia, digits=3)
          relDualGapBoscia = vcat(relDualGapBoscia, rel_gap_boscia)
          df_dim[!, :relGapPajarito] = round.(rel_gap_pajarito,digits=3)
          relDualGapPajarito = vcat(relDualGapPajarito, rel_gap_pajarito)
         # df[!, :relGapHypatiaCont] = rel_gap_hypatia_cont
         # df[!, :relGapHypatiaLimit] = rel_gap_hypatia_limit
         # df[!, :relGapFrankWolfeCont] = rel_gap_frankwolfe_cont
         # df[!, :relGapFrankWolfeLimit] = rel_gap_frankwolfe_limit
  
          # save csv 
          file_name = joinpath(@__DIR__, "csv/Results/" * criterion * "_" * string(m) * "_optimality_" * type * "_non_grouped.csv")
          CSV.write(file_name, df_dim, append=false)

          println("\n")
        end

        println("Build complete non grouped csv")
        df[!,:seed] = df_condi[!,:seed]
        df[!,:numberOfExperiments] = df_condi[!,:numberOfExperiments]
        df[!,:numberOfParameters] = df_condi[!, :numberOfParameters]
        df[!,:numberOfAllowedEx] = df_condi[!,:numberOfAllowedEx]
        df[!,:frac] = df_condi[!,:frac]
        df[!,:eigmax] = df_condi[!,:eigmax]
        df[!,:eigmin] = df_condi[!,:eigmin]
        df[!,:ratio] = df_condi[!,:ratio]

        df[!,:timeBoscia] = timeBoscia
        df[!,:solutionBoscia] = solutionBoscia
        df[!,:terminationBoscia] = terminationBoscia
        df[!,:dualGapBoscia] = dualGapBoscia
        df[!,:relDualGapBoscia] = relDualGapBoscia
        df[!,:lbBoscia] = lbBoscia
        df[!,:numberNodesBoscia] = numberNodesBoscia

        if criterion in ["AF","DF"]
            df[!,:timeScip] = timeScip
            df[!,:solutionScip] = solutionScip
            df[!,:terminationScip] = terminationScip
            df[!,:dualGapScip] = dualGapScip
            df[!,:relDualGapScip] = relDualGapScip
            df[!,:numberCutsScip] = numberCutsScip
      #  else
      #      df[!,:timeCustomBB] = timeCustomBB
      #      df[!,:solutionCustomBB] = solutionCustomBB
      #      df[!,:terminationCustomBB] = terminationCustomBB
      #      df[!,:dualGapCustomBB] = dualGapCustomBB
      #      df[!,:relDualGapCustomBB] = relDualGapCustomBB
      #      df[!,:numberNodesCustomBB] = numberNodesCustomBB
        end

        df[!,:timeCustomBB] = timeCustomBB
        df[!,:solutionCustomBB] = solutionCustomBB
        df[!,:terminationCustomBB] = terminationCustomBB
        df[!,:dualGapCustomBB] = dualGapCustomBB
        df[!,:relDualGapCustomBB] = relDualGapCustomBB
        df[!,:numberNodesCustomBB] = numberNodesCustomBB

        df[!,:timePajarito] = timePajarito
        df[!,:solutionPajarito] = solutionPajarito
        df[!,:terminationPajarito] = terminationPajarito
        df[!,:dualGapPajarito] = dualGapPajarito
        df[!,:relDualGapPajarito] = relDualGapPajarito
        df[!,:numberCutsPajarito] = numberCutsPajarito
        df[!,:numberIterPajarito] = numberIterPajarito


        if criterion in ["AF", "DF"]
            df[!,:minimumTime] = min.(df[!,:timeBoscia], df[!,:timeScip], df[!,:timePajarito])
      #  else
         #   df[!,:minimumTime] = min.(df[!,:timeBoscia], df[!,:timePajarito], df[!,:timeCustomBB])
        end
        df[!,:minimumTime] = min.(df[!,:timeBoscia], df[!,:timePajarito], df[!,:timeCustomBB])

        file_name = joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_" * type * "_non_grouped.csv")
        CSV.write(file_name, df, append=false)
        println("\n")
    end
end

function build_grouped_csv(corr = true)

    function geo_mean(group)
        prod = 1.0
        n = 0
        for element in group
            # @show element
            if element != Inf
                prod = prod * element
                n += 1
            end
        end
        @show prod, n
        if n == 0
            return Inf
        end
        return prod^(1/n)
    end

    function custom_mean(group)
        sum = 0.0
        n = 0
        dash = false
        for element in group
            if element == "-"
                dash = true
                continue
            end
            if element != Inf 
                if typeof(element) == String7 || typeof(element) == String3
                    element = parse(Float64, element)
                end
                sum += element
                n += 1
            end
        end
        if n == 0
            return dash ? "-" : Inf
        end
        return sum/n
    end

   #= function geom_shifted_mean(xs; shift=big"1.0")
        n = length(xs)        
        if n != 0    
            r = prod(xi + shift for xi in xs)
            return Float64(r^(1/n) - shift)
        end
        return Inf
    end=#

    function geom_shifted_mean(xs; shift=big"1.0")
        a = length(xs)  
        n= 0
        prod = 1.0  
        dash = false
        if a != 0 
            for xi in xs
                if xi == "-"
                    dash = true
                    continue
                end
                if xi != Inf && xi != "-"
                    prod = prod*(xi+shift)  
                    n += 1
                end
            end
            return Float64(prod^(1/n) - shift)
        end
        return dash ? "-" : Inf
    end

    #criteria = corr ? ["DF", "AF", "D"] : ["A", "AF", "D", "DF"] # "A" correlated
    #criteria=["A"]
    criteria = ["A","AF","D","DF"]
    for criterion in criteria
        @show criterion
        type = corr ? "correlated" : "independent"
        df = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/" * criterion *"_optimality_" * type * "_non_grouped.csv")))

        println("Clean up Time")
        # set time below 1800
        df[df.timeBoscia.>3600, :timeBoscia] .= 3600
        #df[df.timeFrankWolfeCont.>1800, :timeFrankWolfeCont] .= 1800
        #df[df.timeFrankWolfeLimit.>1800, :timeFrankWolfeLimit] .= 1800
        if criterion == "AF" || criterion == "DF"
            df[df.timeScip.>3600, :timeScip] .= 3600
        else
            df[df.timeCustomBB.>3600, :timeCustomBB] .= 3600
        end
        #df[df.timeHypatiaCont.>1800, :timeHypatiaCont] .= 1800
        #df[df.timeHypatiaLimit.>1800, :timeHypatiaLimit] .= 1800
        df[df.timePajarito.>3600, :timePajarito] .= 3600

        println("Combine data")
        if criterion == "AF" || criterion == "DF"
        gdf = combine(
            groupby(df, [:numberOfExperiments, :numberOfParameters, :frac]),
            :eigmax => custom_mean, :eigmin => custom_mean, :ratio => custom_mean, 
            :timeBoscia => geom_shifted_mean, :terminationBoscia => sum,
            :relDualGapBoscia => custom_mean, :dualGapBoscia => custom_mean, # maybe normal mean?
            :numberNodesBoscia => custom_mean,
           # :timeFrankWolfeCont => geom_shifted_mean, :terminationFrankWolfeCont => sum,
           # :relGapFrankWolfeCont => custom_mean, :dualGapFrankWolfeCont => custom_mean,
           # :intFoundFrankWolfeCont => sum,
           # :timeFrankWolfeLimit => geom_shifted_mean, :terminationFrankWolfeLimit => sum,
           # :relGapFrankWolfeLimit => custom_mean, :dualGapFrankWolfeLimit => custom_mean,
           # :intFoundFrankWolfeLimit => sum,
            :timeScip => geom_shifted_mean, :terminationScip => sum,
            :relDualGapScip => custom_mean, :dualGapScip => custom_mean,
            :numberCutsScip => custom_mean,
            :timeCustomBB => geom_shifted_mean, :terminationCustomBB => sum,
            :relDualGapCustomBB  => custom_mean, :dualGapCustomBB => custom_mean,
            :numberNodesCustomBB => custom_mean,
            #:timeHypatiaCont => geom_shifted_mean, :terminationHypatiaCont => sum,
            #:relGapHypatiaCont => custom_mean, :dualGapHypatiaCont => custom_mean,
            #:intFoundHypatiaCont => sum,
            #:timeHypatiaLimit => geom_shifted_mean, :terminationHypatiaLimit => sum,
            #:relGapHypatiaLimit => custom_mean, :dualGapHypatiaLimit => custom_mean,
            #:intFoundHypatiaLimit => sum,
            :timePajarito => geom_shifted_mean, :terminationPajarito => sum,
            :relDualGapPajarito => custom_mean, :dualGapPajarito => custom_mean, 
            :numberCutsPajarito => custom_mean, :numberIterPajarito => custom_mean,
            nrow => :NumInstances, renamecols=false
        )
        else
            gdf = combine(
                groupby(df, [:numberOfExperiments, :numberOfParameters, :frac]),
                :eigmax => custom_mean, :eigmin => custom_mean, :ratio => custom_mean, 
                :timeBoscia => geom_shifted_mean, :terminationBoscia => sum,
                :relDualGapBoscia => custom_mean, :dualGapBoscia => custom_mean, # maybe normal mean?
                :numberNodesBoscia => custom_mean,
                #:timeFrankWolfeCont => geom_shifted_mean, :terminationFrankWolfeCont => sum,
                #:relGapFrankWolfeCont => custom_mean, :dualGapFrankWolfeCont => custom_mean,
                #:intFoundFrankWolfeCont => sum,
                #:timeFrankWolfeLimit => geom_shifted_mean, :terminationFrankWolfeLimit => sum,
                #:relGapFrankWolfeLimit => custom_mean, :dualGapFrankWolfeLimit => custom_mean,
                #:intFoundFrankWolfeLimit => sum,
                :timeCustomBB => geom_shifted_mean, :terminationCustomBB => sum,
                :relDualGapCustomBB  => custom_mean, :dualGapCustomBB => custom_mean,
                :numberNodesCustomBB => custom_mean,
                #:timeHypatiaCont => geom_shifted_mean, :terminationHypatiaCont => sum,
                #:relGapHypatiaCont => custom_mean, :dualGapHypatiaCont => custom_mean,
                #:intFoundHypatiaCont => sum,
                #:timeHypatiaLimit => geom_shifted_mean, :terminationHypatiaLimit => sum,
                #:relGapHypatiaLimit => custom_mean, :dualGapHypatiaLimit => custom_mean,
                #:intFoundHypatiaLimit => sum,
                :timePajarito => geom_shifted_mean, :terminationPajarito => sum,
                :relDualGapPajarito => custom_mean, :dualGapPajarito => custom_mean, 
                :numberCutsPajarito => custom_mean, :numberIterPajarito => custom_mean,
                nrow => :NumInstances, renamecols=false
            )
        end

        println("Clean up and rounding")
        gdf[!,:ratio] = round.(gdf[!,:ratio], digits=2)

        gdf[!,:timeBoscia] = convert.(Int64,round.(gdf[!,:timeBoscia]))
        #gdf[!,:timeFrankWolfeCont] = convert.(Int64,round.(gdf[!,:timeFrankWolfeCont]))
        #gdf[!,:timeFrankWolfeLimit] = convert.(Int64,round.(gdf[!,:timeFrankWolfeLimit]))
        if criterion == "AF" || criterion == "DF"
            gdf[!,:timeScip] = convert.(Int64,round.(gdf[!,:timeScip]))
      #  else
       #     gdf[!,:timeCustomBB] = convert.(Int64, round.(gdf[!,:timeCustomBB]))
        end
        gdf[!,:timeCustomBB] = convert.(Int64, round.(gdf[!,:timeCustomBB]))
        #gdf[!,:timeHypatiaCont] = convert.(Int64,round.(gdf[!,:timeHypatiaCont]))
        #gdf[!,:timeHypatiaLimit] = convert.(Int64,round.(gdf[!,:timeHypatiaLimit]))
        gdf[!,:timePajarito] = convert.(Int64,round.(gdf[!,:timePajarito]))

        # relative instances solved
        gdf[!,:terminationBosciaRel] = gdf[!,:terminationBoscia]./gdf[!,:NumInstances]*100
       # gdf[!,:terminationFWContRel] = gdf[!,:terminationFrankWolfeCont]./gdf[!,:NumInstances]*100
       # gdf[!,:terminationFWLimitRel] = gdf[!,:terminationFrankWolfeLimit]./gdf[!,:NumInstances]*100
        if criterion == "AF" || criterion == "DF"
            gdf[!,:terminationScipRel] = gdf[!,:terminationScip]./gdf[!,:NumInstances]*100
      #  else
       #     gdf[!,:terminationCustomBBRel] = gdf[!,:terminationCustomBB]./gdf[!,:NumInstances]*100
        end
        gdf[!,:terminationCustomBBRel] = gdf[!,:terminationCustomBB]./gdf[!,:NumInstances]*100
        #gdf[!,:terminationHypatiaContRel] = gdf[!,:terminationHypatiaCont]./gdf[!,:NumInstances]*100
        #gdf[!,:terminationHypatiaLimitRel] = gdf[!,:terminationHypatiaLimit]./gdf[!,:NumInstances]*100
        gdf[!,:terminationPajaritoRel] = gdf[!,:terminationPajarito]./gdf[!,:NumInstances]*100

        # relative integer solutions found 
        #gdf[!, :intFoundHypatiaContRel] = gdf[!,:intFoundHypatiaCont]./gdf[!,:NumInstances]*100
        #gdf[!, :intFoundHypatiaLimitRel] = gdf[!,:intFoundHypatiaLimit]./gdf[!,:NumInstances]*100
        #gdf[!, :intFoundFrankWolfeContRel] = gdf[!,:intFoundFrankWolfeCont]./gdf[!,:NumInstances]*100
        #gdf[!, :intFoundFrankWolfeLimitRel] = gdf[!,:intFoundFrankWolfeLimit]./gdf[!,:NumInstances]*100

        # parse to int
        gdf[!,:terminationBosciaRel] = convert.(Int64,round.(gdf[!,:terminationBosciaRel]))
       # gdf[!,:terminationFWContRel] = convert.(Int64,round.(gdf[!,:terminationFWContRel]))
       # gdf[!,:terminationFWLimitRel] = convert.(Int64,round.(gdf[!,:terminationFWLimitRel]))
        if criterion == "AF" || criterion == "DF"
            gdf[!,:terminationScipRel] = convert.(Int64,round.(gdf[!,:terminationScipRel]))
     #   else
      #      gdf[!,:terminationCustomBBRel] = convert.(Int64,round.(gdf[!,:terminationCustomBBRel]))
        end
        gdf[!,:terminationCustomBBRel] = convert.(Int64,round.(gdf[!,:terminationCustomBBRel]))
        #gdf[!,:terminationHypatiaContRel] = convert.(Int64,round.(gdf[!,:terminationHypatiaContRel]))
        #gdf[!,:terminationHypatiaLimitRel] = convert.(Int64,round.(gdf[!,:terminationHypatiaLimitRel]))
        gdf[!,:terminationPajaritoRel] = convert.(Int64,round.(gdf[!,:terminationPajaritoRel]))

        # Additional Data like number of nodes/custom_mean
        non_inf_entries = findall(isfinite, gdf[!,:numberNodesBoscia])
        gdf[non_inf_entries,:numberNodesBoscia] = convert.(Int64, round.(gdf[non_inf_entries, :numberNodesBoscia]))
        if criterion in ["AF","DF"]
            non_inf_entries = findall(isfinite, gdf[!,:numberCutsScip])
            gdf[non_inf_entries,:numberCutsScip] = convert.(Int64, round.(gdf[non_inf_entries, :numberCutsScip]))
       # else
        #    non_inf_entries = findall(isfinite, gdf[!,:numberNodesCustomBB])
         #   gdf[non_inf_entries,:numberNodesCustomBB] = convert.(Int64, round.(gdf[non_inf_entries, :numberNodesCustomBB]))
        end
        non_inf_entries = findall(isfinite, gdf[!,:numberNodesCustomBB])
        gdf[non_inf_entries,:numberNodesCustomBB] = convert.(Int64, round.(gdf[non_inf_entries, :numberNodesCustomBB]))
        non_inf_entries = findall(isfinite, gdf[!,:numberCutsPajarito])
        gdf[non_inf_entries,:numberCutsPajarito] = convert.(Int64, round.(gdf[non_inf_entries, :numberCutsPajarito]))
        non_inf_entries = findall(isfinite, gdf[!,:numberIterPajarito])
        gdf[non_inf_entries,:numberIterPajarito] = convert.(Int64, round.(gdf[non_inf_entries, :numberIterPajarito]))

        #gdf[!,:dualGapFrankWolfeCont] = round.(gdf[!,:dualGapFrankWolfeCont], digits=2)
        #gdf[!,:dualGapFrankWolfeLimit] = round.(gdf[!,:dualGapFrankWolfeLimit], digits=2)
        #gdf[!,:dualGapHypatiaCont] = round.(gdf[!,:dualGapHypatiaCont], digits=2)
        #gdf[!,:dualGapHypatiaLimit] = round.(gdf[!,:dualGapHypatiaLimit], digits=2)

        #gdf[!,:intFoundFrankWolfeContRel] = round.(gdf[!,:intFoundFrankWolfeContRel], digits=2)
        #gdf[!,:intFoundFrankWolfeLimitRel] = round.(gdf[!,:intFoundFrankWolfeLimitRel], digits=2)
        #gdf[!,:intFoundHypatiaContRel] = round.(gdf[!,:intFoundHypatiaContRel], digits=2)
        #gdf[!,:intFoundHypatiaLimitRel] = round.(gdf[!,:intFoundHypatiaLimitRel], digits=2)

        non_inf_entries = findall(isfinite, gdf[!, :relDualGapBoscia])
        gdf[non_inf_entries, :relDualGapBoscia] = convert.(Int64, round.(gdf[non_inf_entries, :relDualGapBoscia]*100))
        non_inf_entries = findall(isfinite, gdf[!, :relDualGapPajarito])
        gdf[non_inf_entries, :relDualGapPajarito] = convert.(Int64, round.(gdf[non_inf_entries, :relDualGapPajarito]*100))
        #non_inf_entries = findall(x-> x!= "-" && x!=Inf, gdf[!, :relGapHypatiaCont])
        #gdf[non_inf_entries, :relGapHypatiaCont] = convert.(Int64, round.(gdf[non_inf_entries, :relGapHypatiaCont]*100))
        #non_inf_entries = findall(x->x!="-" && x!=Inf, gdf[!, :relGapHypatiaLimit])
        #gdf[non_inf_entries, :relGapHypatiaLimit] = convert.(Int64, round.(gdf[non_inf_entries, :relGapHypatiaLimit]*100))
        #non_inf_entries = findall(x->x!="-" && x!=Inf, gdf[!, :relGapFrankWolfeCont])
        #gdf[non_inf_entries, :relGapFrankWolfeCont] = convert.(Int64, round.(gdf[non_inf_entries, :relGapFrankWolfeCont]*100))
        #non_inf_entries = findall(x->x!="-" && x!=Inf, gdf[!, :relGapFrankWolfeLimit])
        #gdf[non_inf_entries, :relGapFrankWolfeLimit] = convert.(Int64, round.(gdf[non_inf_entries, :relGapFrankWolfeLimit]*100))
        if criterion == "AF" || criterion == "DF"
            non_inf_entries = findall(isfinite, gdf[!, :relDualGapScip])
            gdf[non_inf_entries, :relDualGapScip] = convert.(Int64, round.(gdf[non_inf_entries, :relDualGapScip]*100))
      #  else
       #     non_inf_entries = findall(isfinite, gdf[!, :relDualGapCustomBB])
        #    gdf[non_inf_entries, :relDualGapCustomBB] = convert.(Int64, round.(gdf[non_inf_entries, :relDualGapCustomBB]*100))
        end

        non_inf_entries = findall(isfinite, gdf[!, :relDualGapCustomBB])
            gdf[non_inf_entries, :relDualGapCustomBB] = convert.(Int64, round.(gdf[non_inf_entries, :relDualGapCustomBB]*100))

        file_name = joinpath(@__DIR__, "csv/Results/" * criterion * "_optimality_" * type * "_grouped.csv")
        CSV.write(file_name, gdf, append=false)
    end
    return true
end


# #of instances #solved #avg time total #avg time solved #Nodes/Cuts #rel gap in case not solved
function build_summary_by_criterion()
    function geo_mean(group)
        prod = 1.0
        n = 0
        if isempty(group)
            return -1
        end
        for element in group
            # @show element
            if element != Inf
                prod = prod * abs(element)
                n += 1
            end
        end
        @show prod, n
        if n == 0
            return Inf
        end
        return prod^(1/n)
    end

    function geom_shifted_mean(xs; shift=big"1.0")
        a = length(xs)  
        n= 0
        prod = 1.0  
        if a != 0 
            for xi in xs
                if xi != Inf 
                    prod = prod*(xi+shift)  
                    n += 1
                end
            end
            return Float64(prod^(1/n) - shift)
        end
        return Inf
    end

    function custom_mean(group)
        sum = 0.0
        n = 0
        dash = false

        if isempty(group)
            return -1
        end
        for element in group
            if element == "-"
                dash = true
                continue
            end
            if element != Inf 
                if typeof(element) == String7 || typeof(element) == String3
                    element = parse(Float64, element)
                end
                sum += element
                n += 1
            end
        end
        if n == 0
            return dash ? "-" : Inf
        end
        return sum/n
    end

    df = DataFrame()
    NumInstances = 50

        # load data
        df_AF_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/AF_optimality_independent_non_grouped.csv")))
        df_A_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/A_optimality_independent_non_grouped.csv")))
        df_DF_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/DF_optimality_independent_non_grouped.csv")))
        df_D_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/D_optimality_independent_non_grouped.csv")))

        df_AF_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/AF_optimality_correlated_non_grouped.csv")))
        df_A_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/A_optimality_correlated_non_grouped.csv")))
        df_DF_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/DF_optimality_correlated_non_grouped.csv")))
        df_D_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/D_optimality_correlated_non_grouped.csv")))
    
        # we are only interested in entries with termination
    
        ind_term_AF_ind_Boscia = findall(x-> x==1, df_AF_ind[!,:terminationBoscia])
        ind_term_AF_ind_Scip = findall(x-> x==1, df_AF_ind[!,:terminationScip])
        ind_term_AF_ind_Pajarito = findall(x-> x==1, df_AF_ind[!,:terminationPajarito])
        ind_term_AF_ind_CustomBB = findall(x->x==1, df_AF_ind[!,:terminationCustomBB])
    
        ind_term_A_ind_Boscia = findall(x-> x==1, df_A_ind[!,:terminationBoscia])
        ind_term_A_ind_Pajarito = findall(x-> x==1, df_A_ind[!,:terminationPajarito])
        ind_term_A_ind_CustomBB = findall(x->x==1, df_A_ind[!,:terminationCustomBB])
    
        ind_term_DF_ind_Boscia = findall(x-> x==1, df_DF_ind[!,:terminationBoscia])
        ind_term_DF_ind_Scip = findall(x-> x==1, df_DF_ind[!,:terminationScip])
        ind_term_DF_ind_Pajarito = findall(x-> x==1, df_DF_ind[!,:terminationPajarito])
        ind_term_DF_ind_CustomBB = findall(x->x==1, df_DF_ind[!,:terminationCustomBB])
    
        ind_term_D_ind_Boscia = findall(x-> x==1, df_D_ind[!,:terminationBoscia])
        ind_term_D_ind_Pajarito = findall(x-> x==1, df_D_ind[!,:terminationPajarito])
        ind_term_D_ind_CustomBB = findall(x->x==1, df_D_ind[!,:terminationCustomBB])

        ind_term_AF_corr_Boscia = findall(x-> x==1, df_AF_corr[!,:terminationBoscia])
        ind_term_AF_corr_Scip = findall(x-> x==1, df_AF_corr[!,:terminationScip])
        ind_term_AF_corr_Pajarito = findall(x-> x==1, df_AF_corr[!,:terminationPajarito])
        ind_term_AF_corr_CustomBB = findall(x->x==1, df_AF_corr[!,:terminationCustomBB])
    
        ind_term_A_corr_Boscia = findall(x-> x==1, df_A_corr[!,:terminationBoscia])
        ind_term_A_corr_Pajarito = findall(x-> x==1, df_A_corr[!,:terminationPajarito])
        ind_term_A_corr_CustomBB = findall(x->x==1, df_A_corr[!,:terminationCustomBB])
    
        ind_term_DF_corr_Boscia = findall(x-> x==1, df_DF_corr[!,:terminationBoscia])
        ind_term_DF_corr_Scip = findall(x-> x==1, df_DF_corr[!,:terminationScip])
        ind_term_DF_corr_Pajarito = findall(x-> x==1, df_DF_corr[!,:terminationPajarito])
        ind_term_DF_corr_CustomBB = findall(x->x==1, df_DF_corr[!,:terminationCustomBB])
    
        ind_term_D_corr_Boscia = findall(x-> x==1, df_D_corr[!,:terminationBoscia])
        ind_term_D_corr_Pajarito = findall(x-> x==1, df_D_corr[!,:terminationPajarito])
        ind_term_D_corr_CustomBB = findall(x->x==1, df_D_corr[!,:terminationCustomBB])

        boscia_term = []
        scip_term = []
        pajarito_term = []
        custombb_term = []

    
        boscia_time = []
        scip_time = []
        pajarito_time = []
        custombb_time = []

    
        boscia_time_all = []
        scip_time_all = []
        pajarito_time_all = []
        custombb_time_all = []

    
        boscia_term_rel = []
        scip_term_rel = []
        pajarito_term_rel = []
        custombb_term_rel = []


        boscia_nodes = []
        custombb_nodes = []
        scip_cuts = []
        pajarito_cuts= []
        pajarito_iters = []


        boscia_rel_gap = []
        custombb_rel_gap = []
        scip_rel_gap = []
        pajarito_rel_gap = []

    
        ## termination
        println("Termination")
        push!(boscia_term, length(ind_term_AF_ind_Boscia))
        push!(scip_term, length(ind_term_AF_ind_Scip))
        push!(pajarito_term, length(ind_term_AF_ind_Pajarito))
      #  push!(custombb_term, 0)
        push!(custombb_term, length(ind_term_AF_ind_CustomBB))

        push!(boscia_term, length(ind_term_AF_corr_Boscia))
        push!(scip_term, length(ind_term_AF_corr_Scip))
        push!(pajarito_term, length(ind_term_AF_corr_Pajarito))
        #push!(custombb_term, 0)
        push!(custombb_term, length(ind_term_AF_corr_CustomBB))
    

        push!(boscia_term, length(ind_term_A_ind_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_term_A_corr_Pajarito))
        push!(custombb_term, length(ind_term_A_corr_CustomBB))

        push!(boscia_term, length(ind_term_A_corr_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_term_A_corr_Pajarito))
        push!(custombb_term, length(ind_term_A_corr_CustomBB))
    

        push!(boscia_term, length(ind_term_DF_ind_Boscia))
        push!(scip_term, length(ind_term_DF_ind_Scip))
        push!(pajarito_term, length(ind_term_DF_ind_Pajarito))
        #push!(custombb_term, 0)
        push!(custombb_term, length(ind_term_DF_ind_CustomBB))

        push!(boscia_term, length(ind_term_DF_corr_Boscia))
        push!(scip_term, length(ind_term_DF_corr_Scip))
        push!(pajarito_term, length(ind_term_DF_corr_Pajarito))
        #push!(custombb_term, 0)
        push!(custombb_term, length(ind_term_DF_corr_CustomBB))
    

        push!(boscia_term, length(ind_term_D_ind_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_term_D_ind_Pajarito))
        push!(custombb_term, length(ind_term_D_ind_CustomBB))

        push!(boscia_term, length(ind_term_D_corr_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_term_D_corr_Pajarito))
        push!(custombb_term, length(ind_term_D_corr_CustomBB))
    
        ## relative termination
        println("Relative termination")
        boscia_term_rel = convert.(Int64, round.(boscia_term./NumInstances*100))
        scip_term_rel = convert.(Int64, round.(scip_term./NumInstances*100))
        pajarito_term_rel = convert.(Int64, round.(pajarito_term./NumInstances*100))
        custombb_term_rel = convert.(Int64, round.(custombb_term./NumInstances*100))
    
        ## time for solved instances 
        println("Time solved instances")
        push!(boscia_time, geom_shifted_mean(df_AF_ind[ind_term_AF_ind_Boscia,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_AF_ind[ind_term_AF_ind_Scip,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_AF_ind[ind_term_AF_ind_Pajarito,:timePajarito]))
       # push!(custombb_time, 0)
        push!(custombb_time, geom_shifted_mean(df_AF_ind[ind_term_AF_ind_CustomBB,:timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_AF_corr[ind_term_AF_corr_Boscia,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_AF_corr[ind_term_AF_corr_Scip,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_AF_corr[ind_term_AF_corr_Pajarito,:timePajarito]))
      #  push!(custombb_time, 0)
      push!(custombb_time, geom_shifted_mean(df_AF_corr[ind_term_AF_corr_CustomBB,:timeCustomBB]))
            

        push!(boscia_time, geom_shifted_mean(df_A_ind[ind_term_A_ind_Boscia,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_A_ind[ind_term_A_ind_Pajarito,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_A_ind[ind_term_A_ind_CustomBB, :timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_A_corr[ind_term_A_corr_Boscia,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_A_corr[ind_term_A_corr_Pajarito,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_A_corr[ind_term_A_corr_CustomBB, :timeCustomBB]))

    
        push!(boscia_time, geom_shifted_mean(df_DF_ind[ind_term_DF_ind_Boscia,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_DF_ind[ind_term_DF_ind_Scip,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_DF_ind[ind_term_DF_ind_Pajarito,:timePajarito]))
       # push!(custombb_time, 0)
        push!(custombb_time, geom_shifted_mean(df_DF_ind[ind_term_DF_ind_CustomBB,:timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_DF_corr[ind_term_DF_corr_Boscia,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_DF_corr[ind_term_DF_corr_Scip,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_DF_corr[ind_term_DF_corr_Pajarito,:timePajarito]))
       # push!(custombb_time, 0)
       push!(custombb_time, geom_shifted_mean(df_DF_corr[ind_term_DF_corr_CustomBB,:timeCustomBB]))

    
        push!(boscia_time, geom_shifted_mean(df_D_ind[ind_term_D_ind_Boscia,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_D_ind[ind_term_D_ind_Pajarito,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_D_ind[ind_term_D_ind_CustomBB, :timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_D_corr[ind_term_D_ind_Boscia,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_D_corr[ind_term_D_ind_Pajarito,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_D_corr[ind_term_D_ind_CustomBB, :timeCustomBB]))
    
        boscia_time = round.(boscia_time, digits=3)
        pajarito_time = round.(pajarito_time, digits=3)
        scip_time = round.(scip_time, digits=3)
        custombb_time = round.(custombb_time, digits=3)

        ## number nodes and Cuts
        println("Cuts and Nodes")
        push!(boscia_nodes, custom_mean(df_AF_ind[!, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_AF_ind[!,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_AF_ind[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_AF_ind[!,:numberIterPajarito]))
       # push!(custombb_nodes, 0)
       push!(custombb_nodes, custom_mean(df_AF_ind[!,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_AF_corr[!, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_AF_corr[!,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_AF_corr[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_AF_corr[!,:numberIterPajarito]))
      #  push!(custombb_nodes, 0)
        push!(custombb_nodes, custom_mean(df_AF_corr[!,:numberNodesCustomBB]))


        push!(boscia_nodes, custom_mean(df_A_ind[!, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_A_ind[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_A_ind[!,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_A_ind[!,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_A_corr[!, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_A_corr[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_A_corr[!,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_A_corr[!,:numberNodesCustomBB]))


        push!(boscia_nodes, custom_mean(df_DF_ind[!, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_DF_ind[!,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_DF_ind[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_DF_ind[!,:numberIterPajarito]))
      #  push!(custombb_nodes, 0)
        push!(custombb_nodes, custom_mean(df_DF_ind[!,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_DF_corr[!, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_DF_corr[!,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_DF_corr[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_DF_corr[!,:numberIterPajarito]))
       # push!(custombb_nodes, 0)
        push!(custombb_nodes, custom_mean(df_DF_ind[!,:numberNodesCustomBB]))


        push!(boscia_nodes, custom_mean(df_D_ind[!, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_D_ind[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_D_ind[!,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_D_ind[!,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_D_corr[!, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_D_corr[!,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_D_corr[!,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_D_corr[!,:numberNodesCustomBB]))

        # rounding
        boscia_nodes = convert.(Int64, round.(boscia_nodes))
        scip_cuts = convert.(Int64, round.(scip_cuts))
        pajarito_cuts = convert.(Int64, round.(pajarito_cuts))
        pajarito_iters = convert.(Int64, round.(pajarito_iters))
        custombb_nodes = convert.(Int64, round.(custombb_nodes))


        ## relative gap 
        println("Relative gap")
        push!(boscia_rel_gap, geo_mean(df_AF_ind[!,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_AF_ind[!,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_AF_ind[!,:relDualGapPajarito]))
     #   push!(custombb_rel_gap, 0)
        push!(custombb_rel_gap, geo_mean(df_AF_ind[!,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_AF_corr[!,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_AF_corr[!,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_AF_corr[!,:relDualGapPajarito]))
     #   push!(custombb_rel_gap, 0)
     push!(custombb_rel_gap, geo_mean(df_AF_corr[!,:relDualGapCustomBB]))


        push!(boscia_rel_gap, geo_mean(df_A_ind[!,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_A_ind[!,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_A_ind[!,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_A_corr[!,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_A_corr[!,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_A_corr[!,:relDualGapCustomBB]))


        push!(boscia_rel_gap, geo_mean(df_DF_ind[!,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_DF_ind[!,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_DF_ind[!,:relDualGapPajarito]))
    #    push!(custombb_rel_gap, 0)
        push!(custombb_rel_gap, geo_mean(df_DF_ind[!,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_DF_corr[!,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_DF_corr[!,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_DF_corr[!,:relDualGapPajarito]))
     #   push!(custombb_rel_gap, 0)
        push!(custombb_rel_gap, geo_mean(df_DF_corr[!,:relDualGapCustomBB]))


        push!(boscia_rel_gap, geo_mean(df_D_ind[!,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_D_ind[!,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_D_ind[!,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_D_corr[!,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_D_corr[!,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_D_corr[!,:relDualGapCustomBB]))

        # rounding
        boscia_rel_gap = round.(boscia_rel_gap, digits=3)
        pajarito_rel_gap = round.(pajarito_rel_gap, digits=3)
        scip_rel_gap = round.(scip_rel_gap, digits=3)
        custombb_rel_gap = round.(custombb_rel_gap, digits=3)

        ## time for all instances 
        println("Time over all")
        push!(boscia_time_all, geom_shifted_mean(df_AF_ind[!,:timeBoscia]))
        push!(scip_time_all, geom_shifted_mean(df_AF_ind[!,:timeScip]))
        push!(pajarito_time_all, geom_shifted_mean(df_AF_ind[!,:timePajarito]))
      #  push!(custombb_time_all, 0)
        push!(custombb_time_all, geom_shifted_mean(df_AF_ind[!,:timeCustomBB]))

        push!(boscia_time_all, geom_shifted_mean(df_AF_corr[!,:timeBoscia]))
        push!(scip_time_all, geom_shifted_mean(df_AF_corr[!,:timeScip]))
        push!(pajarito_time_all, geom_shifted_mean(df_AF_corr[!,:timePajarito]))
        #push!(custombb_time_all, 0)
        push!(custombb_time_all, geom_shifted_mean(df_AF_corr[!,:timeCustomBB]))

    
        push!(scip_time_all, 0)
        push!(boscia_time_all, geom_shifted_mean(df_A_ind[!,:timeBoscia]))
        push!(pajarito_time_all, geom_shifted_mean(df_A_ind[!,:timePajarito]))
        push!(custombb_time_all, geom_shifted_mean(df_A_ind[!,:timeCustomBB]))

        push!(scip_time_all, 0)
        push!(boscia_time_all, geom_shifted_mean(df_A_corr[!,:timeBoscia]))
        push!(pajarito_time_all, geom_shifted_mean(df_A_corr[!,:timePajarito]))
        push!(custombb_time_all, geom_shifted_mean(df_A_corr[!,:timeCustomBB]))

    
        push!(boscia_time_all, geom_shifted_mean(df_DF_ind[!,:timeBoscia]))
        push!(scip_time_all, geom_shifted_mean(df_DF_ind[!,:timeScip]))
        push!(pajarito_time_all, geom_shifted_mean(df_DF_ind[!,:timePajarito]))
        #push!(custombb_time_all, 0)
        push!(custombb_time_all, geom_shifted_mean(df_DF_ind[!,:timeCustomBB]))

        push!(boscia_time_all, geom_shifted_mean(df_DF_corr[!,:timeBoscia]))
        push!(scip_time_all, geom_shifted_mean(df_DF_corr[!,:timeScip]))
        push!(pajarito_time_all, geom_shifted_mean(df_DF_corr[!,:timePajarito]))
        #push!(custombb_time_all, 0)
        push!(custombb_time_all, geom_shifted_mean(df_DF_corr[!,:timeCustomBB]))


        push!(boscia_time_all, geom_shifted_mean(df_D_ind[!,:timeBoscia]))
        push!(scip_time_all, 0)
        push!(pajarito_time_all, geom_shifted_mean(df_D_ind[!,:timePajarito]))
        push!(custombb_time_all, geom_shifted_mean(df_D_ind[!,:timeCustomBB]))

        push!(boscia_time_all, geom_shifted_mean(df_D_corr[!,:timeBoscia]))
        push!(scip_time_all, 0)
        push!(pajarito_time_all, geom_shifted_mean(df_D_corr[!,:timePajarito]))
        push!(custombb_time_all, geom_shifted_mean(df_D_corr[!,:timeCustomBB]))
    
        boscia_time_all = round.(boscia_time_all, digits=3)
        pajarito_time_all = round.(pajarito_time_all, digits=3)
        scip_time_all = round.(scip_time_all, digits=3)
        custombb_time_all = round.(custombb_time_all, digits=3)

    
        ## fill data into dataframe
        df[!,:Problem] = vcat(fill("AF",2),fill("A",2),fill("DF",2),fill("D",2))
        df[!,:Type] = repeat(["IND","CORR"],outer=4)

        df[!,:BosciaTerm] = boscia_term
        df[!,:ScipTerm] = scip_term
        df[!,:PajaritoTerm] = pajarito_term
        df[!,:CustomBBTerm] = custombb_term
    
        df[!,:BosciaTime] = boscia_time
        df[!,:ScipTime] = scip_time
        df[!,:PajaritoTime] = pajarito_time
        df[!,:CustomBBTime] = custombb_time
    
        df[!,:BosciaTermRel] = boscia_term_rel
        df[!,:ScipTermRel] = scip_term_rel
        df[!,:PajaritoTermRel] = pajarito_term_rel
        df[!,:CustomBBTermRel] = custombb_term_rel
    
        df[!,:BosciaTimeAll] = boscia_time_all
        df[!,:ScipTimeAll] = scip_time_all
        df[!,:PajaritoTimeAll] = pajarito_time_all
        df[!,:CustomBBTimeAll] = custombb_time_all

        df[!,:BosciaRelGap] = boscia_rel_gap
        df[!,:ScipRelGap] = scip_rel_gap
        df[!,:PajaritoRelGap] = pajarito_rel_gap
        df[!,:CustomBBRelGap] = custombb_rel_gap

        df[!,:BosciaNodes] = boscia_nodes
        df[!,:ScipCuts] = scip_cuts
        df[!,:PajaritoCuts] = pajarito_cuts
        df[!,:PajaritoIter] = pajarito_iters
        df[!,:CustomBBNodes] = custombb_nodes
    
    
        file_name = joinpath(@__DIR__, "csv/summary_by_criterion.csv")
        CSV.write(file_name, df, append=false)
        return true
end

# #of instances #solved #avg time total #avg time solved #Nodes/Cuts
function build_summary_by_dimension()
    function geo_mean(group)
        prod = 1.0
        n = 0
        if isempty(group)
            return -1
        end
        for element in group
            # @show element
            if element != Inf
                prod = prod * abs(element)
                n += 1
            end
        end
        @show prod, n
        if n == 0
            return Inf
        end
        return prod^(1/n)
    end

    function geom_shifted_mean(xs; shift=big"1.0")
        a = length(xs)  
        n= 0
        prod = 1.0  
        if a != 0 
            for xi in xs
                if xi != Inf 
                    prod = prod*(xi+shift)  
                    n += 1
                end
            end
            return Float64(prod^(1/n) - shift)
        end
        return Inf
    end

    function custom_mean(group)
        sum = 0.0
        n = 0
        dash = false

        if isempty(group)
            return -1
        end
        for element in group
            if element == "-"
                dash = true
                continue
            end
            if element != Inf 
                if typeof(element) == String7 || typeof(element) == String3
                    element = parse(Float64, element)
                end
                sum += element
                n += 1
            end
        end
        if n == 0
            return dash ? "-" : Inf
        end
        return sum/n
    end

    df = DataFrame()
    NumInstances = 40

    # storage
    boscia_term = []
    scip_term = []
    pajarito_term = []
    custombb_term = []

    boscia_time = []
    scip_time = []
    pajarito_time = []
    custombb_time = []

    boscia_time_all = []
    scip_time_all = []
    pajarito_time_all = []
    custombb_time_all = []
    
    boscia_term_rel = []
    scip_term_rel = []
    pajarito_term_rel = []
    custombb_term_rel = []

    boscia_nodes = []
    custombb_nodes = []
    scip_cuts = []
    pajarito_cuts= []
    pajarito_iters = []

    boscia_rel_gap = []
    custombb_rel_gap = []
    scip_rel_gap = []
    pajarito_rel_gap = []

    for m in [50,60,80,100,120]
        @show m

        # load data
        df_AF_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/AF_" * string(m) * "_optimality_independent_non_grouped.csv")))
        df_A_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/A_" * string(m) * "_optimality_independent_non_grouped.csv")))
        df_DF_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/DF_" * string(m) * "_optimality_independent_non_grouped.csv")))
        df_D_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/D_" * string(m) * "_optimality_independent_non_grouped.csv")))

        df_AF_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/AF_" * string(m) * "_optimality_correlated_non_grouped.csv")))
        df_A_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/A_" * string(m) * "_optimality_correlated_non_grouped.csv")))
        df_DF_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/DF_" * string(m) * "_optimality_correlated_non_grouped.csv")))
        df_D_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/D_" * string(m) * "_optimality_correlated_non_grouped.csv")))
    
        # we are only interested in entries with termination
    
        ind_term_AF_ind_Boscia = findall(x-> x==1, df_AF_ind[!,:terminationBoscia])
        ind_term_AF_ind_Scip = findall(x-> x==1, df_AF_ind[!,:terminationScip])
        ind_term_AF_ind_Pajarito = findall(x-> x==1, df_AF_ind[!,:terminationPajarito])
    
        ind_term_A_ind_Boscia = findall(x-> x==1, df_A_ind[!,:terminationBoscia])
        ind_term_A_ind_Pajarito = findall(x-> x==1, df_A_ind[!,:terminationPajarito])
        ind_term_A_ind_CustomBB = findall(x->x==1, df_A_ind[!,:terminationCustomBB])
    
        ind_term_DF_ind_Boscia = findall(x-> x==1, df_DF_ind[!,:terminationBoscia])
        ind_term_DF_ind_Scip = findall(x-> x==1, df_DF_ind[!,:terminationScip])
        ind_term_DF_ind_Pajarito = findall(x-> x==1, df_DF_ind[!,:terminationPajarito])
    
        ind_term_D_ind_Boscia = findall(x-> x==1, df_D_ind[!,:terminationBoscia])
        ind_term_D_ind_Pajarito = findall(x-> x==1, df_D_ind[!,:terminationPajarito])
        ind_term_D_ind_CustomBB = findall(x->x==1, df_D_ind[!,:terminationCustomBB])

        ind_term_AF_corr_Boscia = findall(x-> x==1, df_AF_corr[!,:terminationBoscia])
        ind_term_AF_corr_Scip = findall(x-> x==1, df_AF_corr[!,:terminationScip])
        ind_term_AF_corr_Pajarito = findall(x-> x==1, df_AF_corr[!,:terminationPajarito])
    
        ind_term_A_corr_Boscia = findall(x-> x==1, df_A_corr[!,:terminationBoscia])
        ind_term_A_corr_Pajarito = findall(x-> x==1, df_A_corr[!,:terminationPajarito])
        ind_term_A_corr_CustomBB = findall(x->x==1, df_A_corr[!,:terminationCustomBB])
    
        ind_term_DF_corr_Boscia = findall(x-> x==1, df_DF_corr[!,:terminationBoscia])
        ind_term_DF_corr_Scip = findall(x-> x==1, df_DF_corr[!,:terminationScip])
        ind_term_DF_corr_Pajarito = findall(x-> x==1, df_DF_corr[!,:terminationPajarito])
    
        ind_term_D_corr_Boscia = findall(x-> x==1, df_D_corr[!,:terminationBoscia])
        ind_term_D_corr_Pajarito = findall(x-> x==1, df_D_corr[!,:terminationPajarito])
        ind_term_D_corr_CustomBB = findall(x->x==1, df_D_corr[!,:terminationCustomBB])

        
        ## termination
        println("Termination")
        push!(boscia_term, length(ind_term_AF_ind_Boscia) + length(ind_term_A_ind_Boscia) + length(ind_term_DF_ind_Boscia) +length(ind_term_D_ind_Boscia))
        push!(scip_term, length(ind_term_AF_ind_Scip) + length(ind_term_DF_ind_Scip))
        push!(pajarito_term, length(ind_term_AF_ind_Pajarito) + length(ind_term_A_ind_Pajarito) + length(ind_term_DF_ind_Pajarito) +length(ind_term_D_ind_Pajarito))
        push!(custombb_term, length(ind_term_A_ind_CustomBB) + length(ind_term_D_ind_CustomBB))

        push!(boscia_term, length(ind_term_AF_corr_Boscia) + length(ind_term_A_corr_Boscia) + length(ind_term_DF_corr_Boscia) +length(ind_term_D_corr_Boscia))
        push!(scip_term, length(ind_term_AF_corr_Scip) + length(ind_term_DF_corr_Scip))
        push!(pajarito_term, length(ind_term_AF_corr_Pajarito) + length(ind_term_A_corr_Pajarito) + length(ind_term_DF_corr_Pajarito) +length(ind_term_D_corr_Pajarito))
        push!(custombb_term, length(ind_term_A_corr_CustomBB) + length(ind_term_D_corr_CustomBB))
    

        ## time for solved instances 
        println("Time solved instances")
        timeB = vcat(df_AF_ind[ind_term_AF_ind_Boscia,:timeBoscia], df_A_ind[ind_term_A_ind_Boscia,:timeBoscia],df_DF_ind[ind_term_DF_ind_Boscia,:timeBoscia],df_D_ind[ind_term_D_ind_Boscia,:timeBoscia])
        timeP = vcat(df_AF_ind[ind_term_AF_ind_Pajarito,:timePajarito], df_A_ind[ind_term_A_ind_Pajarito,:timePajarito], df_DF_ind[ind_term_DF_ind_Pajarito,:timePajarito], df_D_ind[ind_term_D_ind_Pajarito,:timePajarito])
        timeC = vcat(df_A_ind[ind_term_A_ind_CustomBB, :timeCustomBB], df_D_ind[ind_term_D_ind_CustomBB, :timeCustomBB])
        timeS = vcat(df_AF_ind[ind_term_AF_ind_Scip,:timeScip],df_DF_ind[ind_term_DF_ind_Scip,:timeScip])

        push!(boscia_time, geom_shifted_mean(timeB))
        push!(scip_time, geom_shifted_mean(timeS))
        push!(pajarito_time, geom_shifted_mean(timeP))
        push!(custombb_time, geom_shifted_mean(timeC))


        timeB = vcat(df_AF_corr[ind_term_AF_corr_Boscia,:timeBoscia], df_A_corr[ind_term_A_corr_Boscia,:timeBoscia],df_DF_corr[ind_term_DF_corr_Boscia,:timeBoscia],df_D_corr[ind_term_D_corr_Boscia,:timeBoscia])
        timeP = vcat(df_AF_corr[ind_term_AF_corr_Pajarito,:timePajarito], df_A_corr[ind_term_A_corr_Pajarito,:timePajarito], df_DF_corr[ind_term_DF_corr_Pajarito,:timePajarito], df_D_corr[ind_term_D_corr_Pajarito,:timePajarito])
        timeC = vcat(df_A_corr[ind_term_A_corr_CustomBB, :timeCustomBB], df_D_corr[ind_term_D_corr_CustomBB, :timeCustomBB])
        timeS = vcat(df_AF_corr[ind_term_AF_corr_Scip,:timeScip],df_DF_corr[ind_term_DF_corr_Scip,:timeScip])

        push!(boscia_time, geom_shifted_mean(timeB))
        push!(scip_time, geom_shifted_mean(timeS))
        push!(pajarito_time, geom_shifted_mean(timeP))
        push!(custombb_time, geom_shifted_mean(timeC))
    

        ## number nodes and Cuts
        println("Cuts and Nodes")
        nodesB = vcat(df_AF_ind[!, :numberNodesBoscia], df_A_ind[!, :numberNodesBoscia], df_DF_ind[!, :numberNodesBoscia], df_D_ind[!, :numberNodesBoscia])
        cutsS = vcat(df_AF_ind[!,:numberCutsScip],df_DF_ind[!,:numberCutsScip])
        iterP = vcat(df_AF_ind[!,:numberIterPajarito], df_A_ind[!,:numberIterPajarito], df_DF_ind[!,:numberIterPajarito],df_D_ind[!,:numberIterPajarito])
        cutsP = vcat(df_AF_ind[!,:numberCutsPajarito],df_A_ind[!,:numberCutsPajarito],df_DF_ind[!,:numberCutsPajarito],df_D_ind[!,:numberCutsPajarito])
        nodesC = vcat(df_A_ind[!,:numberNodesCustomBB], df_D_ind[!,:numberNodesCustomBB])

        push!(boscia_nodes, custom_mean(nodesB))
        push!(scip_cuts, custom_mean(cutsS))
        push!(pajarito_cuts, custom_mean(cutsP))
        push!(pajarito_iters, custom_mean(iterP))
        push!(custombb_nodes, custom_mean(nodesC))


        nodesB = vcat(df_AF_corr[!, :numberNodesBoscia], df_A_corr[!, :numberNodesBoscia], df_DF_corr[!, :numberNodesBoscia], df_D_corr[!, :numberNodesBoscia])
        cutsS = vcat(df_AF_corr[!,:numberCutsScip],df_DF_corr[!,:numberCutsScip])
        iterP = vcat(df_AF_corr[!,:numberIterPajarito], df_A_corr[!,:numberIterPajarito], df_DF_corr[!,:numberIterPajarito],df_D_corr[!,:numberIterPajarito])
        cutsP = vcat(df_AF_corr[!,:numberCutsPajarito],df_A_corr[!,:numberCutsPajarito],df_DF_corr[!,:numberCutsPajarito],df_D_corr[!,:numberCutsPajarito])
        nodesC = vcat(df_A_corr[!,:numberNodesCustomBB], df_D_corr[!,:numberNodesCustomBB])

        push!(boscia_nodes, custom_mean(nodesB))
        push!(scip_cuts, custom_mean(cutsS))
        push!(pajarito_cuts, custom_mean(cutsP))
        push!(pajarito_iters, custom_mean(iterP))
        push!(custombb_nodes, custom_mean(nodesC))

        ## relative gap 
        println("Relative gap")
        relGapB = vcat(df_AF_ind[!,:relGapBoscia],df_A_ind[!,:relGapBoscia],df_DF_ind[!,:relGapBoscia],df_D_ind[!,:relGapBoscia])
        relGapP = vcat(df_AF_ind[!,:relGapPajarito],df_A_ind[!,:relGapPajarito],df_DF_ind[!,:relGapPajarito],df_D_ind[!,:relGapPajarito])
        relGapS = vcat(df_AF_ind[!,:relGapScip],df_DF_ind[!,:relGapScip])
        relGapC = vcat(df_A_ind[!,:relGapCustomBB],df_D_ind[!,:relGapCustomBB])

        push!(boscia_rel_gap, geo_mean(relGapB))
        push!(scip_rel_gap, geo_mean(relGapS))
        push!(pajarito_rel_gap, geo_mean(relGapP))
        push!(custombb_rel_gap, geo_mean(relGapC))


        relGapB = vcat(df_AF_corr[!,:relGapBoscia],df_A_corr[!,:relGapBoscia],df_DF_corr[!,:relGapBoscia],df_D_corr[!,:relGapBoscia])
        relGapP = vcat(df_AF_corr[!,:relGapPajarito],df_A_corr[!,:relGapPajarito],df_DF_corr[!,:relGapPajarito],df_D_corr[!,:relGapPajarito])
        relGapS = vcat(df_AF_corr[!,:relGapScip],df_DF_corr[!,:relGapScip])
        relGapC = vcat(df_A_corr[!,:relGapCustomBB],df_D_corr[!,:relGapCustomBB])

        push!(boscia_rel_gap, geo_mean(relGapB))
        push!(scip_rel_gap, geo_mean(relGapS))
        push!(pajarito_rel_gap, geo_mean(relGapP))
        push!(custombb_rel_gap, geo_mean(relGapC))

        ## time for all instances 
        println("Time over all")
        timeB = vcat(df_AF_ind[!,:timeBoscia], df_A_ind[!,:timeBoscia],df_DF_ind[!,:timeBoscia],df_D_ind[!,:timeBoscia])
        timeP = vcat(df_AF_ind[!,:timePajarito], df_A_ind[!,:timePajarito], df_DF_ind[!,:timePajarito], df_D_ind[!,:timePajarito])
        timeC = vcat(df_A_ind[!, :timeCustomBB], df_D_ind[!, :timeCustomBB])
        timeS = vcat(df_AF_ind[!,:timeScip],df_DF_ind[!,:timeScip])

        push!(boscia_time_all, geom_shifted_mean(timeB))
        push!(scip_time_all, geom_shifted_mean(timeS))
        push!(pajarito_time_all, geom_shifted_mean(timeP))
        push!(custombb_time_all, geom_shifted_mean(timeC))


        timeB = vcat(df_AF_corr[!,:timeBoscia], df_A_corr[!,:timeBoscia],df_DF_corr[!,:timeBoscia],df_D_corr[!,:timeBoscia])
        timeP = vcat(df_AF_corr[!,:timePajarito], df_A_corr[!,:timePajarito], df_DF_corr[!,:timePajarito], df_D_corr[!,:timePajarito])
        timeC = vcat(df_A_corr[!, :timeCustomBB], df_D_corr[!, :timeCustomBB])
        timeS = vcat(df_AF_corr[!,:timeScip],df_DF_corr[!,:timeScip])

        push!(boscia_time_all, geom_shifted_mean(timeB))
        push!(scip_time_all, geom_shifted_mean(timeS))
        push!(pajarito_time_all, geom_shifted_mean(timeP))
        push!(custombb_time_all, geom_shifted_mean(timeC))

    end

    ## relative termination
    println("Relative termination")
    boscia_term_rel = convert.(Int64, round.(boscia_term./NumInstances*100))
    scip_term_rel = convert.(Int64, round.(scip_term./NumInstances*100))
    pajarito_term_rel = convert.(Int64, round.(pajarito_term./NumInstances*100))
    custombb_term_rel = convert.(Int64, round.(custombb_term./NumInstances*100))

    # rounding time
    boscia_time = round.(boscia_time, digits=3)
    pajarito_time = round.(pajarito_time, digits=3)
    scip_time = round.(scip_time, digits=3)
    custombb_time = round.(custombb_time, digits=3)

    # rounding nodes/cuts
    boscia_nodes = convert.(Int64, round.(boscia_nodes))
    scip_cuts = convert.(Int64, round.(scip_cuts))
    pajarito_cuts = convert.(Int64, round.(pajarito_cuts))
    pajarito_iters = convert.(Int64, round.(pajarito_iters))
    custombb_nodes = convert.(Int64, round.(custombb_nodes))
    
    # rounding total time
    boscia_time_all = round.(boscia_time_all, digits=3)
    pajarito_time_all = round.(pajarito_time_all, digits=3)
    scip_time_all = round.(scip_time_all, digits=3)
    custombb_time_all = round.(custombb_time_all, digits=3)

    # rounding rel gap
    boscia_rel_gap = round.(boscia_rel_gap, digits=3)
    pajarito_rel_gap = round.(pajarito_rel_gap, digits=3)
    scip_rel_gap = round.(scip_rel_gap, digits=3)
    custombb_rel_gap = round.(custombb_rel_gap, digits=3)

    
    ## fill data into dataframe
    df[!,:Dimension] = vcat(fill(50,2),fill(60,2),fill(80,2),fill(100,2), fill(120,2))
    df[!,:Type] = repeat(["IND","CORR"],outer=5)

    df[!,:BosciaTerm] = boscia_term
    df[!,:ScipTerm] = scip_term
    df[!,:PajaritoTerm] = pajarito_term
    df[!,:CustomBBTerm] = custombb_term
    
    df[!,:BosciaTime] = boscia_time
    df[!,:ScipTime] = scip_time
    df[!,:PajaritoTime] = pajarito_time
    df[!,:CustomBBTime] = custombb_time
    
    df[!,:BosciaTermRel] = boscia_term_rel
    df[!,:ScipTermRel] = scip_term_rel
    df[!,:PajaritoTermRel] = pajarito_term_rel
    df[!,:CustomBBTermRel] = custombb_term_rel
    
    df[!,:BosciaTimeAll] = boscia_time_all
    df[!,:ScipTimeAll] = scip_time_all
    df[!,:PajaritoTimeAll] = pajarito_time_all
    df[!,:CustomBBTimeAll] = custombb_time_all

    df[!,:BosciaRelGap] = boscia_rel_gap
    df[!,:ScipRelGap] = scip_rel_gap
    df[!,:PajaritoRelGap] = pajarito_rel_gap
    df[!,:CustomBBRelGap] = custombb_rel_gap

    df[!,:BosciaNodes] = boscia_nodes
    df[!,:ScipCuts] = scip_cuts
    df[!,:PajaritoCuts] = pajarito_cuts
    df[!,:PajaritoIter] = pajarito_iters
    df[!,:CustomBBNodes] = custombb_nodes
    
    file_name = joinpath(@__DIR__, "csv/summary_by_dimension.csv")
    CSV.write(file_name, df, append=false)
    return true
end

# #of instances #solved #avg time total #avg time solved #Nodes/Cuts #rel gap in case not solved
function build_summary_by_difficulty()
    function geo_mean(group)
        prod = 1.0
        n = 0
        if isempty(group)
            return -1
        end
        for element in group
            # @show element
            if element != Inf
                prod = prod * abs(element)
                n += 1
            end
        end
        if n == 0
            return Inf
        end
        return prod^(1/n)
    end

    function geom_shifted_mean(xs; shift=big"1.0")
        a = length(xs)  
        n= 0
        prod = 1.0  
        if a != 0 
            for xi in xs
                if xi != Inf 
                    prod = prod*(xi+shift)  
                    n += 1
                end
            end
            return Float64(prod^(1/n) - shift)
        end
        return Inf
    end

    function custom_mean(group)
        sum = 0.0
        n = 0
        dash = false

        if isempty(group)
            return -1
        end
        for element in group
            if element == "-"
                dash = true
                continue
            end
            if element != Inf 
                if typeof(element) == String7 || typeof(element) == String3
                    element = parse(Float64, element)
                end
                sum += element
                n += 1
            end
        end
        if n == 0
            return dash ? "-" : Inf
        end
        return sum/n
    end

    df = DataFrame()
    NumInstances = 50

    # load data
    df_AF_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/AF_optimality_independent_non_grouped.csv")))
    df_A_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/A_optimality_independent_non_grouped.csv")))
    df_DF_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/DF_optimality_independent_non_grouped.csv")))
    df_D_ind = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/D_optimality_independent_non_grouped.csv")))

    df_AF_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/AF_optimality_correlated_non_grouped.csv")))
    df_A_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/A_optimality_correlated_non_grouped.csv")))
    df_DF_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/DF_optimality_correlated_non_grouped.csv")))
    df_D_corr = DataFrame(CSV.File(joinpath(@__DIR__, "csv/Results/D_optimality_correlated_non_grouped.csv")))
    
    # we are only interested in entries with termination
    
    ind_term_AF_ind_Boscia = findall(x-> x==1, df_AF_ind[!,:terminationBoscia])
    ind_term_AF_ind_Scip = findall(x-> x==1, df_AF_ind[!,:terminationScip])
    ind_term_AF_ind_Pajarito = findall(x-> x==1, df_AF_ind[!,:terminationPajarito])
    ind_term_AF_ind_CustomBB = findall(x->x==1, df_AF_ind[!,:terminationCustomBB])
    
    ind_term_A_ind_Boscia = findall(x-> x==1, df_A_ind[!,:terminationBoscia])
    ind_term_A_ind_Pajarito = findall(x-> x==1, df_A_ind[!,:terminationPajarito])
    ind_term_A_ind_CustomBB = findall(x->x==1, df_A_ind[!,:terminationCustomBB])
    
    ind_term_DF_ind_Boscia = findall(x-> x==1, df_DF_ind[!,:terminationBoscia])
    ind_term_DF_ind_Scip = findall(x-> x==1, df_DF_ind[!,:terminationScip])
    ind_term_DF_ind_Pajarito = findall(x-> x==1, df_DF_ind[!,:terminationPajarito])
    ind_term_DF_ind_CustomBB = findall(x->x==1, df_DF_ind[!,:terminationCustomBB])
    
    ind_term_D_ind_Boscia = findall(x-> x==1, df_D_ind[!,:terminationBoscia])
    ind_term_D_ind_Pajarito = findall(x-> x==1, df_D_ind[!,:terminationPajarito])
    ind_term_D_ind_CustomBB = findall(x->x==1, df_D_ind[!,:terminationCustomBB])

    ind_term_AF_corr_Boscia = findall(x-> x==1, df_AF_corr[!,:terminationBoscia])
    ind_term_AF_corr_Scip = findall(x-> x==1, df_AF_corr[!,:terminationScip])
    ind_term_AF_corr_Pajarito = findall(x-> x==1, df_AF_corr[!,:terminationPajarito])
    ind_term_AF_corr_CustomBB = findall(x->x==1, df_AF_corr[!,:terminationCustomBB])
    
    ind_term_A_corr_Boscia = findall(x-> x==1, df_A_corr[!,:terminationBoscia])
    ind_term_A_corr_Pajarito = findall(x-> x==1, df_A_corr[!,:terminationPajarito])
    ind_term_A_corr_CustomBB = findall(x->x==1, df_A_corr[!,:terminationCustomBB])
    
    ind_term_DF_corr_Boscia = findall(x-> x==1, df_DF_corr[!,:terminationBoscia])
    ind_term_DF_corr_Scip = findall(x-> x==1, df_DF_corr[!,:terminationScip])
    ind_term_DF_corr_Pajarito = findall(x-> x==1, df_DF_corr[!,:terminationPajarito])
    ind_term_DF_corr_CustomBB = findall(x->x==1, df_DF_corr[!,:terminationCustomBB])
    
    ind_term_D_corr_Boscia = findall(x-> x==1, df_D_corr[!,:terminationBoscia])
    ind_term_D_corr_Pajarito = findall(x-> x==1, df_D_corr[!,:terminationPajarito])
    ind_term_D_corr_CustomBB = findall(x->x==1, df_D_corr[!,:terminationCustomBB])

    boscia_term = []
    scip_term = []
    pajarito_term = []
    custombb_term = []

    
    boscia_time = []
    scip_time = []
    pajarito_time = []
    custombb_time = []

    
    boscia_term_rel = []
    scip_term_rel = []
    pajarito_term_rel = []
    custombb_term_rel = []


    boscia_nodes = []
    custombb_nodes = []
    scip_cuts = []
    pajarito_cuts= []
    pajarito_iters = []


    boscia_rel_gap = []
    custombb_rel_gap = []
    scip_rel_gap = []
    pajarito_rel_gap = []

    num_instances = []

    for time in [0,10,100,1000,2000]
        println("\n")
        @show time

        ## find the time indices
        ind_min_time_AF_ind = findall(x-> x > time, df_AF_ind[!,:minimumTime])
        ind_min_time_A_ind = findall(x-> x > time, df_A_ind[!,:minimumTime])
        ind_min_time_DF_ind = findall(x-> x > time, df_DF_ind[!,:minimumTime])
        ind_min_time_D_ind = findall(x-> x > time, df_D_ind[!,:minimumTime])

        ind_min_time_AF_corr = findall(x-> x > time, df_AF_corr[!,:minimumTime])
        ind_min_time_A_corr = findall(x-> x > time, df_A_corr[!,:minimumTime])
        ind_min_time_DF_corr = findall(x-> x > time, df_DF_corr[!,:minimumTime])
        ind_min_time_D_corr = findall(x-> x > time, df_D_corr[!,:minimumTime])

        push!(num_instances, length(ind_min_time_AF_ind))
        push!(num_instances, length(ind_min_time_AF_corr))

        push!(num_instances, length(ind_min_time_A_ind))
        push!(num_instances, length(ind_min_time_A_corr))

        push!(num_instances, length(ind_min_time_DF_ind))
        push!(num_instances, length(ind_min_time_DF_corr))

        push!(num_instances, length(ind_min_time_D_ind))
        push!(num_instances, length(ind_min_time_D_corr))

        ## Intersection between solved instance within the time limit 
        println("Intersection solved minimum time and termination")
        ind_intersect_AF_ind_Boscia = intersect(ind_min_time_AF_ind, ind_term_AF_ind_Boscia)
        ind_intersect_AF_ind_Scip = intersect(ind_min_time_AF_ind, ind_term_AF_ind_Scip)
        ind_intersect_AF_ind_Pajarito = intersect(ind_min_time_AF_ind, ind_term_AF_ind_Pajarito)
        ind_intersect_AF_ind_CustomBB = intersect(ind_min_time_AF_ind, ind_term_AF_ind_CustomBB)

        ind_intersect_A_ind_Boscia = intersect(ind_min_time_A_ind, ind_term_A_ind_Boscia)
        ind_intersect_A_ind_CustomBB = intersect(ind_min_time_A_ind, ind_term_A_ind_CustomBB)
        ind_intersect_A_ind_Pajarito = intersect(ind_min_time_A_ind, ind_term_A_ind_Pajarito)

        ind_intersect_DF_ind_Boscia = intersect(ind_min_time_DF_ind, ind_term_DF_ind_Boscia)
        ind_intersect_DF_ind_Scip = intersect(ind_min_time_DF_ind, ind_term_DF_ind_Scip)
        ind_intersect_DF_ind_Pajarito = intersect(ind_min_time_DF_ind, ind_term_DF_ind_Pajarito)
        ind_intersect_DF_ind_CustomBB = intersect(ind_min_time_DF_ind, ind_term_DF_ind_CustomBB)

        ind_intersect_D_ind_Boscia = intersect(ind_min_time_D_ind, ind_term_D_ind_Boscia)
        ind_intersect_D_ind_CustomBB = intersect(ind_min_time_D_ind, ind_term_D_ind_CustomBB)
        ind_intersect_D_ind_Pajarito = intersect(ind_min_time_D_ind, ind_term_D_ind_Pajarito)

        ind_intersect_AF_corr_Boscia = intersect(ind_min_time_AF_corr, ind_term_AF_corr_Boscia)
        ind_intersect_AF_corr_Scip = intersect(ind_min_time_AF_corr, ind_term_AF_corr_Scip)
        ind_intersect_AF_corr_Pajarito = intersect(ind_min_time_AF_corr, ind_term_AF_corr_Pajarito)
        ind_intersect_AF_corr_CustomBB = intersect(ind_min_time_AF_corr, ind_term_AF_corr_CustomBB)

        ind_intersect_A_corr_Boscia = intersect(ind_min_time_A_corr, ind_term_A_corr_Boscia)
        ind_intersect_A_corr_CustomBB = intersect(ind_min_time_A_corr, ind_term_A_corr_CustomBB)
        ind_intersect_A_corr_Pajarito = intersect(ind_min_time_A_corr, ind_term_A_corr_Pajarito)

        ind_intersect_DF_corr_Boscia = intersect(ind_min_time_DF_corr, ind_term_DF_corr_Boscia)
        ind_intersect_DF_corr_Scip = intersect(ind_min_time_DF_corr, ind_term_DF_corr_Scip)
        ind_intersect_DF_corr_Pajarito = intersect(ind_min_time_DF_corr, ind_term_DF_corr_Pajarito)
        ind_intersect_DF_corr_CustomBB = intersect(ind_min_time_DF_corr, ind_term_DF_corr_CustomBB)

        ind_intersect_D_corr_Boscia = intersect(ind_min_time_D_corr, ind_term_D_corr_Boscia)
        ind_intersect_D_corr_CustomBB = intersect(ind_min_time_D_corr, ind_term_D_corr_CustomBB)
        ind_intersect_D_corr_Pajarito = intersect(ind_min_time_D_corr, ind_term_D_corr_Pajarito)

        ## termination
        println("Termination")
        push!(boscia_term, length(ind_intersect_AF_ind_Boscia))
        push!(scip_term, length(ind_intersect_AF_ind_Scip))
        push!(pajarito_term, length(ind_intersect_AF_ind_Pajarito))
       # push!(custombb_term, 0)
        push!(custombb_term, length(ind_intersect_AF_ind_CustomBB))

        push!(boscia_term, length(ind_intersect_AF_corr_Boscia))
        push!(scip_term, length(ind_intersect_AF_corr_Scip))
        push!(pajarito_term, length(ind_intersect_AF_corr_Pajarito))
       # push!(custombb_term, 0)
        push!(custombb_term, length(ind_intersect_AF_corr_CustomBB))
    

        push!(boscia_term, length(ind_intersect_A_ind_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_intersect_A_ind_Pajarito))
        push!(custombb_term, length(ind_intersect_A_ind_CustomBB))

        push!(boscia_term, length(ind_intersect_A_corr_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_intersect_A_corr_Pajarito))
        push!(custombb_term, length(ind_intersect_A_corr_CustomBB))
    

        push!(boscia_term, length(ind_intersect_DF_ind_Boscia))
        push!(scip_term, length(ind_intersect_DF_ind_Scip))
        push!(pajarito_term, length(ind_intersect_DF_ind_Pajarito))
       # push!(custombb_term, 0)
        push!(custombb_term, length(ind_intersect_DF_ind_CustomBB))

        push!(boscia_term, length(ind_intersect_DF_corr_Boscia))
        push!(scip_term, length(ind_intersect_DF_corr_Scip))
        push!(pajarito_term, length(ind_intersect_DF_corr_Pajarito))
       # push!(custombb_term, 0)
        push!(custombb_term, length(ind_intersect_DF_corr_CustomBB))
    

        push!(boscia_term, length(ind_intersect_D_ind_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_intersect_D_ind_Pajarito))
        push!(custombb_term, length(ind_intersect_D_ind_CustomBB))

        push!(boscia_term, length(ind_intersect_D_corr_Boscia))
        push!(scip_term, 0)
        push!(pajarito_term, length(ind_intersect_D_corr_Pajarito))
        push!(custombb_term, length(ind_intersect_D_corr_CustomBB))
    
        ## time for solved instances 
        println("Time solved instances")
        push!(boscia_time, geom_shifted_mean(df_AF_ind[ind_min_time_AF_ind,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_AF_ind[ind_min_time_AF_ind,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_AF_ind[ind_min_time_AF_ind,:timePajarito]))
     #   push!(custombb_time, 0)
        push!(custombb_time, geom_shifted_mean(df_AF_ind[ind_min_time_AF_ind,:timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_AF_corr[ind_min_time_AF_corr,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_AF_corr[ind_min_time_AF_corr,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_AF_corr[ind_min_time_AF_corr,:timePajarito]))
       # push!(custombb_time, 0)
        push!(custombb_time, geom_shifted_mean(df_AF_corr[ind_min_time_AF_corr,:timeCustomBB]))
            

        push!(boscia_time, geom_shifted_mean(df_A_ind[ind_min_time_A_ind,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_A_ind[ind_min_time_A_ind,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_A_ind[ind_min_time_A_ind, :timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_A_corr[ind_min_time_A_corr,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_A_corr[ind_min_time_A_corr,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_A_corr[ind_min_time_A_corr, :timeCustomBB]))

    
        push!(boscia_time, geom_shifted_mean(df_DF_ind[ind_min_time_DF_ind,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_DF_ind[ind_min_time_DF_ind,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_DF_ind[ind_min_time_DF_ind,:timePajarito]))
     #   push!(custombb_time, 0)
        push!(custombb_time, geom_shifted_mean(df_DF_ind[ind_min_time_DF_ind,:timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_DF_corr[ind_min_time_DF_corr,:timeBoscia]))
        push!(scip_time, geom_shifted_mean(df_DF_corr[ind_min_time_DF_corr,:timeScip]))
        push!(pajarito_time, geom_shifted_mean(df_DF_corr[ind_min_time_DF_corr,:timePajarito]))
       # push!(custombb_time, 0)
        push!(custombb_time, geom_shifted_mean(df_DF_corr[ind_min_time_DF_corr,:timeCustomBB]))

    
        push!(boscia_time, geom_shifted_mean(df_D_ind[ind_min_time_D_ind,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_D_ind[ind_min_time_D_ind,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_D_ind[ind_min_time_D_ind, :timeCustomBB]))

        push!(boscia_time, geom_shifted_mean(df_D_corr[ind_min_time_D_corr,:timeBoscia]))
        push!(scip_time, 0)
        push!(pajarito_time, geom_shifted_mean(df_D_corr[ind_min_time_D_corr,:timePajarito]))
        push!(custombb_time, geom_shifted_mean(df_D_corr[ind_min_time_D_corr, :timeCustomBB]))


        ## number nodes and Cuts
        println("Cuts and Nodes")
        push!(boscia_nodes, custom_mean(df_AF_ind[ind_min_time_AF_ind, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_AF_ind[ind_min_time_AF_ind,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_AF_ind[ind_min_time_AF_ind,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_AF_ind[ind_min_time_AF_ind,:numberIterPajarito]))
      #  push!(custombb_nodes, 0)
        push!(custombb_nodes,  custom_mean(df_AF_ind[ind_min_time_AF_ind,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_AF_corr[ind_min_time_AF_corr, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_AF_corr[ind_min_time_AF_corr,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_AF_corr[ind_min_time_AF_corr,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_AF_corr[ind_min_time_AF_corr,:numberIterPajarito]))
        #push!(custombb_nodes, 0)
        push!(custombb_nodes,  custom_mean(df_AF_corr[ind_min_time_AF_corr,:numberNodesCustomBB]))


        push!(boscia_nodes, custom_mean(df_A_ind[ind_min_time_A_ind, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_A_ind[ind_min_time_A_ind,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_A_ind[ind_min_time_A_ind,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_A_ind[ind_min_time_A_ind,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_A_corr[ind_min_time_A_corr, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_A_corr[ind_min_time_A_corr,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_A_corr[ind_min_time_A_corr,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_A_corr[ind_min_time_A_corr,:numberNodesCustomBB]))


        push!(boscia_nodes, custom_mean(df_DF_ind[ind_min_time_DF_ind, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_DF_ind[ind_min_time_DF_ind,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_DF_ind[ind_min_time_DF_ind,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_DF_ind[ind_min_time_DF_ind,:numberIterPajarito]))
       # push!(custombb_nodes, 0)
        push!(custombb_nodes,  custom_mean(df_DF_ind[ind_min_time_DF_ind,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_DF_corr[ind_min_time_DF_corr, :numberNodesBoscia]))
        push!(scip_cuts, custom_mean(df_DF_corr[ind_min_time_DF_corr,:numberCutsScip]))
        push!(pajarito_cuts, custom_mean(df_DF_corr[ind_min_time_DF_corr,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_DF_corr[ind_min_time_DF_corr,:numberIterPajarito]))
       # push!(custombb_nodes, 0)
       push!(custombb_nodes,  custom_mean(df_DF_corr[ind_min_time_DF_corr,:numberNodesCustomBB]))


        push!(boscia_nodes, custom_mean(df_D_ind[ind_min_time_D_ind, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_D_ind[ind_min_time_D_ind,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_D_ind[ind_min_time_D_ind,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_D_ind[ind_min_time_D_ind,:numberNodesCustomBB]))

        push!(boscia_nodes, custom_mean(df_D_corr[ind_min_time_D_corr, :numberNodesBoscia]))
        push!(scip_cuts, 0)
        push!(pajarito_cuts, custom_mean(df_D_corr[ind_min_time_D_corr,:numberCutsPajarito]))
        push!(pajarito_iters, custom_mean(df_D_corr[ind_min_time_D_corr,:numberIterPajarito]))
        push!(custombb_nodes, custom_mean(df_D_corr[ind_min_time_D_corr,:numberNodesCustomBB]))

        ## relative gap 
        println("Relative gap")
        push!(boscia_rel_gap, geo_mean(df_AF_ind[ind_min_time_AF_ind,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_AF_ind[ind_min_time_AF_ind,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_AF_ind[ind_min_time_AF_ind,:relDualGapPajarito]))
     #   push!(custombb_rel_gap, 0)
        push!(custombb_rel_gap, geo_mean(df_AF_ind[ind_min_time_AF_ind,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_AF_corr[ind_min_time_AF_corr,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_AF_corr[ind_min_time_AF_corr,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_AF_corr[ind_min_time_AF_corr,:relDualGapPajarito]))
        #push!(custombb_rel_gap, 0)
        push!(custombb_rel_gap, geo_mean(df_AF_corr[ind_min_time_AF_corr,:relDualGapCustomBB]))


        push!(boscia_rel_gap, geo_mean(df_A_ind[ind_min_time_A_ind,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_A_ind[ind_min_time_A_ind,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_A_ind[ind_min_time_A_ind,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_A_corr[ind_min_time_A_corr,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_A_corr[ind_min_time_A_corr,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_A_corr[ind_min_time_A_corr,:relDualGapCustomBB]))


        push!(boscia_rel_gap, geo_mean(df_DF_ind[ind_min_time_DF_ind,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_DF_ind[ind_min_time_DF_ind,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_DF_ind[ind_min_time_DF_ind,:relDualGapPajarito]))
       # push!(custombb_rel_gap, 0)
        push!(custombb_rel_gap, geo_mean(df_DF_ind[ind_min_time_DF_ind,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_DF_corr[ind_min_time_DF_corr,:relDualGapBoscia]))
        push!(scip_rel_gap, geo_mean(df_DF_corr[ind_min_time_DF_corr,:relDualGapScip]))
        push!(pajarito_rel_gap, geo_mean(df_DF_corr[ind_min_time_DF_corr,:relDualGapPajarito]))
      #  push!(custombb_rel_gap, 0)
      push!(custombb_rel_gap, geo_mean(df_DF_corr[ind_min_time_DF_corr,:relDualGapCustomBB]))


        push!(boscia_rel_gap, geo_mean(df_D_ind[ind_min_time_D_ind,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_D_ind[ind_min_time_D_ind,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_D_ind[ind_min_time_D_ind,:relDualGapCustomBB]))

        push!(boscia_rel_gap, geo_mean(df_D_corr[ind_min_time_D_corr,:relDualGapBoscia]))
        push!(scip_rel_gap, 0)
        push!(pajarito_rel_gap, geo_mean(df_D_corr[ind_min_time_D_corr,:relDualGapPajarito]))
        push!(custombb_rel_gap, geo_mean(df_D_corr[ind_min_time_D_corr,:relDualGapCustomBB]))
        end

        ## relative termination
        println("\n")
        println("Relative termination")
        boscia_term_rel = []
        scip_term_rel = []
        pajarito_term_rel = []
        custombb_term_rel = []
        for (i,_) in enumerate(num_instances)
            if num_instances[i] != 0.0
                push!(boscia_term_rel, convert(Int64, round(boscia_term[i]/num_instances[i]*100)))
                push!(scip_term_rel, convert(Int64, round(scip_term[i]/num_instances[i]*100)))
                push!(pajarito_term_rel, convert(Int64, round(pajarito_term[i]/num_instances[i]*100)))
                push!(custombb_term_rel, convert(Int64, round(custombb_term[i]/num_instances[i]*100)))
            else
                push!(boscia_term_rel, NaN)
                push!(scip_term_rel, NaN)
                push!(pajarito_term_rel, NaN)
                push!(custombb_term_rel, NaN)
            end
        end
        #boscia_term_rel = convert.(Int64, round.(boscia_term./num_instances*100))
        #scip_term_rel = convert.(Int64, round.(scip_term./num_instances*100))
        #pajarito_term_rel = convert.(Int64, round.(pajarito_term./num_instances*100))
        #custombb_term_rel = convert.(Int64, round.(custombb_term./num_instances*100))

        # rounding
        boscia_time = round.(boscia_time, digits=3)
        pajarito_time = round.(pajarito_time, digits=3)
        scip_time = round.(scip_time, digits=3)
        custombb_time = round.(custombb_time, digits=3)

        # rounding
        boscia_nodes = convert.(Int64, round.(boscia_nodes))
        scip_cuts = convert.(Int64, round.(scip_cuts))
        pajarito_cuts = convert.(Int64, round.(pajarito_cuts))
        pajarito_iters = convert.(Int64, round.(pajarito_iters))
        custombb_nodes = convert.(Int64, round.(custombb_nodes))

        # rounding
        boscia_rel_gap = round.(boscia_rel_gap, digits=3)
        pajarito_rel_gap = round.(pajarito_rel_gap, digits=3)
        scip_rel_gap = round.(scip_rel_gap, digits=3)
        custombb_rel_gap = round.(custombb_rel_gap, digits=3)

    
        ## fill data into dataframe
        df[!,:Problem] = repeat(vcat(fill("AF",2),fill("A",2),fill("DF",2),fill("D",2)), outer=5)
        df[!,:Correlation] = repeat(["no","yes"],outer=20)
        df[!,:minTime] = vcat(fill(0,8), fill(10,8), fill(100,8), fill(1000,8), fill(2000,8))
        df[!,:numInstances] = num_instances

        df[!,:BosciaTime] = boscia_time
        df[!,:BosciaTerm] = boscia_term
        df[!,:BosciaTermRel] = boscia_term_rel
        df[!,:BosciaRelGap] = boscia_rel_gap
        df[!,:BosciaNodes] = boscia_nodes

        df[!,:PajaritoTime] = pajarito_time
        df[!,:PajaritoTerm] = pajarito_term
        df[!,:PajaritoTermRel] = pajarito_term_rel
        df[!,:PajaritoRelGap] = pajarito_rel_gap
        df[!,:PajaritoCuts] = pajarito_cuts
        df[!,:PajaritoIter] = pajarito_iters

        df[!,:ScipTime] = scip_time
        df[!,:ScipTerm] = scip_term
        df[!,:ScipTermRel] = scip_term_rel
        df[!,:ScipRelGap] = scip_rel_gap
        df[!,:ScipCuts] = scip_cuts

        @show length(custombb_time)
        @show length(custombb_term)
        @show length(custombb_term_rel)
        @show length(custombb_rel_gap)
        @show length(custombb_nodes)
        df[!,:CustomBBTime] = custombb_time
        df[!,:CustomBBTerm] = custombb_term
        df[!,:CustomBBTermRel] = custombb_term_rel
        df[!,:CustomBBRelGap] = custombb_rel_gap
        df[!,:CustomBBNodes] = custombb_nodes
    
    
        file_name = joinpath(@__DIR__, "csv/summary_by_difficulty.csv")
        CSV.write(file_name, df, append=false)
        return true
end

function build_summary(corr)

    function geo_mean(group)
        prod = 1.0
        n = 0
        for element in group
            # @show element
            if element != Inf
                prod = prod * element
                n += 1
            end
        end
        @show prod, n
        if n == 0
            return Inf
        end
        return prod^(1/n)
    end

    function geom_shifted_mean(xs; shift=big"1.0")
        a = length(xs)  
        n= 0
        prod = 1.0  
        if a != 0 
            for xi in xs
                if xi != Inf 
                    prod = prod*(xi+shift)  
                    n += 1
                end
            end
            return Float64(prod^(1/n) - shift)
        end
        return Inf
    end

    type = corr ? "correlated" : "independent"
    df = DataFrame()
    NumInstances = 50

    # load data
    df_AF = DataFrame(CSV.File(joinpath(@__DIR__, "csv/AF_optimality_" * type * "_non_grouped.csv")))
    df_A = DataFrame(CSV.File(joinpath(@__DIR__, "csv/A_optimality_" * type * "_non_grouped.csv")))
    df_DF = DataFrame(CSV.File(joinpath(@__DIR__, "csv/DF_optimality_" * type * "_non_grouped.csv")))
    df_D = DataFrame(CSV.File(joinpath(@__DIR__, "csv/D_optimality_" * type * "_non_grouped.csv")))

    # we are only interested in entries with termination

    ind_term_AF_Boscia = findall(x-> x==1, df_AF[!,:terminationBoscia])
    ind_term_AF_Scip = findall(x-> x==1, df_AF[!,:terminationScip])
    ind_term_AF_Pajarito = findall(x-> x==1, df_AF[!,:terminationPajarito])
    ind_term_AF_Hypatia_cont = findall(x-> x==1, df_AF[!,:terminationHypatiaCont])
    ind_term_AF_Hypatia_limit = findall(x-> x==1, df_AF[!,:terminationHypatiaLimit])
    ind_term_AF_FrankWolfe_cont = findall(x-> x==1, df_AF[!,:terminationFrankWolfeCont])
    ind_term_AF_FrankWolfe_limit = findall(x-> x==1, df_AF[!,:terminationFrankWolfeLimit])

    ind_term_A_Boscia = findall(x-> x==1, df_A[!,:terminationBoscia])
    ind_term_A_Pajarito = findall(x-> x==1, df_A[!,:terminationPajarito])
    ind_term_A_Hypatia_cont = findall(x-> x==1, df_A[!,:terminationHypatiaCont])
    ind_term_A_Hypatia_limit = findall(x-> x==1, df_A[!,:terminationHypatiaLimit])
    ind_term_A_FrankWolfe_cont = findall(x-> x==1, df_A[!,:terminationFrankWolfeCont])
    ind_term_A_FrankWolfe_limit = findall(x-> x==1, df_A[!,:terminationFrankWolfeLimit])

    ind_term_DF_Boscia = findall(x-> x==1, df_DF[!,:terminationBoscia])
    ind_term_DF_Scip = findall(x-> x==1, df_DF[!,:terminationScip])
    ind_term_DF_Pajarito = findall(x-> x==1, df_DF[!,:terminationPajarito])
    ind_term_DF_Hypatia_cont = findall(x-> x==1, df_DF[!,:terminationHypatiaCont])
    ind_term_DF_Hypatia_limit = findall(x-> x==1, df_DF[!,:terminationHypatiaLimit])
    ind_term_DF_FrankWolfe_cont = findall(x-> x==1, df_DF[!,:terminationFrankWolfeCont])
    ind_term_DF_FrankWolfe_limit = findall(x-> x==1, df_DF[!,:terminationFrankWolfeLimit])

    ind_term_D_Boscia = findall(x-> x==1, df_D[!,:terminationBoscia])
    ind_term_D_Pajarito = findall(x-> x==1, df_D[!,:terminationPajarito])
    ind_term_D_Hypatia_cont = findall(x-> x==1, df_D[!,:terminationHypatiaCont])
    ind_term_D_Hypatia_limit = findall(x-> x==1, df_D[!,:terminationHypatiaLimit])
    ind_term_D_FrankWolfe_cont = findall(x-> x==1, df_D[!,:terminationFrankWolfeCont])
    ind_term_D_FrankWolfe_limit = findall(x-> x==1, df_D[!,:terminationFrankWolfeLimit])

    df[!,:criteria] = ["AF", "A", "DF", "D"]
    boscia_term = []
    scip_term = []
    pajarito_term = []
    hypatia_cont_term = []
    hypatia_limit_term = []
    frankwolfe_cont_term = []
    frankwolfe_limit_term = []

    boscia_time = []
    scip_time = []
    pajarito_time = []
    hypatia_cont_time = []
    hypatia_limit_time = []
    frankwolfe_cont_time = []
    frankwolfe_limit_time = []

    boscia_time_all = []
    scip_time_all = []
    pajarito_time_all = []
    hypatia_cont_time_all = []
    hypatia_limit_time_all = []
    frankwolfe_cont_time_all = []
    frankwolfe_limit_time_all = []

    boscia_term_rel = []
    scip_term_rel = []
    pajarito_term_rel = []
    hypatia_cont_term_rel = []
    hypatia_limit_term_rel = []
    frankwolfe_cont_term_rel = []
    frankwolfe_limit_term_rel = []

    ## termination
    push!(boscia_term, length(ind_term_AF_Boscia))
    push!(scip_term, length(ind_term_AF_Scip))
    push!(pajarito_term, length(ind_term_AF_Pajarito))
    push!(hypatia_cont_term, length(ind_term_AF_Hypatia_cont))
    push!(hypatia_limit_term, length(ind_term_AF_Hypatia_limit))
    push!(frankwolfe_cont_term, length(ind_term_AF_FrankWolfe_cont))
    push!(frankwolfe_limit_term, length(ind_term_AF_FrankWolfe_limit))

    push!(boscia_term, length(ind_term_A_Boscia))
    push!(scip_term, 0)
    push!(pajarito_term, length(ind_term_A_Pajarito))
    push!(hypatia_cont_term, length(ind_term_A_Hypatia_cont))
    push!(hypatia_limit_term, length(ind_term_A_Hypatia_limit))
    push!(frankwolfe_cont_term, length(ind_term_A_FrankWolfe_cont))
    push!(frankwolfe_limit_term, length(ind_term_A_FrankWolfe_limit))

    push!(boscia_term, length(ind_term_DF_Boscia))
    push!(scip_term, length(ind_term_DF_Scip))
    push!(pajarito_term, length(ind_term_DF_Pajarito))
    push!(hypatia_cont_term, length(ind_term_DF_Hypatia_cont))
    push!(hypatia_limit_term, length(ind_term_DF_Hypatia_limit))
    push!(frankwolfe_cont_term, length(ind_term_DF_FrankWolfe_cont))
    push!(frankwolfe_limit_term, length(ind_term_DF_FrankWolfe_limit))

    push!(boscia_term, length(ind_term_D_Boscia))
    push!(scip_term, 0)
    push!(pajarito_term, length(ind_term_D_Pajarito))
    push!(hypatia_cont_term, length(ind_term_D_Hypatia_cont))
    push!(hypatia_limit_term, length(ind_term_D_Hypatia_limit))
    push!(frankwolfe_cont_term, length(ind_term_D_FrankWolfe_cont))
    push!(frankwolfe_limit_term, length(ind_term_D_FrankWolfe_limit))


    ## relative termination
    boscia_term_rel = convert.(Int64, round.(boscia_term./NumInstances*100))
    scip_term_rel = convert.(Int64, round.(scip_term./NumInstances*100))
    pajarito_term_rel = convert.(Int64, round.(pajarito_term./NumInstances*100))
    hypatia_cont_term_rel = convert.(Int64, round.(hypatia_cont_term./NumInstances*100))
    hypatia_limit_term_rel = convert.(Int64, round.(hypatia_limit_term./NumInstances*100))
    frankwolfe_cont_term_rel = convert.(Int64, round.(frankwolfe_cont_term./NumInstances*100))
    frankwolfe_limit_term_rel = convert.(Int64, round.(frankwolfe_limit_term./NumInstances*100))

    ## time for solved instances 
    println("Time solved instances A")
    push!(boscia_time, geom_shifted_mean(df_AF[ind_term_AF_Boscia,:timeBoscia]))
    push!(scip_time, geom_shifted_mean(df_AF[ind_term_AF_Scip,:timeScip]))
    push!(pajarito_time, geom_shifted_mean(df_AF[ind_term_AF_Pajarito,:timePajarito]))
    push!(hypatia_cont_time, geom_shifted_mean(df_AF[ind_term_AF_Hypatia_cont,:timeHypatiaCont]))
    push!(hypatia_limit_time, geom_shifted_mean(df_AF[ind_term_AF_Hypatia_limit,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time, geom_shifted_mean(df_AF[ind_term_AF_FrankWolfe_cont,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time, geom_shifted_mean(df_AF[ind_term_AF_FrankWolfe_limit,:timeFrankWolfeLimit]))
        
    push!(boscia_time, geom_shifted_mean(df_A[ind_term_A_Boscia,:timeBoscia]))
    push!(scip_time, 0)
    push!(pajarito_time, geom_shifted_mean(df_A[ind_term_A_Pajarito,:timePajarito]))
    push!(hypatia_cont_time, geom_shifted_mean(df_A[ind_term_A_Hypatia_cont,:timeHypatiaCont]))
    push!(hypatia_limit_time, geom_shifted_mean(df_A[ind_term_A_Hypatia_limit,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time, geom_shifted_mean(df_A[ind_term_A_FrankWolfe_cont,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time, geom_shifted_mean(df_A[ind_term_A_FrankWolfe_limit,:timeFrankWolfeLimit]))

    push!(boscia_time, geom_shifted_mean(df_DF[ind_term_DF_Boscia,:timeBoscia]))
    push!(scip_time, geom_shifted_mean(df_DF[ind_term_DF_Scip,:timeScip]))
    push!(pajarito_time, geom_shifted_mean(df_DF[ind_term_DF_Pajarito,:timePajarito]))
    push!(hypatia_cont_time, geom_shifted_mean(df_DF[ind_term_DF_Hypatia_cont,:timeHypatiaCont]))
    push!(hypatia_limit_time, geom_shifted_mean(df_DF[ind_term_DF_Hypatia_limit,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time, geom_shifted_mean(df_DF[ind_term_DF_FrankWolfe_cont,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time, geom_shifted_mean(df_DF[ind_term_DF_FrankWolfe_limit,:timeFrankWolfeLimit]))

    push!(boscia_time, geom_shifted_mean(df_D[ind_term_D_Boscia,:timeBoscia]))
    push!(scip_time, 0)
    push!(pajarito_time, geom_shifted_mean(df_D[ind_term_D_Pajarito,:timePajarito]))
    push!(hypatia_cont_time, geom_shifted_mean(df_D[ind_term_D_Hypatia_cont,:timeHypatiaCont]))
    push!(hypatia_limit_time, geom_shifted_mean(df_D[ind_term_D_Hypatia_limit,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time, geom_shifted_mean(df_D[ind_term_D_FrankWolfe_cont,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time, geom_shifted_mean(df_D[ind_term_D_FrankWolfe_limit,:timeFrankWolfeLimit]))

    boscia_time = round.(boscia_time, digits=3)
    pajarito_time = round.(pajarito_time, digits=3)
    scip_time = round.(scip_time, digits=3)
    hypatia_cont_time = round.(hypatia_cont_time, digits=3)
    hypatia_limit_time = round.(hypatia_limit_time, digits=3)
    frankwolfe_cont_time = round.(frankwolfe_cont_time, digits=3)
    frankwolfe_limit_time = round.(frankwolfe_limit_time, digits=3)

    ## time for all instances 
    push!(boscia_time_all, geom_shifted_mean(df_AF[!,:timeBoscia]))
    push!(scip_time_all, geom_shifted_mean(df_AF[!,:timeScip]))
    push!(pajarito_time_all, geom_shifted_mean(df_AF[!,:timePajarito]))
    push!(hypatia_cont_time_all, geom_shifted_mean(df_AF[!,:timeHypatiaCont]))
    push!(hypatia_limit_time_all, geom_shifted_mean(df_AF[!,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time_all, geom_shifted_mean(df_AF[!,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time_all, geom_shifted_mean(df_AF[!,:timeFrankWolfeLimit]))

    push!(scip_time_all, 0)
    push!(boscia_time_all, geom_shifted_mean(df_A[!,:timeBoscia]))
    push!(pajarito_time_all, geom_shifted_mean(df_A[!,:timePajarito]))
    push!(hypatia_cont_time_all, geom_shifted_mean(df_A[!,:timeHypatiaCont]))
    push!(hypatia_limit_time_all, geom_shifted_mean(df_A[!,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time_all, geom_shifted_mean(df_A[!,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time_all, geom_shifted_mean(df_A[!,:timeFrankWolfeLimit]))

    push!(boscia_time_all, geom_shifted_mean(df_DF[!,:timeBoscia]))
    push!(scip_time_all, geom_shifted_mean(df_DF[!,:timeScip]))
    push!(pajarito_time_all, geom_shifted_mean(df_DF[!,:timePajarito]))
    push!(hypatia_cont_time_all, geom_shifted_mean(df_DF[!,:timeHypatiaCont]))
    push!(hypatia_limit_time_all, geom_shifted_mean(df_DF[!,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time_all, geom_shifted_mean(df_DF[!,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time_all, geom_shifted_mean(df_DF[!,:timeFrankWolfeLimit]))

    push!(boscia_time_all, geom_shifted_mean(df_D[!,:timeBoscia]))
    push!(scip_time_all, 0)
    push!(pajarito_time_all, geom_shifted_mean(df_D[!,:timePajarito]))
    push!(hypatia_cont_time_all, geom_shifted_mean(df_D[!,:timeHypatiaCont]))
    push!(hypatia_limit_time_all, geom_shifted_mean(df_D[!,:timeHypatiaLimit]))
    push!(frankwolfe_cont_time_all, geom_shifted_mean(df_D[!,:timeFrankWolfeCont]))
    push!(frankwolfe_limit_time_all, geom_shifted_mean(df_D[!,:timeFrankWolfeLimit]))

    boscia_time_all = round.(boscia_time_all, digits=3)
    pajarito_time_all = round.(pajarito_time_all, digits=3)
    scip_time_all = round.(scip_time_all, digits=3)
    hypatia_cont_time_all = round.(hypatia_cont_time_all, digits=3)
    hypatia_limit_time_all = round.(hypatia_limit_time_all, digits=3)
    frankwolfe_cont_time_all = round.(frankwolfe_cont_time_all, digits=3)
    frankwolfe_limit_time_all = round.(frankwolfe_limit_time_all, digits=3)

    ## integer solutions found
    int_hypatia_cont = []
    int_hypatia_limit = []
    int_frankwolfe_cont = []
    int_frankwolfe_limit = []

    push!(int_hypatia_cont, sum(df_AF[!,:intFoundHypatiaCont]))
    push!(int_hypatia_cont, sum(df_A[!,:intFoundHypatiaCont]))
    push!(int_hypatia_cont, sum(df_DF[!,:intFoundHypatiaCont]))
    push!(int_hypatia_cont, sum(df_D[!,:intFoundHypatiaCont]))
    
    push!(int_hypatia_limit, sum(df_AF[!,:intFoundHypatiaLimit]))
    push!(int_hypatia_limit, sum(df_A[!,:intFoundHypatiaLimit]))
    push!(int_hypatia_limit, sum(df_DF[!,:intFoundHypatiaLimit]))
    push!(int_hypatia_limit, sum(df_D[!,:intFoundHypatiaLimit]))
    
    push!(int_frankwolfe_cont, sum(df_AF[!,:intFoundFrankWolfeCont]))
    push!(int_frankwolfe_cont, sum(df_A[!,:intFoundFrankWolfeCont]))
    push!(int_frankwolfe_cont, sum(df_DF[!,:intFoundFrankWolfeCont]))
    push!(int_frankwolfe_cont, sum(df_D[!,:intFoundFrankWolfeCont]))

    push!(int_frankwolfe_limit, sum(df_AF[!,:intFoundFrankWolfeLimit]))
    push!(int_frankwolfe_limit, sum(df_A[!,:intFoundFrankWolfeLimit]))
    push!(int_frankwolfe_limit, sum(df_DF[!,:intFoundFrankWolfeLimit]))
    push!(int_frankwolfe_limit, sum(df_D[!,:intFoundFrankWolfeLimit]))

    int_hypatia_cont_rel = convert.(Int64, round.(int_hypatia_cont./NumInstances*100))
    int_hypatia_limit_rel = convert.(Int64, round.(int_hypatia_limit./NumInstances*100))
    int_frankwolfe_cont_rel = convert.(Int64, round.(int_frankwolfe_cont./NumInstances*100))
    int_frankwolfe_limit_rel = convert.(Int64, round.(int_frankwolfe_limit./NumInstances*100))


    ## fill data into dataframe
    df[!,:BosciaTerm] = boscia_term
    df[!,:ScipTerm] = scip_term
    df[!,:PajaritoTerm] = pajarito_term
    df[!,:HypatiaContTerm] = hypatia_cont_term
    df[!,:HypatiaLimitTerm] = hypatia_limit_term
    df[!,:FrankWolfeContTerm] = frankwolfe_limit_term
    df[!,:FrankWolfeLimitTerm] = frankwolfe_cont_term

    df[!,:BosciaTime] = boscia_time
    df[!,:ScipTime] = scip_time
    df[!,:PajaritoTime] = pajarito_time
    df[!,:HypatiaContTime] = hypatia_cont_time
    df[!,:HypatiaLimitTime] = hypatia_limit_time
    df[!,:FrankWolfeContTime] = frankwolfe_limit_time
    df[!,:FrankWolfeLimitTime] = frankwolfe_cont_time

    df[!,:BosciaTermRel] = boscia_term_rel
    df[!,:ScipTermRel] = scip_term_rel
    df[!,:PajaritoTermRel] = pajarito_term_rel
    df[!,:HypatiaContTermRel] = hypatia_cont_term_rel
    df[!,:HypatiaLimitTermRel] = hypatia_limit_term_rel
    df[!,:FrankWolfeContTermRel] = frankwolfe_limit_term_rel
    df[!,:FrankWolfeLimitTermRel] = frankwolfe_cont_term_rel

    df[!,:BosciaTimeAll] = boscia_time_all
    df[!,:ScipTimeAll] = scip_time_all
    df[!,:PajaritoTimeAll] = pajarito_time_all
    df[!,:HypatiaContTimeAll] = hypatia_cont_time_all
    df[!,:HypatiaLimitTimeAll] = hypatia_limit_time_all
    df[!,:FrankWolfeContTimeAll] = frankwolfe_limit_time_all
    df[!,:FrankWolfeLimitTimeAll] = frankwolfe_cont_time_all

    df[!,:HypatiaContInt] = int_hypatia_cont
    df[!,:HypatiaLimitInt] = int_hypatia_limit
    df[!,:FrankWolfeContInt] = int_frankwolfe_cont
    df[!,:FrankWolfeLimitInt] = int_frankwolfe_limit

    df[!,:HypatiaContIntRel] = int_hypatia_cont_rel
    df[!,:HypatiaLimitIntRel] = int_hypatia_limit_rel
    df[!,:FrankWolfeContIntRel] = int_frankwolfe_cont_rel
    df[!,:FrankWolfeLimitIntRel] = int_frankwolfe_limit_rel


    file_name = joinpath(@__DIR__, "csv/summary_" * type * ".csv")
    CSV.write(file_name, df, append=false)
    return true
end

#=println("FrankWolfe")
        # load FrankWolfe continuous version
        df_fw = DataFrame(CSV.File(joinpath(@__DIR__, "csv/" * direc * "frank_wolfe_" * criterion * "-cont_optimality.csv")))
        df_fw.termination .= replace.(df_fw.termination, "Time limit reached" => "TIME_LIMIT")
        for row in eachrow(df_fw)
            if row.time > 1800
                row.termination = "TIME_LIMIT" 
            end
        end

        termination_fw = [row == "OPTIMAL" ? 1 : 0 for row in df_fw[!,:termination]]
        int_found = [row1 != Inf && row2 == true ? 1 : 0 for (row1, row2) in zip(df_fw[!, :solution_int], df_fw[!,:feasibilty])]

        df_fw[!,:timeFrankWolfeCont] = df_fw[!,:time]
        df_fw[!,:solutionFrankWolfeContInt] = df_fw[!,:solution_int]
        df_fw[!,:solutionFrankWolfeCont] = df_fw[!,:solution_fw]
        df_fw[!,:terminationFrankWolfeCont] = termination_fw
        df_fw[!, :dualGapFrankWolfeCont] = df_fw[!, :solution_int] - lowerBounds
        df_fw[!,:intFoundFrankWolfeCont] = int_found
        df_fw = select(df_fw, [:terminationFrankWolfeCont, :timeFrankWolfeCont, :solutionFrankWolfeCont, :solutionFrankWolfeContInt, :dualGapFrankWolfeCont, :intFoundFrankWolfeCont, :seed, :numberOfExperiments, :numberOfParameters])

        df = outerjoin(df, df_fw, on = [:seed, :numberOfExperiments, :numberOfParameters])

        # load FrankWolfe limit version
        df_fw = DataFrame(CSV.File(joinpath(@__DIR__, "csv/" * direc * "frank_wolfe_" * criterion * "-limit_optimality.csv")))
        df_fw.termination .= replace.(df_fw.termination, "Time limit reached" => "TIME_LIMIT")
        for row in eachrow(df_fw)
            if row.time > 1800
                row.termination = "TIME_LIMIT" 
            end
        end

        termination_fw = [row == "OPTIMAL" ? 1 : 0 for row in df_fw[!,:termination]]
        int_found = [row1 != Inf && row2 == true ? 1 : 0 for (row1, row2) in zip(df_fw[!, :solution_int], df_fw[!,:feasibilty])]

        df_fw[!,:timeFrankWolfeLimit] = df_fw[!,:time]
        df_fw[!,:solutionFrankWolfeLimitInt] = df_fw[!,:solution_int]
        df_fw[!,:solutionFrankWolfeLimit] = df_fw[!,:solution_fw]
        df_fw[!,:terminationFrankWolfeLimit] = termination_fw
        df_fw[!, :dualGapFrankWolfeLimit] = df_fw[!, :solution_int] - lowerBounds 
        df_fw[!, :intFoundFrankWolfeLimit] = int_found
        df_fw = select(df_fw, [:terminationFrankWolfeLimit, :timeFrankWolfeLimit, :solutionFrankWolfeLimit, :solutionFrankWolfeLimitInt, :dualGapFrankWolfeLimit, :intFoundFrankWolfeLimit, :seed, :numberOfExperiments, :numberOfParameters])

        df = outerjoin(df, df_fw, on = [:seed, :numberOfExperiments, :numberOfParameters])=#

        #=  # load hypatia cont
        df_hyp = DataFrame(CSV.File(joinpath(@__DIR__, "csv/" * direc * "hypatia_" * criterion * "-cont_optimality.csv")))
        #df_hyp.termination .= replace.(df_hyp.termination, "Time limit reached" => "TIME_LIMIT")
        for row in eachrow(df_hyp)
            if row.time > 1800
                row.termination = "TIME_LIMIT" 
            end
        end

        termination_hyp = [row == "OPTIMAL" ? 1 : 0 for row in df_hyp[!,:termination]]
        int_found = [row1 != Inf && row2 == true ? 1 : 0 for (row1, row2) in zip(df_hyp[!, :solution_int], df_hyp[!,:feasible])]

        df_hyp[!,:timeHypatiaCont] = df_hyp[!,:time]
        df_hyp[!,:solutionHypatiaCont] = df_hyp[!,:solution]
        df_hyp[!,:solutionHypatiaContInt] = df_hyp[!,:solution_int]
        df_hyp[!,:terminationHypatiaCont] = termination_hyp
        df_hyp[!,:dualGapHypatiaCont] = df_hyp[!, :solution_int] - lowerBounds
        df_hyp[!,:intFoundHypatiaCont] = int_found
        df_hyp = select(df_hyp, [:terminationHypatiaCont, :timeHypatiaCont, :solutionHypatiaCont, :solutionHypatiaContInt, :dualGapHypatiaCont, :intFoundHypatiaCont, :seed, :numberOfExperiments, :numberOfParameters])

        df = outerjoin(df, df_hyp, on = [:seed, :numberOfExperiments, :numberOfParameters])

        # load hypatia limit
        df_hyp = DataFrame(CSV.File(joinpath(@__DIR__, "csv/" * direc * "hypatia_" * criterion * "-limit_optimality.csv")))
        #df_hyp.termination .= replace.(df_hyp.termination, "Time limit reached" => "TIME_LIMIT")
        for row in eachrow(df_hyp)
            if row.time > 1800
                row.termination = "TIME_LIMIT" 
            end
        end

        termination_hyp = [row == "OPTIMAL" ? 1 : 0 for row in df_hyp[!,:termination]]
        int_found = [row1 != Inf && row2 == true ? 1 : 0 for (row1, row2) in zip(df_hyp[!, :solution_int], df_hyp[!,:feasible])]

        df_hyp[!,:timeHypatiaLimit] = df_hyp[!,:time]
        df_hyp[!,:solutionHypatiaLimit] = df_hyp[!,:solution]
        df_hyp[!,:solutionHypatiaLimitInt] = df_hyp[!,:solution_int]
        df_hyp[!,:terminationHypatiaLimit] = termination_hyp
        df_hyp[!,:dualGapHypatiaLimit] = df_hyp[!, :solution_int] - lowerBounds
        df_hyp[!,:intFoundHypatiaLimit] = int_found
        df_hyp = select(df_hyp, [:terminationHypatiaLimit, :timeHypatiaLimit, :solutionHypatiaLimit, :solutionHypatiaLimitInt, :dualGapHypatiaLimit, :intFoundHypatiaLimit, :seed, :numberOfExperiments, :numberOfParameters])

        df = outerjoin(df, df_hyp, on = [:seed, :numberOfExperiments, :numberOfParameters]) =#

         #=  if row.intFoundHypatiaCont == 0
                push!(rel_gap_hypatia_cont, "-")
              elseif min(abs(row.solutionHypatiaContInt), abs(row.lbBoscia)) == 0
                push!(rel_gap_hypatia_cont, round(row.solutionHypatiaContInt - row.lbBoscia, digits=3))
              elseif sign(row.lbBoscia) != sign(row.solutionHypatiaContInt)
                push!(rel_gap_hypatia_cont, Inf)
              else 
                push!(rel_gap_hypatia_cont, round((row.solutionHypatiaContInt - row.lbBoscia)/min(abs(row.solutionHypatiaContInt), abs(row.lbBoscia)), digits=3))
              end
              if row.intFoundHypatiaLimit == 0
                push!(rel_gap_hypatia_limit, "-")
              elseif  min(abs(row.solutionHypatiaLimitInt), abs(row.lbBoscia)) == 0
                push!(rel_gap_hypatia_limit, round(row.solutionHypatiaLimitInt - row.lbBoscia, digits=3))
              elseif sign(row.lbBoscia) != sign(row.solutionHypatiaLimitInt)
                push!(rel_gap_hypatia_limit, Inf)
              else 
                push!(rel_gap_hypatia_limit, round((row.solutionHypatiaLimitInt - row.lbBoscia)/min(abs(row.solutionHypatiaLimitInt), abs(row.lbBoscia)), digits=3))
              end

              if row.intFoundFrankWolfeCont == 0
                push!(rel_gap_frankwolfe_cont, "-")
              elseif  min(abs(row.solutionFrankWolfeContInt), abs(row.lbBoscia)) == 0
                push!(rel_gap_frankwolfe_cont, round(row.solutionFrankWolfeContInt - row.lbBoscia, digits=3))
              elseif sign(row.lbBoscia) != sign(row.solutionFrankWolfeContInt)
                push!(rel_gap_frankwolfe_cont, Inf)
              else 
                push!(rel_gap_frankwolfe_cont, round((row.solutionFrankWolfeContInt - row.lbBoscia)/min(abs(row.solutionFrankWolfeContInt), abs(row.lbBoscia)), digits=3))
              end
              if row.intFoundFrankWolfeLimit == 0
                push!(rel_gap_frankwolfe_limit, "-")
              elseif min(abs(row.solutionFrankWolfeLimitInt), abs(row.lbBoscia)) == 0
                push!(rel_gap_frankwolfe_limit, round(row.solutionFrankWolfeLimitInt - row.lbBoscia, digits=3))
              elseif sign(row.lbBoscia) != sign(row.solutionFrankWolfeLimitInt)
                push!(rel_gap_frankwolfe_limit, Inf)
              else 
                push!(rel_gap_frankwolfe_limit, round((row.solutionFrankWolfeLimitInt - row.lbBoscia)/min(abs(row.solutionFrankWolfeLimitInt), abs(row.lbBoscia)), digits=3))
              end=#
