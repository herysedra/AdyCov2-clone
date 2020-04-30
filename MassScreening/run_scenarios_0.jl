push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./MassScreening")
    include("../MS_functions.jl")

###########
r_R₀=2.5; δ=.05
ϵ = .3; σ = .2; γ = 1/2.5; β = r_R₀*γ/(δ + ϵ*(1-δ))
println("\nR₀,ϵ,δ=",r_R₀,",",ϵ,",",δ,"\n")

###########
session=06
scenarios=[i  for i=1:90]
n_traj=50
nb_months=3
testing_strategies=[zeros(370)  for i=1:90]
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i][day]=1e4;      end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+10][day]=2.5e4;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+20][day]=5e4;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+30][day]=7.5e4;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+40][day]=1e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+50][day]=2.5e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+60][day]=5e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+70][day]=7.5e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     testing_strategies[i+80][day]=1e6;    end end
    #for i=1:10      for day=(i-1)*30+1:(i-1)*30+30     print("\ti=",i,"/day=",day);  end; println(); end

folder="./MassScreening/session"*string(session)*"_"*string(n_traj)*"sims/";
run_set_scenarios(folder,session,scenarios,ϵ,σ,γ,δ,β,r_R₀,n_traj,.5,testing_strategies)
