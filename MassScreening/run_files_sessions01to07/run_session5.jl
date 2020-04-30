push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./MassScreening")
    include("MS_functions.jl")

###########
β = 2.5
ϵ=.5
session=05
scenarios=[i  for i=1:90]#[8]
n_traj=100
nb_months=1
MS_strategies=[zeros(370)  for i=1:90]
    #for day=30*6+1:30*7      MS_strategies[1][day]=1e6; end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i][day]=1e4;      end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+10][day]=2.5e4;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+20][day]=5e4;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+30][day]=7.5e4;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+40][day]=1e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+50][day]=2.5e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+60][day]=5e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+70][day]=7.5e5;    end end
    for i=1:10      for day=(i-1)*30+1:(i-1)*30+30*nb_months     MS_strategies[i+80][day]=1e6;    end end
    #for i=1:10      for day=(i-1)*30+1:(i-1)*30+30     print("\ti=",i,"/day=",day);  end; println(); end
folder="./MassScreening/session"*string(session)*"_"*string(n_traj)*"sims/";
run_set_scenarios(folder,session,scenarios,β,ϵ,n_traj,MS_strategies)
