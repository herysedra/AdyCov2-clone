include("plots_MS.jl")
    push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./MassScreening")
    include("MS_functions.jl")

sc=01;sims=string(1);
    session1=string(sc);session2=string(sc+1);session3=string(sc+2);session4=string(sc+3);
##########  Verify curves
folders=["./MassScreening/session"*session1*"_"*sims*"sims/"]#,"./MassScreening/session"*session2*"_"*sims*"sims/"]
    plot_folder="./MassScreening/session"*session1*"_"*sims*"sims/testcurves/";
    curves1(folders,plot_folder)

####### Mass screening
    folders=["./MassScreening/session"*session1*"_"*sims*"sims/"]
    plot_folder=folders[1]*"ComparisonPlots/"
    #plot_folder="C:/Users/rabia/Documents/JUNO/KenyaCoV-master_2020-04-23_MassScreening/MassScreening/session1_50sims/ComparisonPlots/"
    cumIs,cumIs_4=make_plots_massScreening(folders,plot_folder)
