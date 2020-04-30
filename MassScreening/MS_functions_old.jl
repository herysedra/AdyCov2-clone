
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV_MS
using LinearAlgebra:eigen
using Statistics: median, quantile

###########
function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end
function output_daily_and_final_incidence(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end] #time step = 1 day
    cumIs = [sum(sol(t)[:,:,9:11],dims = 2)[:,1,:]  for t in times] # only cumulative A, M and V
    #Q=[[sum(sol(t)[wa,:,5],dims = 1)[:,1,:] #=+ sum(sol(t)[wa,:,11],dims = 1)[:,1,:]=#  for wa=1:20]  for t in times]      # Q = Q ##+ QS
    I=[[sum(sol(t)[wa,:,4:6],dims = 1)[:,1,:] #=+ sum(sol(t)[wa,:,5],dims = 1)[:,1,:] + sum(sol(t)[wa,:,6],dims = 1)[:,1,:]=#  for wa=1:20]  for t in times]      # I = A+M+V
    return [cumIs,I,sol[end][:,:,9:11]],false # save z (time series with only incidence with no age structure), and save the final distribution (end) age/space but no time
end

function run_set_scenarios(folder,session,scenarios,R₀,n_traj,dt,MS_strategies)
    if !isdir(folder)   mkdir(folder)   end
    global sims_vector=[]
    for i=1:size(scenarios,1)
        sc_nb=session*100+scenarios[i]
        u0,P,P_dest = KenyaCoV_MS.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        u0[KenyaCoV_MS.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
        P.dt = dt;     #P.ext_inf_rate = 0.;    P.ϵ = ϵ;
        P.δ = .1#δ;      P.γ = γ;      P.σ = σ;        #P.β = β;       P.μ₁=μ₁;                P.τ=1/2.;
        P.β=R₀*P.γ/(P.δ + P.ϵ*(1-P.δ))
        println("P.β=",P.β)
        P.MS_strategy=MS_strategies[i]
        prob = KenyaCoV_MS.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=[]
        for i=1:size(sims.u,1)
            push!(sims_vector, sims.u[i])
        end
        @save folder*"sims_sc"*string(sc_nb)*".jld2" sims_vector
    end
    plot([sum(sims_vector[1][1][t])     for t = 1:1:365])
end
