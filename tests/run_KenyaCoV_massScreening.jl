#push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./MassScreening")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
import KenyaCoV_MS
using LinearAlgebra:eigen,normalize


"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV_MS.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
N = sum(u0[:,:,1],dims = 1)


@load "data/agemixingmatrix_china.jld2" M_China
@load "data/agemixingmatrix_Kenya_norestrictions.jld2" M_Kenya
@load "data/agemixingmatrix_Kenya_homeonly.jld2" M_Kenya_ho


"""
Can adjust β to match a desired R₀ value by evaluating the leading eigenvalue
The idea is to match to the chinese epidemic R₀ -- it will be different in Kenya
"""

P.ϵ = 0.1
P.χ = ones(KenyaCoV_MS.n_a)
#R₀_scale = KenyaCoV_MS.calculate_R₀_scale(P)
#P.χ = ones(KenyaCoV_MS.n_a)/R₀_scale
P.β = 2*P.γ

#KenyaCoV_MS.calculate_R₀(P)
#KenyaCoV_MS.calculate_R₀_homeonly(P)

P.dt = 0.25

"""
Run model

"""

u0[KenyaCoV_MS.ind_nairobi_as,5,4] = 5#10 diseased
P.ϵ_D = 1
function ramp_down(t)
    if t < 60.
        return (1-t/60) + 0.5*t/60
    else
        return 0.5
    end
end
P.c_t = ramp_down
P.c_t = t -> 1.
prob = KenyaCoV_MS.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
@time sol = solve(prob,FunctionMap(),dt = P.dt)

cum_inc = [sum(sol(t)[:,:,9:11]) for t = 0.:1.:365]
plot(cum_inc)
plot(diff(cum_inc))
#plot!(diff(cum_inc).+1,yscale = :log10)


prob_ode = KenyaCoV_MS.create_KenyaCoV_ode_prob(u0,(0.,365.),P)
@time sol_ode = solve(prob_ode,Tsit5())

times = 0:1:365
I_area = zeros(Int64,20,length(sol.t))
for i = 1:20,(ind,t) in enumerate(sol.t)
    I_area[i,ind] = sum(sol(t)[i,:,3:4])
end
# I = [sum(sol(t)[:,:,3:4]) for t in sol.t]
plt = plot(sol.t,I_area[1,:], lab = 1,xlims = (0.,30),ylims = (0.,100));
for i = 2:20
    plot!(plt,sol.t,I_area[i,:],lab = i);
end
display(plt)
