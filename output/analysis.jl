
push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,KenyaCoV,DifferentialEquations,StatsPlots,JLD2,FileIO,MAT,RecursiveArrayTools

u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                             "data/optimal_transition_matrix.jld2",
                                                            "data/optimal_movement_matrix.jld2",
                                                            "data/flight_numbers.csv",
                                                            "data/projected_global_prevelance.csv")
"""
Output functions: Target the peak timing for each county, cases by each county,
timing of total peak and total cases
"""
function output_infecteds_by_county(sol,i)
    [sum(reshape(u,n,n_s,2)[:,3:4,1:2,],dims=[2,3])[:,1,1] for u in sol.u],false
end

function output_infecteds_and_cum_by_county(sol,i)
    z = [sum(reshape(u,n,n_s,2)[:,[3,4,7,8],1:2,],dims=[3])[:,:,1] for u in sol.u]
    l = length(sol.u)
    (sol.t[1:10:l],z[1:10:l]),false#Thin output to daily
end

"""
SCENARIO 1: R/J estimates no intervention

Current estimates from Read and Jewell --- assuming that ascertainment is a decent estimate for symptomatic rate

scenario where asymptomatics are just as infectious as symptomatics, this broadly matches Jewell/Read where
undetected infecteds are just as infectious. No need for real-time growth rate matching.
After the first infected assume no more introductions
"""

P.τ = 0.;#No treatment
P.σ = 1/4; #Latent period mean 4 days
P.γ = 1/1.61; #Fast infectious duration
P.β = 1.94;
P.ϵ = 1.;
P.ext_inf_rate = 0.;
u0[30,3,1] += 1#One asymptomatic in Nairobi
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,180.),P)

CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func = output_infecteds_and_cum_by_county)
sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.1,trajectories = 1000)



times = sim[4][1]
epidemic = VectorOfArray(sim[4][2])
size(epidemic)
nai_infs = sum(epidemic[30,1:2,:],dims = 1)[:]
plot!(times,nai_infs)
xlabel!("time (days)")
u0[30,1,1:2]
nai_peak = zeros(500)
for i = 1:500
    epidemic = VectorOfArray(sim[i][2])
    nai_infs = sum(epidemic[30,1:2,:],dims = 1)[:]
    nai_peak[i] = times[argmax(nai_infs)]
end
median(nai_peak)
boxplot(nai_peak)



for i = 1:500
    y = sum(forecasts[i][:,:],dims =1)[:]
    total_peaktimes[i] = times[argmax(y)]
end
@save "output/total_peaktimes.jld2" total_peaktimes
matwrite("output/total_peaktimes.mat",Dict("total_peaktimes"=>total_peaktimes))

peaktimes_by_county = zeros(500,47)
for i = 1:500,j=1:47
    y = forecasts[i][j,:]
    peaktimes_by_county[i,j] = times[argmax(y)]
end
@save "output/peaktimes_by_county.jld2" peaktimes_by_county
matwrite("output/peaktimes_by_county.mat",Dict("peaktimes_by_county"=>peaktimes_by_county))



"""
SCENARIO 2: Liu et al estimates no intervention (https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1)
    No treatment and everyone equally infectious
    Randomised parameters
"""
function randomise_params(jumpprob,i,repeat)
    P = jumpprob.prob.p
    P.σ = 1/exp(rand(Normal(log(5.2),0.35))) #1/5.2; #Latent period mean prediction 5.2 days
    P.γ = 1/rand(Uniform(2,3)); #Fast infectious duration
    R₀ = rand(Gamma(100,2.2/100))
    P.β = R₀*P.γ
    return jumpprob
end

P.τ = 1/7.;#No treatment
P.σ = 1/exp(rand(Normal(log(5.2),0.35))) #1/5.2; #Latent period mean prediction 5.2 days
P.γ = 1/rand(Uniform(2,3)); #Fast infectious duration
R₀ = rand(Gamma(100,2.2/100))
P.β = R₀*P.γ
P.ϵ = 1.;
P.ext_inf_rate = 0.;
randomise_params(jump_prob_tl,1,1)
u0[30,3,1]  += 1#One asymptomatic in Nairobi
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,180.),P)

CoV_ens_prob = EnsembleProblem(jump_prob_tl,
                output_func = output_infecteds_and_cum_by_county,
                prob_func = randomise_params)
sim_liu_treatment = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.1,trajectories = 1000)

@save "/Users/Sam/Documents/sim_liu_treatment.jld2" sim_liu_treatment
times = sim_liu_no_treatment[1][1]

sim_liu_no_treatment = 0
"""
SCENARIO 3: Liu et al estimates no intervention (https://www.biorxiv.org/content/10.1101/2020.01.25.919787v1)
    No treatment and everyone equally infectious
    central estiamted parameters
"""

P.τ = 0.;#No treatment
P.σ = 1/5.2 #1/5.2; #Latent period mean prediction 5.2 days
P.γ = 1/2.5; #Fast infectious duration
R₀ =  2.2
P.β = R₀*P.γ
P.ϵ = 1.;
P.ext_inf_rate = 0.;

u0[30,3,1] += 1#One asymptomatic in Nairobi
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)

CoV_ens_prob = EnsembleProblem(jump_prob_tl,
                output_func = output_infecteds_and_cum_by_county)
sim_liu_no_treatment = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.1,trajectories = 1)
sim = sim_liu_no_treatment.u
y = [x[2] for x in sim]




times = sim_liu_no_treatment[1][1]
peaks = zeros(100,47)
for i = 1:100
    epidemic = VectorOfArray(sim_liu_no_treatment[i][2])
    for j = 1:47
        y = sum(epidemic[j,1:2,:],dims = 1)[:]
        peaks[i,j] = times[argmax(y)]
    end
end



median(nai_peak)
boxplot(nai_peak)


P.τ = 1/7.;#No treatment
P.δ = 0.05
CoV_ens_prob = EnsembleProblem(jump_prob_tl,
                output_func = output_infecteds_and_cum_by_county)

sim_liu_with_treatment = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 100)


#
@load "/Users/Sam/Documents/sim_liu_treatment.jld2" sim_liu_treatment
@load "/Users/Sam/Documents/sim_liu_no_treatment.jld2" sim_liu_no_treatment
total_cases = [sum(s[2][end][:,3:4]) for s in sim_liu_no_treatment]
total_cases_treatment = [sum(s[2][end][:,3:4]) for s in sim_liu_treatment]

fig = boxplot(total_cases,lab = "No treatment")
boxplot!(fig, total_cases_treatment,lab = "7 days post-symptoms isolation",
                ylims = (0.,6e7),xticks = [],
                ylabel = "total cases")
plot!([1,2],[sum(u0),sum(u0)])
savefig(fig,"totalcases.png")

VectorOfArray(sim_liu_no_treatment[][2])
times_no_t = sim_liu_no_treatment[1][1]
forecasts_no_t = [VectorOfArray(s[2])[:,3:4,:] for s in sim_liu_no_treatment]

total_epidemic = zeros(length(times_no_t),1000)
for i = 1:1000
    total_epidemic[:,i] = sum(forecasts_no_t[i],dims = [1,2])[:]
end
plot(times_no_t[2:end],diff(total_epidemic[:,2]),lab = "",ylabel = "cum. infs")




for i = 1:500
    y = sum(forecasts[i][:,:],dims =1)[:]
    total_peaktimes[i] = times[argmax(y)]
end
@save "output/total_peaktimes.jld2" total_peaktimes
matwrite("output/total_peaktimes.mat",Dict("total_peaktimes"=>total_peaktimes))

peaktimes_by_county = zeros(500,47)
for i = 1:500,j=1:47
    y = forecasts[i][j,:]
    peaktimes_by_county[i,j] = times[argmax(y)]
end
@save "output/peaktimes_by_county.jld2" peaktimes_by_county
matwrite("output/peaktimes_by_county.mat",Dict("peaktimes_by_county"=>peaktimes_by_county))
