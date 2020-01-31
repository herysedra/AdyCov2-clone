push!(LOAD_PATH, "./src")
using Plots,Parameters,Distributions
using KenyaCoV
#Load data and completely susceptible Population
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv",0.01)
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl = solve_KenyaCoV_prob(u0,(0.,60.),P,0.25)



ũ = [reshape(u,n,n_s,2) for u in sol_tl.u]

susceptibles = [sum(u[:,1,:])/sum(u) for u in ũ ]
infecteds_A = [sum(u[:,3,:]) for u in ũ ]
infecteds_D = [sum(u[:,4,:]) for u in ũ ]
recovereds = [sum(u[:,5,:])/sum(u) for u in ũ ]
cum_infecteds =  [sum(u[:,7:8,:])/sum(u) for u in ũ ]

plot(sol_tl.t,susceptibles,lab="S")
plot(sol_tl.t,infecteds_A,lab ="I_A")
plot!(sol_tl.t,infecteds_D,lab ="I_D")
plot!(sol_tl.t,recovereds,lab="R")
plot!(sol_tl.t,infecteds_D,lab ="I_D")
plot!(sol_tl.t,cum_infecteds,lab ="cum. I")
