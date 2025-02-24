

@with_kw mutable struct CoVParameters
    β::Float64 = 2.5/3.6
    γ::Float64 = 1/3.6
    σ::Float64 = 1/2.
    δ::Float64 = 0.05#Proportion of symptomatic/diseased vs non-symptomatic cases
    τ::Float64 = 1/15. #treatment/isolation rate for symptomatics
    μ₁::Float64 = 0.01#Excess mortality due to disease
    ϵ::Float64 = 0.1 #Relative infectiousness of undetectable infecteds
    ρ::Vector{Float64} = [0.01 for i in 1:n] #spatial coupling
    T::Matrix{Float64} = zeros(n,n)#transmission matrix
    ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    into_mom::Vector{Int}#Number of people flying into Mombassa each day
    into_nai::Vector{Int}#Number of people flying into Nairobi each day
    global_prev::Vector{Float64}
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Vector{Float64} = zeros(n) #For inplace calculations
    N̂::Vector{Float64} = zeros(n)#For inplace calculations
    λ_urb::Vector{Float64} = zeros(n)#For inplace calculations
    λ_rur::Vector{Float64} = zeros(n)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n*n_t)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n*n_t)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n*n_s*2,n*n_t))#For inplace calculations
end

"""
    mutable struct CoVParameters_AS

Struct for containing relevant epidemilogical parameters for the age-structured version of KenyaCoV
"""
@with_kw mutable struct CoVParameters_AS
    #Epidemiological parameters and social contact/mixing rates
    β::Float64 = 1. #Basic transmission probability per contact OR infectious contact rate (can be greater than one depending on scaling)
    c_t::Function = t -> 1. #Time varying basic contact rate
    γ::Float64 = 1/2. #recovery rate for mild and asymptomatic cases (Consensus estimate mean 4 days infectious, default is 1/σ₂ + 1/γ = 4 days)
    σ₁::Float64 = 1/5. #end of latency rate (Consensus estimate mean 5 days)
    σ₂::Float64 = 1/2. #end of pre-symptomatic rate (Consensus estimate mean 4 days infectious, default is 1/σ₂ + 1/γ = 4 days)
    δ::Float64 = 0.9#Proportion of over 80s who get identified
    τ::Float64 = 1/5. #Rate of hospitalisation treatment, conditional on eventually needing it (V category)
    τ_initial::Float64 = 0. # isolation rate for symptomatics at beginning of epidemic
    clear_quarantine = 0. # Average time to end isolation NOT USED IN THIS VERSION
    μ₁::Float64 = 0.0#Excess mortality due to disease NOT USED IN THIS VERSION
    ϵ::Float64 = 0.1 #Relative infectiousness of undetectable/undetected infecteds both pre-symptomatic and asymptomatic
    ϵ_D::Float64 = 1.#Relative infectiousness of mild infecteds due to social avoidance
    ϵ_V::Float64 = 0.4#Relative infectiousness of mild, then severe, infecteds due to social avoidance
    rel_detection_rate::Vector{Float64} = ones(n_a) #relative symptomatic rate
    hₐ::Vector{Float64} = zeros(n_a)   #proportion of severe cases if symptomatic
    ICUₐ::Vector{Float64} = zeros(n_a) #proportion of hospitalised cases that become critical
    χ::Vector{Float64} = ones(n_a) #relative susceptibility
    ρ::Vector{Float64} = [0.01 for i in 1:n] #Time spent outside of area
    T::Matrix{Float64} = zeros(n,n)#Probability distributions of location for mobile individuals
    M::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix
    M_ho::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix --- if only home contacts are made
    ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    into_mom::Vector{Int}#Number of people flying into Mombassa each day
    into_nai::Vector{Int}#Number of people flying into Nairobi each day
    global_prev::Vector{Float64}#Chance that each person arriving was infected on that day
    #Control variables
    isolating_detecteds::Bool = false #This determines if people are still being isolated
    lockdown::Bool = false  #This determines if social distancing and travel restrictions are still in force
    #Calculation variables
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Matrix{Float64} = zeros(n_wa,n_a) #For inplace calculations
    N̂::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    λ::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    λ_loc::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n_wa*n_a*n_ta)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n_wa*n_a*n_ta)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n_wa*n_a*n_s,n_wa*n_a*n_ta))#For inplace calculations
    du_linear::Vector{Int64} = zeros(Int64,n_wa*n_a*n_s) #for inplace calculations
end
