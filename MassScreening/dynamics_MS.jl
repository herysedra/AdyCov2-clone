"""
Basic representation of the state of the age structured model:
u[wa_index,age_group,disease state]

This is row major unpacked so
first 1,...,n_wa entries are 0-4 year old susceptibles in wider areas 1,...,n_wa
then n_wa+1,...,2n_wa are 5-9 year old susceptibles in wider areas 1,...,n_wa
.
.
.
then n_wa*n_s + 1, ...,  n_wa*n_s + n_wa entries are 0-4 year old exposed in wider areas 1,...,n_wa
.
.
.

States:
1 -> S(usceptibles)
2 -> E(xposed/latent infected)
3 -> P(re-symptomatic infected)
4 -> A(symptomatic)
5 -> M(ild) symptomatics
6 -> first mild then eventually (se)V(ere) symptomatics
7 -> H(ospitalised)
8 -> Recovered
9 -> Cumulative P->A
10-> Cumulative P->M
11-> Cumulative P->V
12 -> Cumulative V->H
13 -> Q(uarantined)
14 -> Q(uarantined) (se)V(ere) Qᵥ

Events for each wider area and age group:

1-> S to E
2-> E to P
3-> P to A
4-> P to M
5-> P to V
6-> V to H
7-> M to R
8-> A to R
9-> E to Q
10->P to Q
11->A to Q
12->M to Q
13->V to Q
14->Q to R
15->Q to H

"""

dc_age = zeros(Int64,n_wa*n_a*n_s,n_ta*n*n_a)

function import_rate_mom(t,into_mom,global_prev)
    if t+1>min(length(into_mom),length(global_prev))
        t_int=min(length(into_mom),length(global_prev))
    else
        t_int=Int(floor(t))+1
    end
    return into_mom[t_int]*global_prev[t_int]
end

function import_rate_nai(t,into_nai,global_prev)
    if t+1>min(length(into_nai),length(global_prev))
        t_int=min(length(into_nai),length(global_prev))
    else
        t_int=Int(floor(t))+1
    end
    return into_nai[t_int]*global_prev[t_int]
end


# asymp_indices = zeros(Bool,n_wa,n_a,n_s)
# asymp_indices[:,:,3] .= true;
# f_asymp_indices = findall(asymp_indices[:])
# diseased_indices = zeros(Bool,n_wa,n_a,n_s)
# diseased_indices[:,:,4] .= true;
# f_diseased_indices = findall(diseased_indices[:])
# asymp_indices = 0;#free memory
# diseased_indices = 0;
"""
    function calculate_infection_rates!(u,p::CoVParameters_AS,t)

Inplace method for calculating the force of infection on susceptibles in each spatial and age group.
"""
function calculate_infection_rates!(u,p::CoVParameters_AS,t)
    I_P = @view u[:,:,3]
    I_A = @view u[:,:,4]
    I_M = @view u[:,:,5]
    I_V = @view u[:,:,6]
    mul!(p.Î,p.T,p.ϵ*I_P .+ p.ϵ*I_A .+ p.ϵ_D*I_M .+ p.ϵ_V*I_V)  #Local infecteds **if** everyone moved around
    p.Î[:,immobile_age_indices] .= p.ϵ*I_P[:,immobile_age_indices] .+ p.ϵ*I_A[:,immobile_age_indices] .+ p.ϵ_D*I_M[:,immobile_age_indices] .+ p.ϵ_V*I_V[:,immobile_age_indices]#This corrects for immobility
    mul!(p.λ_loc,p.β*p.c_t(t).*(p.Î ./p.N̂),p.M)#Local force of infection due to age-mixing --- M is in to (row), from (col) format
    mul!(p.λ,p.T',p.λ_loc)#this accounts for mobile susceptibles contracting away from home
    p.λ[:,immobile_age_indices] .= p.λ_loc[:,immobile_age_indices]#This corrects for immobility of susceptibles
    p.λ[ind_mombasa_as,:] .+= p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ[ind_nairobi_as,:] .+= p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)
    return nothing
end

"""
    function rates(out,u,p::CoVParameters_AS,t)

Inplace method for calculating the rate of each of the 8 events per spatial and age group.
    1-> S to E
    2-> E to P
    3-> P to A
    4-> P to M
    5-> P to V
    6-> V to H
    7-> M to R
    8-> A to R
The out vector is in linear index form.
"""
function rates(out,u,p::CoVParameters_AS,t)
    @unpack λ,γ,σ₁,σ₂,δ,τ,μ₁,χ,rel_detection_rate,clear_quarantine,hₐ = p
    calculate_infection_rates!(u,p,t)
    for k = 1:length(out)
        i,a,eventtype = Tuple(index_as_events[k])
        if eventtype ==1
            out[k] = χ[a]*λ[i,a]*u[i,a,1] #Transmission #S->E
        end
        if eventtype ==2
            out[k] = σ₁*u[i,a,2] #E->P
        end
        if eventtype ==3
            out[k] = σ₂*(1-(δ*rel_detection_rate[a]))*u[i,a,3] #P->A
        end
        if eventtype ==4
            out[k] = σ₂*δ*(1-hₐ[a])*rel_detection_rate[a]*u[i,a,3] #P->M
        end
        if eventtype ==5
            out[k] = σ₂*δ*hₐ[a]*rel_detection_rate[a]*u[i,a,3] #P->V
        end
        if eventtype ==6
            out[k] = τ*u[i,a,6] # V->H
        end
        if eventtype ==7
            out[k] = γ*u[i,a,5] # M->R
        end
        if eventtype ==8
            out[k] = γ*u[i,a,4] # A->R
        end
        if eventtype ==14
            out[k] = γ*u[i,a,13] # Q->R
        end
        if eventtype ==15
            out[k] = τ*u[i,a,14] # Qᵥ->H
        end
    end
end

"""
    function change_matrix(dc)

Inplace method for creating a change matrix, that is a matrix encoding the linear map between frequency of events
    occuring and the change this causes in the state of the epidemic. Both the list of events and the list of states
    are encoded in a linear index.
"""
function change_matrix(dc)
    d1,d2 = size(dc)
    for k = 1:d2 #loop over
        i,a,eventtype = Tuple(index_as_events[k]) #
        if eventtype ==1 #S->E
            ind_S = linear_as[i,a,1]
            ind_E = linear_as[i,a,2]
            dc[ind_S,k] = -1
            dc[ind_E,k] = 1
        end
        if eventtype ==2 #E->P
            ind_E = linear_as[i,a,2]
            ind_P = linear_as[i,a,3]
            dc[ind_E,k] = -1
            dc[ind_P,k] = 1
        end
        if eventtype ==3 # P->A
            ind_P = linear_as[i,a,3]
            ind_A = linear_as[i,a,4]
            ind_cumA = linear_as[i,a,9]
            dc[ind_P,k] = -1
            dc[ind_A,k] = 1
            dc[ind_cumA,k] = 1
        end
        if eventtype ==4 # P->M
            ind_P = linear_as[i,a,3]
            ind_M = linear_as[i,a,5]
            ind_cumM = linear_as[i,a,10]
            dc[ind_P,k] = -1
            dc[ind_M,k] = 1
            dc[ind_cumM,k] = 1
        end
        if eventtype ==5 # P->V
            ind_P = linear_as[i,a,3]
            ind_V = linear_as[i,a,6]
            ind_cumV = linear_as[i,a,11]
            dc[ind_P,k] = -1
            dc[ind_V,k] = 1
            dc[ind_cumV,k] = 1
        end
        if eventtype ==6# V->H
            ind_V = linear_as[i,a,6]
            ind_H = linear_as[i,a,7]
            ind_cumVH = linear_as[i,a,12]
            dc[ind_V,k] = -1
            dc[ind_H,k] = 1
            dc[ind_cumVH,k] = 1
        end
        if eventtype ==7 # M->R
            ind_M = linear_as[i,a,5]
            ind_R = linear_as[i,a,8]
            dc[ind_M,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==8 # A->R
            ind_A = linear_as[i,a,4]
            ind_R = linear_as[i,a,8]
            dc[ind_A,k] = -1
            dc[ind_R,k] = 1
        end
        #Mass screening
        if eventtype ==9# E->Q via screening
            ind_E = linear_as[i,a,2]
            ind_Q = linear_as[i,a,13]
            dc[ind_E,k] = -1
            dc[ind_Q,k] = 1
        end
        if eventtype ==10# P->Q via screening
            ind_P = linear_as[i,a,3]
            ind_Q = linear_as[i,a,13]
            dc[ind_P,k] = -1
            dc[ind_Q,k] = 1
        end
        if eventtype ==11# A->Q via screening
            ind_A = linear_as[i,a,4]
            ind_Q = linear_as[i,a,13]
            ind_cumA = linear_as[i,a,9]
            dc[ind_A,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumA,k] = -1
        end
        if eventtype ==12# M->Q via screening
            ind_M = linear_as[i,a,5]
            ind_Q = linear_as[i,a,13]
            ind_cumM = linear_as[i,a,10]
            dc[ind_M,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumM,k] = -1
        end
        if eventtype ==13# V->Qᵥ via screening
            ind_V = linear_as[i,a,6]
            ind_Qᵥ = linear_as[i,a,14]
            ind_cumV = linear_as[i,a,11]
            dc[ind_V,k] = -1
            dc[ind_Qᵥ,k] = 1
            dc[ind_cumV,k] = -1
        end
        if eventtype ==14# Q->R via screening
            ind_Q = linear_as[i,a,13]
            ind_R = linear_as[i,a,8]
            dc[ind_Q,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==15# Qᵥ->H via screening
            ind_Qᵥ = linear_as[i,a,14]
            ind_H = linear_as[i,a,7]
            ind_cumH = linear_as[i,a,12]
            dc[ind_Qᵥ,k] = -1
            dc[ind_H,k] = 1
            dc[ind_cumH,k] = 1
        end
    end
end

"""
    function PP_drivers(dN::Vector{Int64},rates,p)
This method in-place generates the crude number of each type of event proposed by a Poisson process with
rates calculated by rates(out,u,p::CoVParameters_AS,t).
"""
function PP_drivers(dN::Vector{Int64},rates,p)
    for i = 1:length(dN)
        if rates[i] >= 0.
            dN[i] = rand(Poisson(p.dt*rates[i]))
        else
            dN[i] = 0
        end
    end
end

"""
    function max_change(out,u,p::CoVParameters_AS)
This method in-place modifies the number of each type of event proposed by the Poisson process
    so that non-negativity is respected.
"""
function max_change(out,u,p::CoVParameters_AS)
    @unpack δ,rel_detection_rate,hₐ = p
    for i = 1:n_wa,a = 1:n_a
        ind_trans = linear_as_events[i,a,1]
        ind_EP = linear_as_events[i,a,2]
        ind_PA = linear_as_events[i,a,3]
        ind_PM = linear_as_events[i,a,4]
        ind_PV = linear_as_events[i,a,5]
        ind_VH = linear_as_events[i,a,6]
        ind_MR = linear_as_events[i,a,7]
        ind_AR = linear_as_events[i,a,8]
        ind_EQ = linear_as_events[i,a,9]    #MS
        ind_PQ = linear_as_events[i,a,10]
        ind_AQ = linear_as_events[i,a,11]
        ind_MQ = linear_as_events[i,a,12]
        ind_VQᵥ = linear_as_events[i,a,13]
        ind_QR = linear_as_events[i,a,14]
        ind_QᵥH = linear_as_events[i,a,15]

        out[ind_trans] = min(out[ind_trans],u[i,a,1])
        #out[ind_EP] = min(out[ind_EP],u[i,a,2])
        #out[ind_VH] = min(out[ind_VH],u[i,a,6])
        #out[ind_MR] = min(out[ind_MR],u[i,a,5])
        #out[ind_AR] = min(out[ind_AR],u[i,a,4])
        out[ind_EQ] = min(out[ind_EQ],u[i,a,2])     #MS, priority to quarantine events
        out[ind_PQ] = min(out[ind_PQ],u[i,a,3])     #MS
        out[ind_AQ] = min(out[ind_AQ],u[i,a,4])     #MS
        out[ind_MQ] = min(out[ind_MQ],u[i,a,5])     #MS
        #out[ind_VQᵥ]= min(out[ind_VQᵥ],u[i,a,6])    #MS, will give priority to VH over VQᵥ
        out[ind_QR] = min(out[ind_QR],u[i,a,13])    #MS
        out[ind_QᵥH]= min(out[ind_QᵥH],u[i,a,14])   #MS

        if out[ind_EP] + out[ind_EQ] > u[i,a,2] #MS
            out[ind_EP] = u[i,a,2] - out[ind_EQ]
        end
        #splitting events using multinomial sampling #MS, priority to event PQ
        if out[ind_PA] + out[ind_PM] + out[ind_PV] + out[ind_PQ] > u[i,a,3] #More transitions than actual P so all P individuals transition
            rel_rate_each_event = [1-(δ*rel_detection_rate[a]),δ*(1-hₐ[a])*rel_detection_rate[a],δ*hₐ[a]*rel_detection_rate[a]]
            D = rand(Multinomial(u[i, a, 3] - out[ind_PQ], LinearAlgebra.normalize!(rel_rate_each_event,1) )) #Draw the states of the P individuals after transition
            out[ind_PA] = D[1]
            out[ind_PM] = D[2]
            out[ind_PV] = D[3]
        end
        if out[ind_AR] + out[ind_AQ] > u[i,a,4]     #MS, priority to event AQ
            out[ind_AR] = u[i,a,4] - out[ind_AQ]
        end
        if out[ind_MR] + out[ind_MQ] > u[i,a,5]     #MS, priority to event MQ
            out[ind_MR] = u[i,a,5] - out[ind_MQ]
        end
        if out[ind_VH] + out[ind_VQᵥ] > u[i,a,6]     #MS, priority to event VH
            out[ind_VH] = min(out[ind_VH],u[i,a,6])
            out[ind_VQᵥ] = u[i,a,6] - out[ind_VH]
        end

    end
end

"""
    nonneg_tauleap(du,u,p::CoVParameters_AS,t)

This performs one in-place simulation of a time step
"""
function nonneg_tauleap(du,u,p::CoVParameters_AS,t)
    @unpack dc,dN,poi_rates,du_linear = p
    rates(poi_rates,u,p,t) #calculate rates of underlying Poisson processes
    PP_drivers(dN,poi_rates,p)#Generate Poisson rvs with rates scaled by time step dt
    #max_change(dN,u,p)#Cap the size of the Poisson rvs to maintain non-negativity
    mass_screening(dN,u,p,t)    #MS
    max_change(dN,u,p)          #MS
    mul!(du_linear,dc,dN)#Calculates the effect on the state in the inplace du vector
    du .= reshape(du_linear,n_wa,n_a,n_s)
    du .+= u #Calculates how the state should change
end

using SimpleRandom
"""
    function mass_screening(dN,u,p::CoVParameters_AS,t)

This applies the mass screening intervention
"""
function mass_screening(dN,u,p::CoVParameters_AS,t)
    @unpack MS_strategy,MS_probaₐ,dt,MS_nb_testsₐ,MS_probaₐₛ= p
    wa=4    #tests in Nairobi only
    MS_nb_testsᵢ=Int(floor(MS_strategy[Int(ceil(t+1))]*dt))     #puttin all tests in Nairobi
    if MS_nb_testsᵢ!=0
        for a=1:n_a     MS_probaₐ[a]=(sum(u[wa,a,1:6])+u[wa,a,8])/(sum(u[wa,:,1:6])+sum(u[wa,:,8]));    end
        MS_nb_testsₐ = rand(Multinomial(MS_nb_testsᵢ, MS_probaₐ))#LinearAlgebra.normalize!(MS_probaₐ,1) ))

        #=if t==140
            println("\nMS_probaₐ=",MS_probaₐ)
            println("\nMS_nb_testsₐ=",MS_nb_testsₐ)
        end=#
        for a=1:n_a
            for s∈[1,2,3,4,5,6,8]   MS_probaₐₛ[a,s]=u[wa,a,s]/(sum(u[wa,a,1:6])+u[wa,a,8]); end  #we can test s∈{S,E,P,A,M,V,R} and !{Q,Qᵥ,H}
            MS_nb_testsₛ = rand(Multinomial(MS_nb_testsₐ[a], MS_probaₐₛ[a,:]))
            dN[linear_as_events[wa, a, 9]] = MS_nb_testsₛ[2] #s=E
            dN[linear_as_events[wa, a, 10]]= MS_nb_testsₛ[3] #s=P
            dN[linear_as_events[wa, a, 11]]= MS_nb_testsₛ[4] #s=A
            dN[linear_as_events[wa, a, 12]]= MS_nb_testsₛ[5] #s=M
            dN[linear_as_events[wa, a, 13]]= MS_nb_testsₛ[6] #s=V
            #=if t==140
                println("\n\n\nAge=",a,"     proba_per_s=",proba_per_s)
                println("\nnb_tests_per_s=",nb_tests_per_s)
                s=0; println("\ndN=");for e=9:11     s+=dN[linear_as_events[wa, a, e]];print("    dN[4,",a,",",e,"]=",dN[linear_as_events[wa, a, e]])     end
                println("\n--> Age=",a,"   total sN[9:11]=",s)
            end=#
        end
    end
end

"""
    function ode_model(du,u,p::CoVParameters_AS,t)

This is the vector field of the related KenyaCoV ODE model
"""
function ode_model(du,u,p::CoVParameters_AS,t)
    @unpack λ,γ,σ₁,σ₂,δ,τ,μ₁,χ,rel_detection_rate,hₐ = p
    calculate_infection_rates!(u,p,t)
    for i = 1:n_wa,a = 1:n_a
        du[i,a,1] = (-1)*χ[a]*λ[i,a]*u[i,a,1]
        du[i,a,2] = χ[a]*λ[i,a]*u[i,a,1] - σ₁*u[i,a,2]
        du[i,a,3] = σ₁*u[i,a,2] - σ₂*u[i,a,3]
        du[i,a,4] = σ₂*(1-(δ*rel_detection_rate[a]))*u[i,a,3] - γ*u[i,a,4]
        du[i,a,5] = σ₂*δ*(1-hₐ[a])*rel_detection_rate[a]*u[i,a,3] - γ*u[i,a,5]
        du[i,a,6] = σ₂*δ*hₐ[a]*rel_detection_rate[a]*u[i,a,3] - τ*u[i,a,6]
        du[i,a,7] = τ*u[i,a,6]
        du[i,a,8] = γ*u[i,a,5] + γ*u[i,a,4]
        du[i,a,9] = σ₂*(1-(δ*rel_detection_rate[a]))*u[i,a,3]
        du[i,a,10] = σ₂*δ*(1-hₐ[a])*rel_detection_rate[a]*u[i,a,3]
        du[i,a,11] = σ₂*δ*hₐ[a]*rel_detection_rate[a]*u[i,a,3]
        du[i,a,12] = τ*u[i,a,6] #Same as being hosp. for this version
     end
end
