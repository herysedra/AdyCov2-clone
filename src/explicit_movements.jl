"""
Basic representation of the state of the age structured model with
explicit movements:

u[home_area,current_area,age_group,disease state]

States:
1 -> S
2 -> E
3 -> I_subclinical
4 -> I_diseased
5 -> H(ospitalised)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative Dead

Events for each home + current area pair and age group:

1-> Transmission
2-> Incubation into asymptomatic E->A
3-> Incubation into diseased E->D
4-> Diseased become hospitalised/treated D->H
5-> Hospitalised/treated recover
6-> Diseased recover D->R
7-> Asymptomatics recover A->R
8->Hospitalised/treated->death

"""
function generate_movements(u,p::CoVParameters_AS_EM,t)
    for k in index_as_em
        i,j,a,s = Tuple(k)
        moveprob = 0.
        if a in mobile_age_indices && i!=j
            moveprob = 1 - exp(-p.r[i]*p.dt) #Chance of returning
            p.du_moves[k] = rand(Binomial(u[k],moveprob ))
        end
        if a in mobile_age_indices && i==j
            moveprob = 1 - exp(-p.l[i]*p.dt) #Chance of leaving
            p.du_moves[k] = rand(Binomial(u[k],moveprob ))
        end
end

function calculate_infection_rates!(u,p::CoVParameters_AS_EM,t)
    I_A = @view u[:,:,:,3]
    I_D = @view u[:,:,:,4]
    p.Î .= sum(p.ϵ*I_A .+ I_D,dims=1)[1,:,:]#Local infecteds
    mul!(p.λ,p.β .*(p.Î ./p.N̂),p.M')#Local force of infection due to age-mixing
    p.λ[ind_mombasa_as,:] .+= p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ[ind_nairobi_as,:] .+= p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)
    return nothing
end
