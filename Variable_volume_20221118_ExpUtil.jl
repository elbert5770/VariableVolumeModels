#The model is inspired by the following work: 
# http://www.physics.uoi.gr/assimakopoulos/files/scientific_papers/39_Health_Physics_61_(1991)_245-253.pdf
# However, the model in the paper does not consistently account for flow of liquid from body to the udder to produce
# milk.  This is important because it is unlikely that the cesium wouldn't be secreted into the milk - it is an ion.
# Thus, there needs to be some mass balance where the flow that fills the udder carries cesium with it, as
# is done below (b1).  b2 represents milking - it is either on or off.  Secretion continues even during milking,
# which is biological.  The exchange represented by μ2 and λ2 are in the original model, but have values that
# are too small to affect the results.  They are included for completeness.
using QuadGK
using GLMakie
using DifferentialEquations
using ExponentialUtilities

function variable_volume!(dy,y,p,t) 
    # Equations for ODE solver

    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    # V1 volume of body
    # a volume of milk at time = 0
    # b1 rate of fluid flow from body to udder
    # b2 rate of milk flow out of udder
    # ξ1 rate of cesium clearance volume/time
    # μ2 rate of receptor or channel-mediated transfer of cesium from udder to body volume/time
    # λ2 rate of receptor or channel-mediated transfer of cesium from body to udder volume/time
    # γP rate of absorption of cesium from contaminated grain through intestines Bq/volume/day
    # Note that fluid intake and V1 are constant, so more milk production would be balanced by 
    # less production of urine or feces

    dy[1] = (γP + (μ2*y[3] -(λ2+ξ1+b1)*y[1]))/V1 # Body Cesium concentration Bq/volume/day
    dy[2] = -μ2*y[3] + λ2*y[1] + b1*y[1] +b2*y[3] # Milk Cesium amount Bq/day
    dy[3] = 1/y[4]*dy[2] - y[2]/y[4]^2*(b1+b2) # Milk Cesium concentration Bq/volume/day
    dy[4] = b1+b2   # Rate of change of udder volume
end


function main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)
    
    @show fast_milking
    @show start_with_drain
    @show short_simulation
    #### Constants and settings
    hour_per_day = 24.0 # rate constants given in units of vol/day
    V1 = 33.0 # L Body volume
    a1 = 0.09 # L
    b1 = 2.2/hour_per_day # L/h

    if fast_milking[] == true
        milking_time = 0.1 # h
    else
        milking_time = 5. # h
    end

    milking_period = 12.0 # h

    a2 = a1 + b1*(milking_period-milking_time) # L
    a = a1 # L
    b2 = -b1*(milking_period)/milking_time # L
    ξ1 = 8.6/hour_per_day # L/h
    μ2 = 1e-4/hour_per_day # L/h
    λ2 = 0.1/hour_per_day # L/h

    γP = 351.0/hour_per_day # radioactive cesium absorption, Bq/h

    P = [γP; 0] # cesium inputs, Bq/h



    x0 = [0.0;0.0] # concentrations of cesium in compartments, Bq/vol
    y0 = [0.0;0.0] # amount of cesium in compartments, Bq
    V0 = [V1 0;0 a2] # initial volumes of compartments

    delt_drain = 0.01 # time stepping during draining of the udder (milking)
    delt_fill = 0.1 # time stepping not during milking

    tmax_fill = milking_period-milking_time # max time between milkings for time stepping
    tmax_drain = milking_time # max time during milking

    if start_with_drain[] == true # can start feeding cesium at the start of milking
        x0DE = [0;0;0;a2] # initial conditions for ode
    else # or right after milking
        x0DE = [0;0;0;a1]
    end

    # Set up ODEProblem
    tspan = (0.0, tmax_drain) 
    p = [V1,a,b1,b2,ξ1,μ2,λ2,γP]
    prob = ODEProblem(variable_volume!, x0DE, tspan, p)

    # Will push solutions into arrays
    time = Vector{Float64}([])
    valuex = []
    valuey = []
    volume = Vector{Float64}([])
    timeDE =  Vector{Float64}([])
    valueDE = []
    ####

    #### Calculation of exact solutions and numerical solution
    # Short simulations allow finer  small differences between solutions to be observed
    if short_simulation[] == true
        number_of_days = 1
    else
        number_of_days = 14
    end


    for j = 1:2*number_of_days
        if j > 14
            γP = 0.0 # Cesium contaminated grain is fed for 7 days or 14 milking cycles
            P = [γP,0.0]   
        end
        for k = 0:1
            if k != start_with_drain[] # Settings during milking
                b2 = -b1*(milking_period)/milking_time          
                a = a2  
                V0 = [V1 0;0 a2]
                solution_method = solution_method_drain[]
                tmax = tmax_drain
                delt = delt_drain
                current_time = milking_period*(j-1)+(milking_period - milking_time)*(1-start_with_drain[])
            else # Settings when not milking
                b2 = 0.0
                a = a1
                V0 = [V1 0;0 a1]
                solution_method = solution_method_fill[]
                tmax = tmax_fill
                delt = delt_fill
                current_time = milking_period*(j-1)+milking_time*start_with_drain[]
            end
            p = [V1,a,b1,b2,ξ1,μ2,λ2,γP] # Reset parameters
            
            R = calcR(p) # Reaction rates depend on parameters but not time
            dVdt = calcdVdt(p) # dV/dt depends on paramters but not time
            t = 0.0
            G0 = exponential!(calcintinvVdVdt(t,p) - calcintinvVR(t,p,R)) # Integrating factor at time = 0
        
            for t = 0.0:delt:tmax
                V = calcV(t,p) # Volume of compartments at time = t
                invV = calcinvV(t,p) # Inverse of volume

                ##### Solution calculated from dVx/dt = P + Rx with integrating factor
                if solution_method == 1
                    invG = exponential!(calcintinvVR(t,p,R) - calcintinvVdVdt(t,p))
                    integral, err = quadgk(x -> calcouterint(x,p,invV,R), 0.0, t)     
                    x = invG*integral*P + invG*G0*x0
                    y = V*x
                    #@show t, calcintinvVR(t,p,R), calcintinvVdVdt(t,p), invG, invV, calcouterint(t,p,invV,R)
                end
                #####

                ##### Solution calculated from dy/dt = P + RV^(-1)y with integrating factor
                if solution_method == 2
                    expRintinvV = exponential!(calcRintinvV(t,p,R))
                    integral, err = quadgk(x -> calcexpnegRintinvV(x,p,R), 0.0, t)
                    x = inv(V)*expRintinvV*(integral*P + V0*x0 )
                    y = expRintinvV*(integral*P + y0 )
                    # @show t, calcRintinvV(t,p,R) expRintinvV, calcexpnegRintinvV(t,p,R),integral,inv(V),integral*P 
                end
                #####

                ##### Solution calculated from dVx/dt = P + Rx as separable equation
                if solution_method == 3
                    invG = exponential!((R-dVdt)*calcintinvV(t,p))
                    x = inv(R-dVdt)*(invG*(P+(R-dVdt)*x0) - P)
                    y = V*x
                end
                #####

                push!(time,Float64(t+current_time))
                push!(valuex,x)
                push!(valuey,y)
                push!(volume,Float64(V[2,2]))
            end
            # Set next initial values to last values
            x0 = valuex[end]
            y0 = valuey[end]
            
            # Remake ODE problem and solve
            tspan = (0.0, tmax)
            sol = solve(remake(prob,u0=x0DE,p=p,tspan=tspan), Rosenbrock23(),dtmax = 0.001,saveat=collect(0.0:delt:tmax))
        
            reduce(append!,(timeDE,sol.t .+ current_time))

            push!(valueDE,sol.u)
            x0DE = sol.u[end] # Set next initial values to last value
        end



    end

    # Plotting solutions

    valueplotx = reduce(vcat,transpose.(valuex))
    valueploty = reduce(vcat,transpose.(valuey))
    # timeDE = reduce(hcat,transpose.(timeDE))
    valueplotDE = reduce(vcat,transpose.(reduce(vcat,valueDE)))



    empty!(ax)
    empty!(ax2)
    sc1 = scatter!(ax,time, valueplotx[:,2],color=(:blue,0.5),markersize=3)
    sc2 =scatter!(ax,timeDE, valueplotDE[:,3],color=(:red,0.5),markersize=3)
    sc3 = scatter!(ax,time,volume,color=:orange,markersize=3)
    sa1 = scatter!(ax2,time, valueploty[:,2],color=(:blue,0.5),markersize=3)
    sa2 = scatter!(ax2,timeDE, valueplotDE[:,2],color=(:red,0.5),markersize=3)
    sa3 = scatter!(ax2,time,volume,color=:orange,markersize=3)
    return nothing
end


function create_figure()
    fig = Figure(resolution = (1400, 800))

    # Toggles for boolean variables
    toggle_choices = ["Short simulation","Start with drain","Fast milking"]
    toggles = [Toggle(fig, active=false);
            Toggle(fig,active=true);
            Toggle(fig,active=false)]
    labels = [Label(fig, l) for (t,l) in zip(toggles,toggle_choices)]
    fig[1,3] = grid!(hcat(toggles,labels),tellheight=false)


    # Menus for exact solution methods
    solution_method_drains = [1,2,3]
    solution_method_fills = [1,2,3]
    menu1 = Menu(fig, options = zip(["Conc. Integ. Factor", "Mass Integ. Factor", "Conc. Separable"],solution_method_drains), default = "Conc. Integ. Factor")
    menu2 = Menu(fig, options = zip(["Conc. Integ. Factor", "Mass Integ. Factor", "Conc. Separable"],solution_method_fills), default = "Mass Integ. Factor")
    fig[6, 3] = vgrid!(
    Label(fig, "Drain", width = nothing),
    menu1,
    Label(fig, "Fill", width = nothing),
    menu2;
    tellheight = false, width = 200)

    # Create axes for two plots
    ax = Axis(fig[:,1],xlabel="Time (h)",ylabel="Concentration or volume",title="Concentration of cesium in milk")
    ax2 = Axis(fig[:,2],xlabel="Time (h)",ylabel="Mass or volume",title="'Mass' (i.e. amount) of cesium in milk")
    
    # Create legend
    elem_1 = [MarkerElement(color = :blue, marker=:CIRCLE, markersize = 12)]
    elem_2 = [MarkerElement(color = :red, marker=:CIRCLE, markersize = 12)]
    elem_3 = [MarkerElement(color = :orange, marker=:CIRCLE, markersize = 12)]
    axislegend(ax, [elem_1, elem_2, elem_3], ["Matrix Exponential", "Numerical soln.", "Volume"], position = :lt)
    axislegend(ax2, [elem_1, elem_2,elem_3], ["Matrix Exponential", "Numerical soln.", "Volume"], position = :lt)

    # Show diff eq forms
    Label(fig[2,3],"Diff. Eq. w/concentration")
    Label(fig[3,3],L"V\frac{d\mathbf{C}}{dt}+C\frac{d\mathbf{V}}{dt}=\mathbf{P}+\mathbf{RC}")
    Label(fig[4,3],"Diff. Eq. w/mass")
    Label(fig[5,3],L"\frac{d\mathbf{M}}{dt}=\mathbf{P}+\mathbf{RV}^{-1}\mathbf{M}")
   
    pfig = fig,toggles,labels,menu1,menu2,ax,ax2
    
    return pfig
end

### Functions to calculate terms in solutions

function calcRintinvV(t,p,R)
    # Calculate R*Integrate[inv(V)]
     intinvV = calcintinvV(t,p)
     RintinvV = R*intinvV
     return RintinvV
 end
 
 function calcexpnegRintinvV(t,p,R)
     # Calculate exponential!(-R*Integrate[inv(V)])
     return exponential!(-calcRintinvV(t,p,R))
 end
 
 function calcintinvVR(t,p,R)
    # Calculate Integrate[inv(V)]*R
     intinvV = calcintinvV(t,p)
     intinvVR = intinvV*R
     return intinvVR
 end
 
 function calcouterint(t,p,invV,R)
     # Calculate G*inv(V) for solution method 0
     G = exponential!(calcintinvVdVdt(t,p) - calcintinvVR(t,p,R))
     GinvV = G*invV
     return GinvV
 end

function calcR(p)
    # R = reaction rates (including flows btwn compartments), time invariant
    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    R = [-(λ2+ξ1+b1) μ2;(λ2+b1) -μ2+b2]
    return R
end

function calcV(t,p)
    # V, volumes of compartments, time varying
    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    V = [V1 0;0 a+(b1+b2)*t]
    return V    
end

function calcinvV(t,p)
    # inv(V)
    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    invV = [1/V1 0;0 1/(a+(b1+b2)*t)]
    return invV 
end

function calcdVdt(p)
    # dV/dt
    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    dVdt = [0.0 0.0;0.0 b1+b2]
    return dVdt
end

function calcintinvV(t,p)
    # Integrate[inv(V),t]
    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    intinvV = [t/V1 0;0 log(1+(b1+b2)*t/a)/(b1+b2)]
    return intinvV
end

function calcintinvVdVdt(t,p)
    # Integrate[inv(V)*dV/dt,t]
    V1,a,b1,b2,ξ1,μ2,λ2,γP = p
    intinvVdVdt = [0.0 0.0;0.0 log(1+(b1+b2)*t/a)]
    return intinvVdVdt
end

# function mat_exp(A)
#     # Exact calculation of matrix exponentials for a 2x2 matrix, has no noticable effect on solution
#     # and is not used here
#     a = A[1,1]
#     b = A[1,2]
#     c = A[2,1]
#     d = A[2,2]

#     B = 1/2*sqrt((a-d)^2 + 4*b*c)
#     if B == 0.0
#         RV = [1.0 0.0;0.0 1.0]
#     else
#         RV = exponential!((a+d)/2).*[cosh(B)+(a-d)/(2*B)*sinh(B) b/B*sinh(B);c/B*sinh(B) cosh(B)-(a-d)/(2*B)*sinh(B)]
#     end
#     return RV
# end

#### Form of the exact solution affects accuracy
#### Exact solutions assume that R, P and dV/dt are not functions of t
#### 1 = Solution calculated from dVx/dt = P + Rx with integrating factor
####    Almost identical to numerical soln during drain period, somewhat off for fill period
#### 2 = Solution calculated from dy/dt = P + RV^(-1)y with integrating factor
####    Almost identical to numerical soln for fill period, not close for drain period
#### 3 = Solution calculated from dVx/dt = P + Rx as separable equation
####    Worse than 1 for fill, better than 2 for drain

# Default values
# solution_method_drain = 1
# solution_method_fill = 2
# start_with_drain = true
# short_simulation = false
# fast_milking = false

pfig = create_figure()
fig,toggles,labels,menu1,menu2,ax,ax2 =pfig
display(fig)
short_simulation=toggles[1].active
start_with_drain=toggles[2].active
fast_milking=toggles[3].active
solution_method_drain = menu1.selection
solution_method_fill = menu2.selection


main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)

on(menu1.selection) do s
    solution_method_drain = s
    main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)
end
on(menu2.selection) do s
    solution_method_fill = s
    main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)
end

on(toggles[1].active) do val
    short_simulation=toggles[1].active 
    main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)
end
on(toggles[2].active) do val
    start_with_drain=toggles[2].active 
    main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)
end
on(toggles[3].active) do val
    fast_milking=toggles[3].active 
    main(ax,ax2,solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)
end
