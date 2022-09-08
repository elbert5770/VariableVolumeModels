using QuadGK
using Plots
using DifferentialEquations
gr()

function calcRintinvV(t,p,R)
   # Calculate R*Integrate[inv(V)]
    intinvV = calcintinvV(t,p)
    RintinvV = R*intinvV
    return RintinvV
end

function calcexpnegRintinvV(t,p,R)
    # Calculate exp(-R*Integrate[inv(V)])
    return exp(-calcRintinvV(t,p,R))
end

function mat_exp(A)
    # Exact calculation of matrix exponentials for a 2x2 matrix, has no noticable effect on solution
    a = A[1,1]
    b = A[1,2]
    c = A[2,1]
    d = A[2,2]

    B = 1/2*sqrt((a-d)^2 + 4*b*c)
    if B == 0.0
        RV = [1.0 0.0;0.0 1.0]
    else
        RV = exp((a+d)/2).*[cosh(B)+(a-d)/(2*B)*sinh(B) b/B*sinh(B);c/B*sinh(B) cosh(B)-(a-d)/(2*B)*sinh(B)]
    end
    return RV
end

function calcintinvVR(t,p,R)
   # Calculate Integrate[inv(V)]*R
    intinvV = calcintinvV(t,p)
    intinvVR = intinvV*R
    return intinvVR
end

function calcouterint(t,p,invV,R)
    # Calculate G*inv(V) for solution method 0
    G = exp(calcintinvVdVdt(t,p) - calcintinvVR(t,p,R))
    GinvV = G*invV
    return GinvV
end


function variable_volume!(dy,y, p, t) 
    # Equations for ODE solver

    V1,a,b1,b2,ξ1,μ2,λ2,γP = p

    dy[1] = (γP + (μ2*y[3] -(λ2+ξ1+b1)*y[1]))/V1
    dy[2] = (-μ2*y[3] + λ2*y[1]) + b1*y[1] +b2*y[3]
    dy[3] = 1/y[4]*dy[2] - y[2]/y[4]^2*(b1+b2)
    dy[4] = b1+b2   
end

function main(solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)

    #### Constants and settings
    hour_per_day = 24.0 # rate constants given in units of vol/day
    V1 = 33.0 # L
    a1 = 0.09 # L
    b1 = 2.2/hour_per_day # L/h

    if fast_milking == true
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
    
    if start_with_drain == true # can start feeding cesium at the start of milking
        x0DE = [0;0;0;a2] # initial conditions for ode
    else # or right after milking
        x0DE = [0;0;0;a1]
    end

    # Set up ODEProblem
    tspan = (0.0, tmax_drain) 
    p = [V1,a,b1,b2,ξ1,μ2,λ2,γP]
    prob = ODEProblem(variable_volume!, x0DE, tspan, p)

    # Will push solutions into arrays
    time = []
    valuex = []
    valuey = []
    volume = []
    timeDE = []
    valueDE = []
    ####

    #### Calculation of exact solutions and numerical solution
    # Short simulations allow finer  small differences between solutions to be observed
    if short_simulation == true
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
            if k != start_with_drain # Settings during milking
                b2 = -b1*(milking_period)/milking_time          
                a = a2  
                V0 = [V1 0;0 a2]
                solution_method = solution_method_drain
                tmax = tmax_drain
                delt = delt_drain
                current_time = milking_period*(j-1)+(milking_period - milking_time)*(1-start_with_drain)
            else # Settings when not milking
                b2 = 0.0
                a = a1
                V0 = [V1 0;0 a1]
                solution_method = solution_method_fill
                tmax = tmax_fill
                delt = delt_fill
                current_time = milking_period*(j-1)+milking_time*start_with_drain
            end
            p = [V1,a,b1,b2,ξ1,μ2,λ2,γP] # Reset parameters

            R = calcR(p) # Reaction rates depend on parameters but not time
            dVdt = calcdVdt(p) # dV/dt depends on paramters but not time
            t = 0.0
            G0 = exp(calcintinvVdVdt(t,p) - calcintinvVR(t,p,R)) # Integrating factor at time = 0
        
            for t = 0.0:delt:tmax
                V = calcV(t,p) # Volume of compartments at time = t
                invV = calcinvV(t,p) # Inverse of volume

                ##### Solution calculated from dVx/dt = P + Rx with integrating factor
                if solution_method == 0
                    invG = exp(calcintinvVR(t,p,R) - calcintinvVdVdt(t,p))
                    integral, err = quadgk(x -> calcouterint(x,p,invV,R), 0.0, t)     
                    x = invG*integral*P + invG*G0*x0
                    y = V*x
                   # @show t, calcintinvVR(t,p,R), calcintinvVdVdt(t,p), invG, invV, calcouterint(t,p,invV,R)
                end
                #####

                ##### Solution calculated from dy/dt = P + RV^(-1)y with integrating factor
                if solution_method == 1
                    expRintinvV = exp(calcRintinvV(t,p,R))
                    integral, err = quadgk(x -> calcexpnegRintinvV(x,p,R), 0.0, t)
                    x = inv(V)*expRintinvV*(integral*P + V0*x0 )
                    y = expRintinvV*(integral*P + y0 )
                   # @show t, calcRintinvV(t,p,R) expRintinvV, calcexpnegRintinvV(t,p,R),integral,inv(V),integral*P 
                end
                #####

                ##### Solution calculated from dVx/dt = P + Rx as separable equation
                if solution_method == 2
                    invG = exp((R-dVdt)*calcintinvV(t,p))
                    x = inv(R-dVdt)*(invG*(P+(R-dVdt)*x0) - P)
                    y = V*x
                end
                #####

                push!(time,t+current_time)
                push!(valuex,x)
                push!(valuey,y)
                push!(volume,V[2,2])
            end
            # Set next initial values to last values
            x0 = valuex[end]
            y0 = valuey[end]
            
            # Remake ODE problem and solve
            tspan = (0.0, tmax)
            sol = solve(remake(prob,u0=x0DE,p=p,tspan=tspan), Rosenbrock23(),dtmax = 0.001,saveat=collect(0.0:delt:tmax))
        
            push!(timeDE,sol.t .+ current_time)
            push!(valueDE,sol.u)
            x0DE = sol.u[end] # Set next initial values to last value
        end



    end

    # Plotting solutions

    valueplotx = reduce(vcat,transpose.(valuex))
    valueploty = reduce(vcat,transpose.(valuey))
    timeDE = reduce(hcat,transpose.(timeDE))
    valueplotDE = reduce(vcat,transpose.(reduce(vcat,valueDE)))


    #    p1 = plot(time,valueplotx,legend=false,ylim=(0,31),label=["Cesium blood conc" "Cesium milk conc"])
    #    plot!(p1,time,valueploty)
    #    plot!(p1,time,volume)

    #   p2 = plot(timeDE',valueplotDE,legend=false,ylim=(0,31),label=["Cesium blood conc" "Cesium milk mass" "Cesium milk conc" "udder volume"])
    
    p3 = scatter(time,valueplotx[:,2] ,markersize=1,markerstrokewidth=0,markercolor=:blue,label="Conc (exact)",legendfontsize=5)
    scatter!(p3,timeDE', valueplotDE[:,3] ,markersize=1,markerstrokewidth=0,markercolor=:red,label="Conc (ODEsolve)")
    scatter!(p3,time,volume,markersize=1,markerstrokewidth=0,label="Volume")
    #scatter!(p3,timeDE',valueplotDE[:,4],markersize=1,markerstrokewidth=0)
    # scatter!(p3,time,valueplotx[:,1],markersize=1,markercolor=:yellow,markerstrokewidth=0)
    # scatter!(p3,timeDE',valueplotDE[:,1],markersize=1,markercolor=:green,markerstrokewidth=0)
   xlabel!(p3,"Time (h)")
   ylabel!(p3,"Concentration or volume")
    p4 = scatter(time,valueploty[:,2],markersize=1,markerstrokewidth=0,markercolor=:blue,label="Exact",legendfontsize=6)
    scatter!(p4,timeDE', valueplotDE[:,2],markersize=1,markercolor=:red,markerstrokewidth=0,label="ODEsolve")
    xlabel!(p4,"Time (h)")
    ylabel!(p4,"Amount")
    display(plot(p3,p4,layout=(1,2)))
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

#### Form of the exact solution affects accuracy
#### Exact solutions assume that R, P and dV/dt are not functions of t
#### 0 = Solution calculated from dVx/dt = P + Rx with integrating factor
#### Almost identical to numerical soln for drain period, somewhat off for fill period
#### 1 = Solution calculated from dy/dt = P + RV^(-1)y with integrating factor
#### Almost identical to numerical soln for fill period, not close for drain period
#### 2 = Solution calculated from dVx/dt = P + Rx as separable equation
#### Worse than 0 for fill, better than 1 for drain
solution_method_drain = 0
solution_method_fill = 1
start_with_drain = true
short_simulation = false
fast_milking = false
main(solution_method_drain,solution_method_fill,start_with_drain,short_simulation,fast_milking)

