using LinearAlgebra
using CSV
using DataFrames
using PyPlot

include("getdata.jl")

moonmass=7.3477e22
M0=5.972e24

function simulate(A, B, C)
    #= simulates planet deformation given coefficients A, B, and C, where
    A = kappa0*q0^2/(G0*M0^2)
    B = ka0*L0^3/(G0*M0^2)
    C = kc0*L0^3/(G0*M0^2)
    =#
    N = 10
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 1
    G = 1
    vi = 0.3
    ui = Vector{Float64}()
    for i = 1:N
        push!(ui, cos(2 * pi * i / N))
        push!(ui, sin(2 * pi * i / N))
    end
    for i = 1:N
        push!(ui, vi * -sin(2 * pi * i / N))
        push!(ui, vi * cos(2 * pi * i / N))
    end
    mi = moonmass / (M0 * N)
    qi = 1 / N
    mlist = fill(mi, N)
    qlist = fill(qi, N)
    ka = 1
    kc = 1
    d = 1
    refdist = 0
    #this is at arbitrary time
    #qlist is list with length N that includes list of each charges
    # mlist is list with length N that includes mass of each charges
    #kappa=8.99 x 10^9 N·m^2/C^2
    #maybe implement varying ka
    #G=6.67 x 10^–11 m^3 kg^–1 s^–2
    #ka is the spring constant between each point mass
    #kc is the spring constant between the point mass and the center
    #u is stored as [x1;y1;x2;y2;.......xn;yn;dx1dt;dy1dt;dx2dt;dy2dt;,,,,,,,dxndt;dyndt] meaning this u will be length of 4N
    #For now we are only working on implementing this coloumb forces
    # let's find how many body we have

    u = zeros(Int64(4 * N))
    # I created this ulist array so that we can store u from each time step
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    #let's say we already have Nbodydynamics code built
    for t = 1:Nt
        #Let's start with finding those k1k2k3k4
        #dudt(y)=Nbodydynamicsfunction
        # we don't have time dependence our function has format of Nbodydynamics(u,t,kappa,qlist,mlist)
        # technically we don't even have explicit time dependence meaning it depends on the time in the sense of
        #dydt=Nbodydynamics(y(t))
        k1 = Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k2 = Nbodydynamicssimulator(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k3 = Nbodydynamicssimulator(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k4 = Nbodydynamicssimulator(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)

        # our goal is to end current time step by finding U at next time step
        # println(k1+2*k2+2*k3+k4)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

    println(ulist[1:4, 1]) # initial position
    println(ulist[1:4, Nt]) # final position

    tempx = zeros(Int64(N))
    tempy = zeros(Int64(N))
    # Graphing as pos of time
    for t = 1:50:Nt
        for i = 1:N
            tempx[i] = ulist[2i-1, t]
            tempy[i] = ulist[2i, t]
        end
        clf()
        scatter(tempx, tempy)
        xlim(-10, 10)
        ylim(-10, 10)
        if t == 1
            pause(5)
        else
            pause(0.1)
        end
    end
    dimulist = size(ulist)
    angularv1list = zeros(2, dimulist[2])
    angularv1maglist = zeros(dimulist[2])
    for t = 1:Nt
        pos1 = ulist[1:2, t]
        v1 = ulist[2*N+1:2*N+2, t]
        angularv1 = v1 - dot(pos1, v1) / dot(pos1, pos1) * pos1
        angularv1list[:, t] = angularv1
        angularv1mag = norm(angularv1)
        angularv1maglist[t] = angularv1mag
    end
    # Plot: Angular Velocity 1
    figure()
    plot(1:Nt, angularv1maglist)
    xlabel("t")
    ylabel("Angular Velocity of Point 1")
    title("Angular Velocity vs. Time")
    show()
    # Plot: Trajectory 1
    pos = ulist[1:2, :]  # assuming ulist is a 2D array (e.g., 2 x N)
    figure()
    plot(pos[1, :], pos[2, :])
    xlabel("x")
    ylabel("y")
    title("Trajectory")
    show()
    
    # r is the average radius of the orbit, defined as the average distance between the geometric center of the points and one of the points
    global r = 0
    xiavglist = Float64[]
    yiavglist = Float64[]
    # ulist is stored as [x1(1),y1(1),x2(1),y2(1),...,xn(1),yn(1),...,x1(Nt),y2,(Nt),...,xn(Nt),yn(Nt)]
    for i = 1:N
        xi = ulist[2i-1, :]
        yi = ulist[2i, :]
        for dt = 1:Nt
            sumxi = sum(ulist[1:2:end, dt])
            sumyi = sum(ulist[2:2:end, dt])
            xiavg = sumxi / N
            yiavg = sumyi / N
            if i == 1
                push!(xiavglist, xiavg)
                push!(yiavglist, yiavg)
            end
            center = [xiavg; yiavg]
            tempr = sqrt((center[1] - xi[dt])^2 + (center[2] - yi[dt])^2) # distance between (x1,y1) and center
            global r += tempr
        end
    end
    r = r / (N * Nt)
    println("r=$r")
    # GRAPH: x1-xavg, y1-yavg, as a function of time
    x1 = ulist[1, :]
    y1 = ulist[2, :]
    deltax = x1 .- xiavglist
    deltay = y1 .- yiavglist
    figure()
    plot(1:Nt, deltax, "b-") # blue = x
    plot(1:Nt, deltay, "r-") # red = y
    xlabel("t")
    ylabel("Δ")
    title("Deviation of coordinates from the Average Over Time")
    show()
end

function Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
    #this is at arbitrary time
    #qlist is list with length N that includes list of each charges
    # mlist is list with length N that includes mass of each charges
    #kappa=8.99 x 10^9 N·m^2/C^2
    #maybe implement varying ka
    #G=6.67 x 10^–11 m^3 kg^–1 s^–2
    #ka is the spring constant between each point mass
    #kc is the spring constant between the point mass and the center
    #d is the damping factor
    #refdist is the relaxed length of the spring
    #u is stored as [x1;y1;x2;y2;.......xn;yn;dx1dt;dy1dt;dx2dt;dy2dt;,,,,,,,dxndt;dyndt] meaning this u will be length of 4N
    #For now we are only working on implementing this coloumb forces
    # let's find how many body we have
    N = length(qlist)
    Nacceleration = zeros(2 * N)
    for i = 1:N
        # we are trying to find overall force that results on ith charge
        overallaxi = 0
        overallayi = 0
        xi = u[Int64(2 * (i - 1) + 1)]
        yi = u[Int64(2 * (i))]
        qi = qlist[i]
        mi = mlist[i]
        for j = 1:N
            if i != j
                xj = u[Int64(2 * (j - 1) + 1)]
                yj = u[Int64(2 * (j))]
                #rij is a unit vector in the direction of ij
                rij = [xj - xi; yj - yi]
                rmag = norm(rij)
                rij = rij / rmag
                # how do we find qi and qj? qi=qlist[i] qj=qlist[j]
                # how do we find mi? mi=mlist[i]
                qj = qlist[j]
                mj = mlist[j]
                # electric force
                tempaxielectric = A * -kappa * qi * qj / (rmag^2 * rij[1] * mi)
                tempayielectric = A * -kappa * qi * qj / (rmag^2 * rij[2] * mi)
                overallaxi = overallaxi + tempaxielectric
                overallayi = overallayi + tempayielectric
            end
            # ka spring force
            if i == 1
                xibefore = u[Int64(2 * (N - 1) + 1)]
                yibefore = u[Int64(2 * (N))]
            else
                xibefore = u[Int64(2 * (i - 1))-1]
                yibefore = u[Int64(2 * (i - 1))]
            end
            if i == N
                xinext = u[1]
                yinext = u[2]
            else
                xinext = u[Int64(2 * (i) + 1)]
                yinext = u[Int64(2 * (i) + 2)]
            end
            riinext = [xinext - xi; yinext - yi]
            rmagnext = norm(riinext)
            riinext = riinext / rmagnext
            riibefore = [xibefore - xi; yibefore - yi]
            rmagbefore = norm(riibefore)
            riibefore = riibefore / rmagbefore

            tempaxispring = B * ka * (rmagnext - refdist) * riinext[1] / mi
            tempayispring = B * ka * (rmagnext - refdist) * riinext[2] / mi
            overallaxi = overallaxi + tempaxispring
            overallayi = overallayi + tempayispring
            tempaxispring = B * ka * (rmagbefore - refdist) * riibefore[1] / mi
            tempayispring = B * ka * (rmagbefore - refdist) * riibefore[2] / mi
            overallaxi = overallaxi + tempaxispring
            overallayi = overallayi + tempayispring
        end
        # Calculating the center of the moon
        sumxi = sum(u[1:2:Int64(2*N-1)])
        sumyi = sum(u[2:2:Int64(2*N)])
        xiavg = sumxi / N
        yiavg = sumyi / N

        # kc spring force
        ri11 = [xiavg - xi; yiavg - yi] # center at [xiavg,yiavg]
        ri11mag = norm(ri11)
        ri11 = ri11 / ri11mag

        tempaxispring = C * kc * (ri11mag - refdist) * ri11[1] / mi
        tempayispring = C * kc * (ri11mag - refdist) * ri11[2] / mi
        overallaxi = overallaxi + tempaxispring
        overallayi = overallayi + tempayispring
        # add earth and sun based on data
        m12 = 1 # earth
        m13 = 1.989e30 / M0  # sun
        # Computing relative positions of the earth and the sun according to the center
        ri12 = [earthposrel[round(Int, t)][1] - xi + xiavg; earthposrel[round(Int, t)][2] - yi + yiavg]
        ri13 = [sunposrel[round(Int, t)][1] - xi + xiavg; sunposrel[round(Int, t)][2] - yi + yiavg]
        ri12mag = norm(ri12)
        ri13mag = norm(ri13)
        ri12 = ri12 / ri12mag
        ri13 = ri13 / ri13mag
        tempaxigravity = G * m12 / ri12mag^2 * ri12[1]
        tempayigravity = G * m12 / ri12mag^2 * ri12[2]
        overallaxi = overallaxi + tempaxigravity
        overallayi = overallayi + tempayigravity
        tempaxigravity = G * m13 / ri13mag^2 * ri13[1]
        tempayigravity = G * m13 / ri13mag^2 * ri13[2]
        overallaxi = overallaxi + tempaxigravity
        overallayi = overallayi + tempayigravity
        # damping
        vdotri11 = u[2*N+(2*i-1)] * ri11[1] + u[2*N+(2*i)] * ri11[2]
        overallaxi = overallaxi - d * (vdotri11) * ri11[1]
        overallayi = overallayi - d * (vdotri11) * ri11[2]
        Nacceleration[Int64(2 * (i - 1) + 1)] = overallaxi
        Nacceleration[Int64(2 * (i))] = overallayi
    end
    dudt = [u[Int64(2 * N + 1):end]; Nacceleration]
    #u is stored as [x1;y1;x2;y2;.......xn;yn;dx1dt;dy1dt;dx2dt;dy2dt;,,,,,,,dxndt;dyndt]
    # dudt should be [dx1dt;dy1dt;dx2dt;dy2dt;.......dxndt;dyndt;d2x1dt2;d2y1dt2;.....d2xndt2;d2yndt2]
    return dudt
end

function averageRadius(A, B, C)
    N = 10
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 1
    G = 1
    vi = 0.3
    ui = Vector{Float64}()
    for i = 1:N
        push!(ui, cos(2 * pi * i / N))
        push!(ui, sin(2 * pi * i / N))
    end
    for i = 1:N
        push!(ui, vi * -sin(2 * pi * i / N))
        push!(ui, vi * cos(2 * pi * i / N))
    end
    mi = moonmass / (M0 * N)
    qi = 1 / N
    mlist = fill(mi, N)
    qlist = fill(qi, N)
    ka = 1
    kc = 1
    d = 1
    refdist = 0
    #this is at arbitrary time
    #qlist is list with length N that includes list of each charges
    # mlist is list with length N that includes mass of each charges
    #kappa=8.99 x 10^9 N·m^2/C^2
    #maybe implement varying ka
    #G=6.67 x 10^–11 m^3 kg^–1 s^–2
    #ka is the spring constant between each point mass
    #kc is the spring constant between the point mass and the center
    #u is stored as [x1;y1;x2;y2;.......xn;yn;dx1dt;dy1dt;dx2dt;dy2dt;,,,,,,,dxndt;dyndt] meaning this u will be length of 4N
    #For now we are only working on implementing this coloumb forces
    # let's find how many body we have

    u = zeros(Int64(4 * N))
    # I created this ulist array so that we can store u from each time step
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    #let's say we already have Nbodydynamics code built
    for t = 1:Nt
        #Let's start with finding those k1k2k3k4
        #dudt(y)=Nbodydynamicsfunction
        # we don't have time dependence our function has format of Nbodydynamics(u,t,kappa,qlist,mlist)
        # technically we don't even have explicit time dependence meaning it depends on the time in the sense of
        #dydt=Nbodydynamics(y(t))
        k1 = Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k2 = Nbodydynamicssimulator(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k3 = Nbodydynamicssimulator(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k4 = Nbodydynamicssimulator(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)

        # our goal is to end current time step by finding U at next time step
        # println(k1+2*k2+2*k3+k4)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

    # r is the average radius of the orbit, defined as the average distance between the geometric center of the points and one of the points
    global r = 0
    xiavglist = Float64[]
    yiavglist = Float64[]
    # ulist is stored as [x1(1),y1(1),x2(1),y2(1),...,xn(1),yn(1),...,x1(Nt),y2,(Nt),...,xn(Nt),yn(Nt)]
    for i = 1:N
        xi = ulist[2i-1, :]
        yi = ulist[2i, :]
        for dt = 1:Nt
            sumxi = sum(ulist[1:2:end, dt])
            sumyi = sum(ulist[2:2:end, dt])
            xiavg = sumxi / N
            yiavg = sumyi / N
            if i == 1
                push!(xiavglist, xiavg)
                push!(yiavglist, yiavg)
            end
            center = [xiavg; yiavg]
            tempr = sqrt((center[1] - xi[dt])^2 + (center[2] - yi[dt])^2) # distance between (x1,y1) and center
            global r += tempr
        end
    end
    r = r / (N * Nt)
    return r
end

function averageAngularVelocity(A,B,C)
    N = 10
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 1
    G = 1
    vi = 0.3
    ui = Vector{Float64}()
    for i = 1:N
        push!(ui, cos(2 * pi * i / N))
        push!(ui, sin(2 * pi * i / N))
    end
    for i = 1:N
        push!(ui, vi * -sin(2 * pi * i / N))
        push!(ui, vi * cos(2 * pi * i / N))
    end
    mi = moonmass / (M0 * N)
    qi = 1 / N
    mlist = fill(mi, N)
    qlist = fill(qi, N)
    ka = 1
    kc = 1
    d = 1
    refdist = 0
    #this is at arbitrary time
    #qlist is list with length N that includes list of each charges
    # mlist is list with length N that includes mass of each charges
    #kappa=8.99 x 10^9 N·m^2/C^2
    #maybe implement varying ka
    #G=6.67 x 10^–11 m^3 kg^–1 s^–2
    #ka is the spring constant between each point mass
    #kc is the spring constant between the point mass and the center
    #u is stored as [x1;y1;x2;y2;.......xn;yn;dx1dt;dy1dt;dx2dt;dy2dt;,,,,,,,dxndt;dyndt] meaning this u will be length of 4N
    #For now we are only working on implementing this coloumb forces
    # let's find how many body we have

    u = zeros(Int64(4 * N))
    # I created this ulist array so that we can store u from each time step
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    #let's say we already have Nbodydynamics code built
    for t = 1:Nt
        #Let's start with finding those k1k2k3k4
        #dudt(y)=Nbodydynamicsfunction
        # we don't have time dependence our function has format of Nbodydynamics(u,t,kappa,qlist,mlist)
        # technically we don't even have explicit time dependence meaning it depends on the time in the sense of
        #dydt=Nbodydynamics(y(t))
        k1 = Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k2 = Nbodydynamicssimulator(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k3 = Nbodydynamicssimulator(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)
        k4 = Nbodydynamicssimulator(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist,A,B,C)

        # our goal is to end current time step by finding U at next time step
        # println(k1+2*k2+2*k3+k4)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end
    dimulist = size(ulist)
    angularv1list = zeros(2, dimulist[2])
    angularv1maglist = zeros(dimulist[2])
    for t = 1:Nt
        pos1 = ulist[1:2, t]
        v1 = ulist[2*N+1:2*N+2, t]
        angularv1 = (v1 - dot(pos1, v1) / dot(pos1, pos1) * pos1) / norm(pos1)
        angularv1list[:, t] = angularv1
        angularv1mag = norm(angularv1)
        angularv1maglist[t] = angularv1mag
    end
    return sum(angularv1maglist) / length(angularv1maglist)
end