using LinearAlgebra
using PyPlot
include("resultsfunctions.jl")

# GRAPH: A
Avalues = Float64[]
Arvalues = Float64[]
r = 0
for i = -36:-13 # A = 10^-36 to 10^-10
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
        k1 = Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist, 10.0^i, 2.64e-5, 10^-5)
        k2 = Nbodydynamicssimulator(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist, 10.0^i, 2.64e-5, 10^-5)
        k3 = Nbodydynamicssimulator(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist, 10.0^i, 2.64e-5, 10^-5)
        k4 = Nbodydynamicssimulator(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist, 10.0^i, 2.64e-5, 10^-5)

        # our goal is to end current time step by finding U at next time step
        # println(k1+2*k2+2*k3+k4)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

    tempx = zeros(Int64(N))
    tempy = zeros(Int64(N))

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
    println("A = $(10.0^i) | r = $r")
    push!(Avalues, i)
    push!(Arvalues, r)
    r = 0
end
println("=======================")
figure()
plot(Avalues, Arvalues)
xlabel("log A")
ylabel("r")
title("Average Orbit Radius r vs. log A")
show()

# GRAPH: B
Bvalues = Float64[]
Brvalues = Float64[]
r = 0
for i = 10^-5:10^-5:10^-3 # B = 10^-5 to 10^-3
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
        k1 = Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, i, 10^-5)
        k2 = Nbodydynamicssimulator(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, i, 10^-5)
        k3 = Nbodydynamicssimulator(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, i, 10^-5)
        k4 = Nbodydynamicssimulator(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, i, 10^-5)

        # our goal is to end current time step by finding U at next time step
        # println(k1+2*k2+2*k3+k4)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

    tempx = zeros(Int64(N))
    tempy = zeros(Int64(N))

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
    println("B = $(i) | r = $r")
    push!(Bvalues, i)
    push!(Brvalues, r)
    r = 0
end
println("=======================")
figure()
plot(Bvalues, Brvalues)
xlabel("B")
ylabel("r")
title("Average Orbit Radius r vs. B")
show()

# GRAPH: C
Cvalues = Float64[]
Crvalues = Float64[]
r = 0
for i = 10^-5:10^-5:10^-3 # B = 10^-5 to 10^-3
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
        k1 = Nbodydynamicssimulator(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, 2.64e-5, i)
        k2 = Nbodydynamicssimulator(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, 2.64e-5, i)
        k3 = Nbodydynamicssimulator(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, 2.64e-5, i)
        k4 = Nbodydynamicssimulator(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist, 10^-35, 2.64e-5, i)

        # our goal is to end current time step by finding U at next time step
        # println(k1+2*k2+2*k3+k4)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

    tempx = zeros(Int64(N))
    tempy = zeros(Int64(N))

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
    println("C = $(i) | r = $r")
    push!(Cvalues, i)
    push!(Crvalues, r)
    r = 0
end
println("=======================")
figure()
plot(Cvalues, Crvalues)
xlabel("C")
ylabel("r")
title("Average Orbit Radius r vs. C")
show()

# CONTOUR PLOT: r vs. B, C
figure()
a = 1e-35
B = range(1e-7, 1e-3, length=10)
C = range(1e-7, 1e-3, length=10)
r = [averageRadius(a, b, c) for c in C, b in B];  # (x,y,z)=(B,C,r)
cf = contourf(B, C, r, 20,cmap="rainbow")#levels=10,contour_lines=true,20)
colorbar(cf)
title("Plot of r vs. B and C")
xlabel("B")
ylabel("C")
savefig("contourBCr.pdf",dpi=300)
# CONTOUR PLOT: r vs. A, B
figure()
A = range(1e-35, 1e-32, length=10)
B = range(1e-7, 1e-3, length=10)
c = 1e-5
r = [averageRadius(a, b, c) for b in B, a in A];  # (x,y,z)=(A,B,r)
cf = contourf(A, B, r, 20,cmap="rainbow")#levels=10,contour_lines=true,20)
colorbar(cf)
title("Plot of r vs. A and B")
xlabel("A")
ylabel("B")
savefig("contourABr.pdf",dpi=300)
# CONTOUR PLOT: r vs. A, C
figure()
A = range(1e-35, 1e-32, length=10)
b = 2.64e-5
C = range(1e-7, 1e-3, length=10)
r = [averageRadius(a, b, c) for c in C, a in A];  # (x,y,z)=(A,C,r)
cf = contourf(A, C, r, 20,cmap="rainbow")#levels=10,contour_lines=true,20)
colorbar(cf)
title("Plot of r vs. A and C")
xlabel("A")
ylabel("C")
savefig("contourACr.pdf",dpi=300)

# PLOT: omega vs. A
figure()
e = -36:-13
A = 10.0 .^ e
logA = collect(e)
b = 2.64e-5
c = 1e-5
omega = [averageAngularVelocity(a, b, c) for a in A];
cf = plot(logA, omega)
title("Plot of ω vs. log A")
xlabel("log A")
ylabel("ω")
# PLOT: omega vs. B
figure()
a=1.74e-34
B=10^-5:10^-5:10^-3
c=1e-5
omega = [averageAngularVelocity(a, b, c) for b in B];
cf = plot(B, omega)
title("Plot of ω vs. B")
xlabel("B")
ylabel("ω")
# PLOT: omega vs. C
figure()
a=1.74e-34
b=2.64e-5
C=10^-5:10^-5:10^-3
omega = [averageAngularVelocity(a, b, c) for c in C];
cf = plot(C, omega)
title("Plot of ω vs. C")
xlabel("C")
ylabel("ω")