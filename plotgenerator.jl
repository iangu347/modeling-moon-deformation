using LinearAlgebra
using PyPlot
include("withrefdistfunctions.jl")

# GRAPH: ka
kavalues = Float64[]
karvalues = Float64[]
r = 0
for i = 0:100 # ka = 0 to 100
    N = 4
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 8.99 * 10^9
    G = 6.67 * 10^-11
    ui = [-1; 0; 0; 1; 1; 0; 0; -1; 0; -0.3; -0.3; 0; 0; 0.3; 0.3; 0]
    mlist = [1000, 1000, 1000, 1000]
    qlist = [5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4)]
    ka = i
    kc = 15
    d = 0
    refdist = 0
    u = zeros(Int64(4 * N))
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    for t = 1:Nt
        k1 = Nbodydynamics(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k2 = Nbodydynamics(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k3 = Nbodydynamics(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k4 = Nbodydynamics(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end


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
    println("ka = $(ka) | r = $r")
    push!(kavalues, ka)
    push!(karvalues, r)
    r = 0
end
println("=======================")

figure()
plot(kavalues, karvalues)
xlabel("ka")
ylabel("r")
title("Average Orbit Radius r vs. Spring Constant ka")
show()

# GRAPH: kc
kcvalues = Float64[]
kcrvalues = Float64[]
r = 0
for i = 0:100 # kc = 0 to 100
    N = 4
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 8.99 * 10^9
    G = 6.67 * 10^-11
    ui = [-1; 0; 0; 1; 1; 0; 0; -1; 0; -0.3; -0.3; 0; 0; 0.3; 0.3; 0]
    mlist = [1000, 1000, 1000, 1000]
    qlist = [5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4)]
    ka = 15
    kc = i
    d = 0
    refdist = 0
    u = zeros(Int64(4 * N))
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    for t = 1:Nt
        k1 = Nbodydynamics(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k2 = Nbodydynamics(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k3 = Nbodydynamics(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k4 = Nbodydynamics(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

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
    println("kc = $(kc) | r = $r")
    push!(kcvalues, kc)
    push!(kcrvalues, r)
    r = 0
end
println("=======================")

figure()
plot(kcvalues, kcrvalues)
xlabel("kc")
ylabel("r")
title("Average Orbit Radius r vs. Spring Constant kc")
show()

# GRAPH: d
dvalues = Float64[]
drvalues = Float64[]
r = 0
for i = 0:0.1:10 # d = 0 to 10
    N = 4
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 8.99 * 10^9
    G = 6.67 * 10^-11
    ui = [-1; 0; 0; 1; 1; 0; 0; -1; 0; -0.3; -0.3; 0; 0; 0.3; 0.3; 0]
    mlist = [1000, 1000, 1000, 1000]
    qlist = [5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4)]
    ka = 15
    kc = 15
    d = i
    refdist = 0
    u = zeros(Int64(4 * N))
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    for t = 1:Nt
        k1 = Nbodydynamics(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k2 = Nbodydynamics(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k3 = Nbodydynamics(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k4 = Nbodydynamics(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

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
    println("d = $(d) | r = $r")
    push!(dvalues, d)
    push!(drvalues, r)
    r = 0
end
println("=======================")

figure()
plot(dvalues, drvalues)
xlabel("d")
ylabel("r")
title("Average Orbit Radius r vs. Damping Factor d")
show()

# GRAPH: m
mvalues = Float64[]
mrvalues = Float64[]
r = 0
for i = 100:100:10000 # m = 100 to 10000
    N = 4
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 8.99 * 10^9
    G = 6.67 * 10^-11
    ui = [-1; 0; 0; 1; 1; 0; 0; -1; 0; -0.3; -0.3; 0; 0; 0.3; 0.3; 0]
    mlist = [i, i, i, i]
    qlist = [5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4)]
    ka = 15
    kc = 15
    d = 0
    refdist = 0
    u = zeros(Int64(4 * N))
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    for t = 1:Nt
        k1 = Nbodydynamics(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k2 = Nbodydynamics(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k3 = Nbodydynamics(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k4 = Nbodydynamics(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

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
    println("m = $(mlist[1]) | r = $r")
    push!(mvalues, mlist[1])
    push!(mrvalues, r)
    r = 0
end
println("=======================")

figure()
plot(mvalues, mrvalues)
xlabel("m")
ylabel("r")
title("Average Orbit Radius r vs. Mass m")
show()

# GRAPH: q
qvalues = Float64[]
qrvalues = Float64[]
for i = 10^(-4):10^-5:10^(-3) # q = 10^-4 to 10^-3
    N = 4
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 8.99 * 10^9
    G = 6.67 * 10^-11
    ui = [-1; 0; 0; 1; 1; 0; 0; -1; 0; -0.3; -0.3; 0; 0; 0.3; 0.3; 0]
    mlist = [1000, 1000, 1000, 1000]
    qlist = [i, i, i, i]
    ka = 15
    kc = 15
    d = 0
    refdist = 0
    u = zeros(Int64(4 * N))
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    for t = 1:Nt
        k1 = Nbodydynamics(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k2 = Nbodydynamics(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k3 = Nbodydynamics(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k4 = Nbodydynamics(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

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
    println("q = $(qlist[1]) | r = $r")
    push!(qvalues, qlist[1])
    push!(qrvalues, r)
    r = 0
end
println("=======================")

figure()
plot(qvalues, qrvalues)
xlabel("q")
ylabel("r")
title("Average Orbit Radius r vs. Charges q")
show()

# GRAPH: refdist
refdistvalues = Float64[]
refdistrvalues = Float64[]
r = 0
for i = 0:0.1:10 # refdist = 0 to 10
    N = 4
    T = 100
    Nt = 10000
    dt = T / Nt
    kappa = 8.99 * 10^9
    G = 6.67 * 10^-11
    ui = [-1; 0; 0; 1; 1; 0; 0; -1; 0; -0.3; -0.3; 0; 0; 0.3; 0.3; 0]
    mlist = [1000, 1000, 1000, 1000]
    qlist = [5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4), 5 * 10^(-4)]
    ka = 15
    kc = 15
    d = 0
    refdist = i
    u = zeros(Int64(4 * N))
    ulist = zeros(Int64(4 * N), Nt)
    u .= ui
    for t = 1:Nt
        k1 = Nbodydynamics(u, t, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k2 = Nbodydynamics(u + dt / 2 * k1, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k3 = Nbodydynamics(u + dt / 2 * k2, t + dt / 2, kappa, G, ka, kc, d, qlist, mlist, refdist)
        k4 = Nbodydynamics(u + dt * k3, t + dt, kappa, G, ka, kc, d, qlist, mlist, refdist)
        u .= u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ulist[:, t] .= u
    end

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
    println("refdist = $(refdist) | r = $r")
    push!(refdistvalues, refdist)
    push!(refdistrvalues, r)
    r = 0
end
println("=======================")

figure()
plot(refdistvalues, refdistrvalues)
xlabel("Relaxed Spring Length")
ylabel("r")
title("Average Orbit Radius r vs. Relaxed Spring Length")
show()