# if dydt=f

# y(next time step)=y(previous time step)+ f(previous time step)*dt
#If you go through Numerical Analysis by Burden, you'd see that when people numerically integrate stuff people like to use higher order method than the above which is simple forward euler
# type integration scheme which gives you only first order accuracy
# What im saying is
# y(t+h)=y(t)+dydt(t)*h+d2ydt2(t)*h^2+d3ydt3*h^3+O(h^4) by taylor series expansion
# you'll learn that the method that you implemented at wolfram and also on line 3 leaves you error with
#error= d2ydt2(t)*h^2+d3ydt3*h^3+O(h^4) meaning we have error that's  O(h^2) leaving this method with order of accuracy equal to first order
# that is right. By using other scheme, we want to increase the order of accuracy when we perform this error analysis with taylor series expansion
# If you want, you can use your own time to go over order of accuracy analysis for Runge kutta
#But main idea is that Runge kutta 4 RK4 gives you order of accuracy equal to 4.
#Error s O(h^4)
using LinearAlgebra
using PyPlot
include("withrefdistfunctions.jl")
include("getdata.jl")
# Characteristic Values
L0=3.844e7
G0=6.67e-11
M0=5.972e24
kappa0=8.99e9
q0=1.93e-3 # moon charge, calculated using Q=CV=4pie0RV
ka0=15
kc0=30
# T0=0.1382 days
moonmass=7.3477e22
# Initialization
N=4
T=100
Nt=10000
dt=T/Nt
kappa=1
G=1
vi=0.3
ui=Vector{Float64}()
for i=1:N
	push!(ui,cos(2*pi*i/N))
	push!(ui,sin(2*pi*i/N))
end
for i=1:N
	push!(ui,vi*-sin(2*pi*i/N))
	push!(ui,vi*cos(2*pi*i/N))
end
mi=moonmass/(M0*N)
qi=1/N
mlist=fill(mi,N)
qlist=fill(qi,N)
ka=15
kc=30
d=0.05
refdist=0
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

u=zeros(Int64(4*N))
# I created this ulist array so that we can store u from each time step
ulist=zeros(Int64(4*N),Nt)
u.=ui
#let's say we already have Nbodydynamics code built
for t=1:Nt
#Let's start with finding those k1k2k3k4
#dudt(y)=Nbodydynamicsfunction
# we don't have time dependence our function has format of Nbodydynamics(u,t,kappa,qlist,mlist)
# technically we don't even have explicit time dependence meaning it depends on the time in the sense of
#dydt=Nbodydynamics(y(t))
    k1=Nbodydynamics(u,t,kappa,G,ka,kc,d,qlist,mlist,refdist)
    k2=Nbodydynamics(u+dt/2*k1,t+dt/2,kappa,G,ka,kc,d,qlist,mlist,refdist)
    k3=Nbodydynamics(u+dt/2*k2,t+dt/2,kappa,G,ka,kc,d,qlist,mlist,refdist)
    k4=Nbodydynamics(u+dt*k3,t+dt,kappa,G,ka,kc,d,qlist,mlist,refdist)
	
# our goal is to end current time step by finding U at next time step
	# println(k1+2*k2+2*k3+k4)
    u.=u+dt/6*(k1+2*k2+2*k3+k4)
    ulist[:,t].=u
end

println(ulist[1:4,1]) # initial position
println(ulist[1:4,Nt]) # final position

tempx=zeros(Int64(N))
tempy=zeros(Int64(N))
# Graphing as pos of time
for t=1:50:Nt
	for i=1:N
		tempx[i]=ulist[2i-1,t]
		tempy[i]=ulist[2i,t]
	end
	clf()
	scatter(tempx,tempy)
	xlim(-10,10)
	ylim(-10,10)
	if t==1
		pause(8)
	else
		pause(0.1)
	end
end

pos=ulist[1:2,:]
figure()
plot(pos[1,:],pos[2,:])


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
		if i==1
			push!(xiavglist,xiavg)
			push!(yiavglist,yiavg)
		end
        center = [xiavg; yiavg]
        tempr = sqrt((center[1] - xi[dt])^2 + (center[2] - yi[dt])^2) # distance between (x1,y1) and center
        global r += tempr
    end
end
r = r / (N*Nt)
println("r=$r")
# GRAPH: x1-xavg, y1-yavg, as a function of time
x1=ulist[1,:]
y1=ulist[2,:]
deltax = x1 .- xiavglist
deltay = y1 .- yiavglist
figure()
plot(1:Nt, deltax,"b-") # blue = x
plot(1:Nt, deltay,"r-") # red = y
xlabel("t")
ylabel("Δ")
title("Deviation of coordinates from the Average Over Time")
show()