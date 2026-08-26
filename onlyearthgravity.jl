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
include("onlyearthgravityfunctions.jl")
#= STANDARD CONDITIONS
N=2
T=100
Nt=10000
dt=T/Nt
kappa=8.99*10^9
G=6.67*10^-11
ui=[-1;0;1;0;0;-0.3;0;0.3]
mlist=[1000,1000]
qlist=[5*10^(-4),5*10^(-4)]
ka=15
=#
# All 3 forces acting
N=4
T=100
Nt=10000
dt=T/Nt
kappa=8.99*10^9
G=6.67*10^-11
ui=[-1;0;1;0;0;-1;0;1;0;-0.3;0;0.3;0.3;0;-0.3;0]
mlist=[1000,1000,1000,1000]
qlist=[5*10^(-4),-5*10^(-4),5*10^(-4),-5*10^(-4)]
ka=15
# Circular motion
#ui=[-1;0;1;0;0;-0.3;0;0.3]
# without spring
#ka=0
# without gravity
#G=0
# without electric
#kappa=0


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
    k1=Nbodydynamics(u,t,kappa,G,ka,0,qlist,mlist)
    k2=Nbodydynamics(u+dt/2*k1,t+dt/2,kappa,G,ka,0,qlist,mlist)
    k3=Nbodydynamics(u+dt/2*k2,t+dt/2,kappa,G,ka,0,qlist,mlist)
    k4=Nbodydynamics(u+dt*k3,t+dt,kappa,G,ka,0,qlist,mlist)
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
	xlim(-30,30)
	ylim(-30,30)
	pause(0.1)
end

for i=1:N
	pos=ulist[2i-1:2i,:]
	figure()
	plot(pos[1,:],pos[2,:])
end