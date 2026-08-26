using LinearAlgebra
using CSV
using DataFrames

include("getdata.jl")
nd=(3.844e7)^2/(6.67e-11*5.972e24) # nondimensionalizer

function Nbodydynamics(u,t,kappa,G,ka,kc,d,qlist,mlist,refdist=0)
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
	N=length(qlist)
	Nacceleration=zeros(2*N)
	for i=1:N
	# we are trying to find overall force that results on ith charge
		overallaxi=0
		overallayi=0
		xi=u[Int64(2*(i-1)+1)]
		yi=u[Int64(2*(i))]
		qi=qlist[i]
		mi=mlist[i]
		for j=1:N
			if i!=j
				xj=u[Int64(2*(j-1)+1)]
				yj=u[Int64(2*(j))]
				#rij is a unit vector in the direction of ij
				rij=[xj-xi;yj-yi]
				rmag=norm(rij)
				rij=rij/rmag
				# how do we find qi and qj? qi=qlist[i] qj=qlist[j]
				# how do we find mi? mi=mlist[i]
				qj=qlist[j]
				mj=mlist[j]
				# electric force
				tempaxielectric=-kappa*qi*qj/rmag^2*rij[1]/mi
				tempayielectric=-kappa*qi*qj/rmag^2*rij[2]/mi
				overallaxi=overallaxi+tempaxielectric
				overallayi=overallayi+tempayielectric
			end
			# ka spring force
			if i==1	
				xibefore=u[Int64(2*(N-1)+1)]
				yibefore=u[Int64(2*(N))]
			else
				xibefore=u[Int64(2*(i-1))-1]
				yibefore=u[Int64(2*(i-1))]
			end
			if i==N
				xinext=u[1]
				yinext=u[2]
			else
				xinext=u[Int64(2*(i)+1)]
				yinext=u[Int64(2*(i)+2)]
			end
			riinext=[xinext-xi;yinext-yi]
			rmagnext=norm(riinext)
			riinext=riinext/rmagnext
			riibefore=[xibefore-xi;yibefore-yi]
			rmagbefore=norm(riibefore)
			riibefore=riibefore/rmagbefore

			tempaxispring=ka*(rmagnext-refdist)*riinext[1]/mi
			tempayispring=ka*(rmagnext-refdist)*riinext[2]/mi
			overallaxi=overallaxi+tempaxispring
			overallayi=overallayi+tempayispring
			tempaxispring=ka*rmagbefore*riibefore[1]/mi
			tempayispring=ka*rmagbefore*riibefore[2]/mi
			overallaxi=overallaxi+tempaxispring
			overallayi=overallayi+tempayispring
		end
		# kc spring force
		ri11=[0-xi;0-yi] # center at [0,0]
		ri11mag=norm(ri11)
		ri11=ri11/ri11mag
		tempaxispring=kc*(ri11mag-refdist)*ri11[1]/mi
		tempayispring=kc*(ri11mag-refdist)*ri11[2]/mi
		overallaxi=overallaxi+tempaxispring
		overallayi=overallayi+tempayispring
        # add earth and sun based on data
        m12=1 # earth
        m13=1.989e30/M0 # sun
        ri12=[earthposrel[round(Int, t)][1]-xi;earthposrel[round(Int, t)][2]-yi]
        ri13=[sunposrel[round(Int, t)][1]-xi;sunposrel[round(Int, t)][2]-yi]
        ri12mag=norm(ri12)
        ri13mag=norm(ri13)
        ri12=ri12/ri12mag
        ri13=ri13/ri13mag
        tempaxigravity=G*mi*m12/ri12mag^2*ri12[1]/mi
		tempayigravity=G*mi*m12/ri12mag^2*ri12[2]/mi
		overallaxi=overallaxi+tempaxigravity
		overallayi=overallayi+tempayigravity
        tempaxigravity=G*mi*m13/ri13mag^2*ri13[1]/mi
		tempayigravity=G*mi*m13/ri13mag^2*ri13[2]/mi
		overallaxi=overallaxi+tempaxigravity
		overallayi=overallayi+tempayigravity
		# damping
		vdotri11=u[2*N+(2*i-1)]*ri11[1]+u[2*N+(2*i)]*ri11[2]
		overallaxi=overallaxi-d*(vdotri11)*ri11[1]
		overallayi=overallayi-d*(vdotri11)*ri11[2]
        # before setting to acceleration, multiply by nondimensionalizer nd=L^2/(GM)~=3.70954911 to make dimensionless
		Nacceleration[Int64(2*(i-1)+1)]=overallaxi
		Nacceleration[Int64(2*(i))]=overallayi
	end
	dudt=[u[Int64(2*N+1):end];Nacceleration]
	#u is stored as [x1;y1;x2;y2;.......xn;yn;dx1dt;dy1dt;dx2dt;dy2dt;,,,,,,,dxndt;dyndt]
	# dudt should be [dx1dt;dy1dt;dx2dt;dy2dt;.......dxndt;dyndt;d2x1dt2;d2y1dt2;.....d2xndt2;d2yndt2]
	return dudt
end