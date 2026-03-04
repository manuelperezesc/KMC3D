# Kinetic MonteCarlo in 3D for amorphous films
 
using Distributions
using LinearAlgebra
using Base.Threads 
using Random

# Function definition
function readinput(input_file)

	input_string=readlines(input_file)

	npaths=parse(Int64,input_string[1])
	nwalks=parse(Int64,input_string[2])
	EF=parse.(Float64,split(input_string[3],','))
	target_hops=parse(Int64,input_string[4])
	temperature=parse(Float64,input_string[5])
	reorg=parse(Float64,input_string[6])
	dist=parse(Float64,input_string[7])
	print(split(input_string[8]," ")[1])
	tc=String(split(input_string[8]," ")[1])
	xc=parse.(Float64,split(input_string[8]," ")[2])
	yc=parse.(Float64,split(input_string[8]," ")[3])
	xcmax=parse.(Float64,split(input_string[8]," ")[4])
	xe=parse.(Float64,split(input_string[9]," ")[2])
	ye=parse.(Float64,split(input_string[9]," ")[3])
	xemax=parse.(Float64,split(input_string[9]," ")[4])

	return npaths,nwalks,EF,target_hops,temperature,reorg,dist,tc,xc,yc,xcmax,xe,ye,xemax
end


function MarcusRate(coup,ene,r,temperature)
        #constants
	hbar=6.58211928e-16
	kb=8.6173324e-5 
        #rate calculation
	k=(2.0*pi/hbar)*coup*coup/sqrt(4.0*pi*r*kb*temperature)*exp(-(ene+r)^2.0 / (4*r*kb*temperature))
        return k
end

function lorentzian_from_random(rn,x0,y0,maxval)
	val=x0+y0*tan(pi*(rn-0.5))
	return clamp(val,-maxval,maxval)
end

function normal_from_random(rn,x0,y0,xmax)
	val=x0+y0*rn
	return clamp(val,-xmax,xmax)
	#return x0+y0*rn
end

function double_exponential_from_random(rn, x0, y0, maxval)
    val = x0 - y0 * sign(rn - 0.5) * log(1 - 2 * abs(rn - 0.5))
    return clamp(val, -maxval, maxval)
end

function random_vector_on_sphere(magnitude,rng)
	u3=rand(rng)
	u4=rand(rng)

	x=magnitude*sin(pi*u3)*cos(2*pi*u4)
	y=magnitude*sin(pi*u3)*sin(2*pi*u4)
	z=magnitude*cos(pi*u3)

	#u=2.0*rand(rng)-1
	#phi=2*pi*rand(rng)

	#z=magnitude*u
	#r_xy=magnitude*sqrt(max(0.0,1.0-u^2.0))
	#x=r_xy*cos(phi)
	#y=r_xy*sin(phi)

	return[x,y,z]
end

function run_walk_double_exponential(npaths,efield,target_hops,temperature,reorg,dist,x,y,xmax,xx,yy,xxmax)
	nhops=0
	t=0.0
	displacement=[0.0,0.0,0.0]

	global_seed = rand(RandomDevice(),UInt64)   
	time_seed = UInt64(time_ns())               
	thread_seed = UInt64(threadid())            
	sum_seed = (global_seed + time_seed + thread_seed) % UInt64(typemax(Int))
    	final_seed = Int(sum_seed)                        
   	rng = MersenneTwister(final_seed)
	
	while nhops < target_hops
		
		random_numbers=rand(rng,npaths+2)		

		r1=random_numbers[1]
		r2=random_numbers[2]

		#Positions
		new_positions=Vector{Vector{Float64}}()
		                                                   
		for i in 1:npaths
			vi=random_vector_on_sphere(dist,rng)
			push!(new_positions,vi)
		end	

		rates=Float64[]
		a0=0.0	
	
		#Rates calculation
		for i in 1:npaths
			v=double_exponential_from_random(random_numbers[i+2],x,y,xmax)/1000.0
			ae=(clamp(xx+yy*randn(rng),-xxmax,xxmax))/1000.0+dot(efield.*1e-8,new_positions[i])
			k=MarcusRate(v,ae,reorg,temperature)*1e-12
			a0+=k
			push!(rates,k)
		end

		if a0==0
			break
		end

		#Time step and hopping condition	
		ts=(1.0/a0)*log(1.0/r1) # time step
		a0=a0*r2
		t+=ts

		#Hop evaluation
		sum1=0.0
		for i in 1:npaths
			sum1+=rates[i]
			if sum1 > a0

				nhops+=1.0
				
				xi=displacement[1]
				yi=displacement[2]
				zi=displacement[3]

				xf=xi+=new_positions[i][1]
				yf=yi+=new_positions[i][2]
				zf=zi+=new_positions[i][3]			

				displacement[1]=xf
				displacement[2]=yf
				displacement[3]=zf

				break
			end 
		end 
	end

	mobility=abs(dot(displacement,efield)*10.0^(-8.0)/(t*10.0^(-12.0)*dot(efield,efield)))

	return t, nhops, displacement, mobility
end


function run_walk_lorentzian(npaths,efield,target_hops,temperature,reorg,dist,x,y,xmax,xx,yy,xxmax)
	nhops=0
	t=0.0
	displacement=[0.0,0.0,0.0]

	global_seed = rand(RandomDevice(),UInt64)   
	time_seed = UInt64(time_ns())               
	thread_seed = UInt64(threadid())            
	sum_seed = (global_seed + time_seed + thread_seed) % UInt64(typemax(Int))
    	final_seed = Int(sum_seed)                        
   	rng = MersenneTwister(final_seed)
	
	while nhops < target_hops
		
		random_numbers=rand(rng,npaths+2)		

		r1=random_numbers[1]
		r2=random_numbers[2]

		#Positions
		new_positions=Vector{Vector{Float64}}()
		                                                   
		for i in 1:npaths
			vi=random_vector_on_sphere(dist,rng)
			push!(new_positions,vi)
		end	

		rates=Float64[]
		a0=0.0	
	
		#Rates calculation
		for i in 1:npaths
			v=lorentzian_from_random(random_numbers[i+2],x,y,xmax)/1000.0
			ae=(clamp(xx+yy*randn(rng),-xxmax,xxmax))/1000.0+dot(efield.*1e-8,new_positions[i])
			k=MarcusRate(v,ae,reorg,temperature)*1e-12
			a0+=k
			push!(rates,k)
		end

		if a0==0
			break
		end

		#Time step and hopping condition	
		ts=(1.0/a0)*log(1.0/r1) 
		a0=a0*r2
		t+=ts

		#Hop evaluation
		sum1=0.0
		for i in 1:npaths
			sum1+=rates[i]
			if sum1 > a0

				nhops+=1.0
				
				xi=displacement[1]
				yi=displacement[2]
				zi=displacement[3]

				xf=xi+=new_positions[i][1]
				yf=yi+=new_positions[i][2]
				zf=zi+=new_positions[i][3]			

				displacement[1]=xf
				displacement[2]=yf
				displacement[3]=zf

				break
			end 
		end 
	end

	mobility=abs(dot(displacement,efield)*10.0^(-8.0)/(t*10.0^(-12.0)*dot(efield,efield)))

	return t, nhops, displacement, mobility
end

function run_parallelization(nwalks,npaths,efield,target_time,temperature,reorg,dist,type,x,y,xmax,xx,yy,xxmax)

	times=zeros(nwalks)
	hops=zeros(nwalks)
	mobilities=zeros(nwalks)

	if (type == "l")
		@threads for walk_id in 1:nwalks
			time,nhops,dd,mobility=run_walk_lorentzian(npaths,efield,target_time,temperature,reorg,dist,x,y,xmax,xx,yy,xxmax)
			mobilities[walk_id]=mobility
			hops[walk_id]=nhops
			
		end

	else
		if (type == "d")
			@threads for walk_id in 1:nwalks
				time,nhops,dd,mobility=run_walk_double_exponential(npaths,efield,target_time,temperature,reorg,dist,x,y,xmax,xx,yy,xxmax)
				mobilities[walk_id]=mobility
				hops[walk_id]=nhops
				
			end
		end
	end

	return mobilities, hops

end

#Read input
n_paths,walks,EF,hops,temp,re,d,tc,xc,yc,xcmax,xe,ye,xemax=readinput(ARGS[1])
#Actual KMC calculation
total_mobilities,total_hops=run_parallelization(walks,n_paths,EF,hops,temp,re,d,tc,xc,yc,xcmax,xe,ye,xemax)
#Average mobility
av_mobility=(sum(total_mobilities)/length(total_mobilities))
av_hops=(sum(total_hops)/length(total_hops))

print(av_mobility,"\n")
print(av_hops,"\n")

