#=
Author: Yoav Dekel
=#

using LinearAlgebra, DifferentialEquations

# Orbit Struct Functionality
struct Orbit
    a::Float64
    e::Float64
    ex::Float64
    ey::Float64
    i::Float64
    Ω::Float64
    θ::Float64
    ω::Float64
end

function eMagOrbit(a,e,i,Ω,θ,ω)
    ex = e*cos(ω)
    ey = e*sin(ω)

    return Orbit(a,e,ex,ey,i,Ω,θ,ω)
end

function eVecOrbit(a,ex,ey,i,Ω,θ)
    e = norm([ex,ey])
    ω = atan(ey,ex)

    return Orbit(a,e,ex,ey,i,Ω,θ,ω)
end

# Orbit Functions
function Orbit2Cart(μ,orbit::Orbit)
    ν = orbit.θ - orbit.ω
    p = orbit.a*(1-orbit.e^2)

    # orbital plane to cartesian plane DCM
    DCM = EA2DCM([orbit.Ω,orbit.i,orbit.ω],[3,1,3])

    r = [cos(ν),sin(ν),0]
    r = r*p/(1 + orbit.e*cos(ν))
    r = DCM*r # orbital to cartesian

    v = [-sin(ν),orbit.e + cos(ν),0]
    v = v*sqrt(μ/p)
    v = DCM*v # orbital to cartesian

    return [r;v]
end

function Cart2Orbit(μ,cart)
    r = cart[1:3]
    v = cart[4:6]

    h = cross(r,v)

    e = cross(v,h)/μ - r/norm(r)

    if e != 0
        ê = e/norm(e)
    else
        ê = [1 0 0]
    end

    a = (norm(h)^2/μ)/(1-norm(e)^2)

    ẑ = [0,0,1]

    ĥ = h/norm(h)

    n = cross(ẑ,ĥ)

    if norm(n) != 0
        n̂ = n/norm(n)
    else
        n̂ = [1,0,0]
    end

    i = acos(h[3]/norm(h))
    Ω = atan(ĥ[1],-ĥ[2])

    cosω = ê[1]*n̂[1] + ê[2]*n̂[2]
    sinω = dot(cross(n̂,ê),ĥ)
    ω = atan(sinω,cosω)

    cosν = dot(ê,r)/norm(r)
    sinν = dot(cross(ê,r),ĥ)/norm(r)
    ν = atan(sinν,cosν)

    θ = ω + ν

    return eMagOrbit(a,norm(e),i,Ω,θ,ω)
end

function Orbit2DCM(orbit::Orbit)
    return EA2DCM([orbit.Ω,orbit.i,orbit.θ],[3,1,3])
end

function Cart2DCM(μ,cart)
    orbit = Cart2Orbit(μ,cart)

    return Orbit2DCM(orbit)
end

function Kepler(M,e;ϵ = 1e-12,nᵢ = 1e3)
    # ϵ- tolerance
    # nᵢ- maximum iterations

    M = mod(M,2*pi)
    E = M

    for _ in 1:nᵢ
        f = M - E + e*sin(E)
        if abs(f) < ϵ
            cosν = (cos(E) - e)/(1 - e*cos(E))
            sinν = (sqrt(1 - e^2)*sin(E))/(1 - e*cos(E))
            return mod(atan(sinν,cosν),2pi)
        end

        Δf = -1 + e*cos(E)
        E = E - f/Δf
    end

    println("Error in Kepler(): Max iterations ($nᵢ) reached")
    return nothing # error
end

# Orbit Utils

function EA2DCM(EA,sequence)
    DCM = Matrix(1.0I,3,3)

    for i in 3:-1:1
        DCM = RFunction(EA[i],sequence[i])*DCM
    end

    return DCM
end

function GetPeriod(μ,orbit::Orbit)
    return 2pi*sqrt(orbit.a^3/μ)
end

# Satellite Constellation Utils

# returns an array of Ω and M
function LatticeFlower2D(Ns, Np = 1, Nc = 0)
    # malloc
    constellation = zeros(Ns,2)

    sat = 1

    for i in 1:Np
        for j in 1:Ns/Np
            constellation[sat,1] = mod(2pi*(i-1)/Np,2pi)
            constellation[sat,2] = mod(2pi*Np/Ns*(j-1) - 2pi*Nc/Ns*(i-1),2pi)

            sat += 1;
        end
    end

    return constellation
end

# creating Satellies vector, which holds the Orbit struct of eat satellite
function ConToSats(Ns, constellation, e = 0, i = 0, ω = 0)
    θ = zeros(Ns)

    for idx in 1:Ns
        M = constellation[idx,2]
        ν = Kepler(M,e)
        θ[idx] = ω + ν
    end

    return eMagOrbit.(a,e,i,constellation[:,1],θ,ω);
end

#= Convention for this function is

[DCM] = cart --> cart'

cart' = [DCM]*cart

Where all basis vectors are column vectors

=#

function RFunction(angle,axis)
    if axis == 1
        return [1 0 0;0 cos(angle) -sin(angle);0 sin(angle) cos(angle)]
    elseif axis == 2
        return [cos(angle) 0 sin(angle);0 1 0;-sin(angle) 0 cos(angle)]
    elseif axis == 3
        return [cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1]
    end
end

# Propagators
"""
    This function is an Orbit Propagator for a spacecraft under Keplerian motion

    This must be called with the ODEProblem() function, as it is using the ! syntax with du inside the function

    Inputs: (units are not important as long as consistent)
        du = vector that is derivative of u vector (relavent for calling in a DifferentialEquations.jl environment)
        u = vector with components: [r₁,r₂,r₃,v₁,v₂,v₃], the position and velocity vectors (in x y z Cartesian cordinates)
        μ = magnitude of the Gravitational parameter
        t = time (relavent for calling in a DifferentialEquations.jl environment)
"""

function OrbitPropCartesian!(du,u,μ,t)
    x,y,z = u

    du[1:3] = u[4:6]

    r = sqrt(x^2 + y^2 + z^2)

    du[4] = -μ*x/r^3
    du[5] = -μ*y/r^3
    du[6] = -μ*z/r^3
end

"""
    This function is an Orbit Propagator for a spacecraft under J2 perturbations where the force model
        is using Cartesian coordinates

    This must be called with the ODEProblem() function, as it is using the ! syntax with du inside the function

    Inputs: (units are not important as long as consistent)
        du = vector that is derivative of u vector (relavent for calling in a DifferentialEquations.jl environment)
        u = vector with components: [r₁,r₂,r₃,v₁,v₂,v₃], the position and velocity vectors (in x y z Cartesian cordinates)
        consts = constants in the problem: [μ,J₂,R] which are the Gravitational parameter, J₂ pertubation constant and the radius of the planet
        t = time (relavent for calling in a DifferentialEquations.jl environment)
"""
function OrbitJ2PropCartesian!(du,u,consts::Tuple{Number,Number,Number},t)
    x,y,z = u
    μ,J₂,R = consts

    du[1:3] = u[4:6]

    r = sqrt(x^2 + y^2 + z^2)

    cc = (-3/2)μ*J₂*R^2/r^5 # placeholder

    du[4] = (-μ/r^3 + cc*(1 - 5z^2/r^2))x
    du[5] = (-μ/r^3 + cc*(1 - 5z^2/r^2))y
    du[6] = (-μ/r^3 + cc*(3 - 5z^2/r^2))z
end

"""
    This function is an Orbit Propagator for a spacecraft under J2 perturbations where the force model
        is using Cartesian coordinates

    This must be called with the ODEProblem() function, as it is using the ! syntax with du inside the function

    Inputs: (units are not important as long as consistent)
        du = vector that is derivative of u vector (relavent for calling in a DifferentialEquations.jl environment)
        u = vector with components: [a, e, i, Ω, ω, M], the position and velocity vectors (in x y z Cartesian cordinates)
        consts = constants in the problem: [μ,J₂,R] which are the Gravitational parameter, J₂ pertubation constant and the radius of the planet
        t = time (relavent for calling in a DifferentialEquations.jl environment)
"""
function OrbitJ2PropKeplerian!(du,u,consts,t)
    a,e,i,Ω,ω,M = u
    μ,J₂,R = consts
    
    ν = Kepler(M,e)

    # orbital plane to cartesian plane DCM
    DCM = EA2DCM([Ω,i,ν+ω],[3,1,3])

    p = a*(1 - e^2)
    b = a*sqrt(1 - e^2)
    h = sqrt(p*μ)
    r = p/(1 + e*cos(ν))

    f0 = [zeros(5,1);sqrt(μ/a^3)]

    B = (1/h)*[2e*(a^2)sin(ν) (2a^2)p/r 0;
        p*sin(ν) (p+r)cos(ν)+r*e 0;
        0 0 cos(ν+ω)r;
        0 0 sin(ν+ω)r/sin(i);
        -cos(ν)p/e (p+r)sin(ν)/e -sin(ν+ω)r/tan(i);
        b*cos(ν)p/(a*e)-2b*r/a -(p+r)sin(ν)b/(a*e) 0]

    x,y,z = DCM*[r,0,0] # orbital plane to cartesian

    aJ₂ = zeros(3,1) # malloc

    cc = (-3/2)μ*J₂*R^2/r^5 # placeholder

    aJ₂[1] = cc*(1 - 5z^2/r^2)x
    aJ₂[2] = cc*(1 - 5z^2/r^2)y
    aJ₂[3] = cc*(3 - 5z^2/r^2)z

    du[1:6] = f0 + B*(DCM'*aJ₂) # converting aJ₂ from cartesian to orbital inside
end

# Models

function densityModel(alt)
    # model source : Vallado, Fundamentals of Astrodynamics and Applications
    if alt > 1000
        h₀ = 1000;
        ρ₀ = 3.019e-15;
        H = 268.00;
    elseif alt > 900
        h₀ = 900;
        ρ₀ = 5.245e-15;
        H = 181.05;
    elseif alt > 800
        h₀ = 800;
        ρ₀ = 1.170e-14;
        H = 125.64;
    elseif alt > 700
        h₀ = 700;
        ρ₀ = 3.614e-14;
        H = 88.667;
    elseif alt > 600
        h₀ = 600;
        ρ₀ = 1.454e-13;
        H = 71.835;
    elseif alt > 500
        h₀ = 500;
        ρ₀ = 6.967e-13;
        H = 63.822;
    elseif alt > 450
        h₀ = 450;
        ρ₀ = 1.585e-12;
        H = 60.828;
    elseif alt > 400
        h₀ = 400;
        ρ₀ = 3.725e-12;
        H = 58.515;
    elseif alt > 350
        h₀ = 350;
        ρ₀ = 9.518e-12;
        H = 53.298;
    elseif alt > 300
        h₀ = 300;
        ρ₀ = 2.418e-11;
        H = 53.628;
    elseif alt > 250
        h₀ = 250;
        ρ₀ = 7.248e-11;
        H = 45.546;
    elseif alt > 200
        h₀ = 200;
        ρ₀ = 2.789e-10;
        H = 37.105;
    elseif alt > 180
        h₀ =  180;
        ρ₀ = 5.464e-10;
        H = 29.740;
    elseif alt > 150
        h₀ = 150;
        ρ₀ = 2.07e-9;
        H = 22.523;
    elseif alt > 140
        h₀ = 140;
        ρ₀ = 3.845e-9;
        H = 16.149;
    elseif alt > 130
        h₀ = 130;
        ρ₀ = 8.484e-9;
        H = 12.636;
    elseif alt > 120
        h₀ = 120;
        ρ₀ = 2.438e-8;
        H = 9.473;
    elseif alt > 110
        h₀ = 110;
        ρ₀ = 9.661e-8;
        H = 7.263;
    else
        h₀ = 100;
        ρ₀ = 5.297e-7;
        H = 5.877;
    end

    # calculated density
    ρ =  ρ₀*exp((h₀-alt)/H) # [kg/m^3]

    # returns density in [kg/km^3]
    return ρ*1e9
end