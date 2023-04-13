## Using Makie ECOSYSTEM for plotting
# currently only tested and working on GLMakie backend

# imports
using GLMakie, FileIO, GeometryBasics, Colors

# uses OrbitFunctions
include("OrbitFunctions.jl")

# All functions assume figure and 3D axis already created

# Plots 3D Earth with graphic, returns mesh for manipulation (like spinning)
function EarthPlot(ax)
    # Earth Image
    earth_img = load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"));

    # Earth 3D Plot
    earth_mesh = mesh!(ax,
        Sphere(Point3f0(0), R⨁),
        color = earth_img,
        shading = false,
    )

    return earth_mesh
end

# given a single point in an orbit, plot entire Keplerian orbit for 1 period
# uses DifferentialEquations
function PlotKeplerianOrbit(orbit::Orbit; color=(:blue, 0.25))
    P = GetPeriod(μ,orbit)

    u₀ = Orbit2Cart(μ, orbit)

    tspan = [0., P]

    prob = ODEProblem(OrbitPropCartesian!,u₀,tspan,μ)
    sol = solve(prob,Tsit5(),reltol = 1e-12,abstol = 1e-12)

    # malloc
    cart = zeros(length(sol.t),3) # cartesian positions

    for idx in 1:length(sol.t)
        cart[idx,:] = sol.u[idx][1:3]
    end

    lines!(cart, color = color);
end

function PlotSatellite(orbit::Orbit; color=(:blue, 0.25))
    cart = Orbit2Cart(μ, orbit)

    scatter!(cart[1], cart[2], cart[3], color = color);
end

# Rotates Earth by timestep [seconds]
function RotateEarth(earth_mesh, time)
    GLMakie.rotate!(earth_mesh, Vec3f(0,0,1),ω⨁*time)
end

# Plots 3D Satellite Constellation
# Optional timestep
function Plot3DSatCon(satellites; title_str = "3D Satellite Constellation")
    f = Figure(resolution=(1920,1080))
    ax = Axis3(
        f[1,1],
        title = title_str,
        xlabel = "x [km]",
        ylabel = "y [km]",
        zlabel = "z [km]",
        aspect = :data,
        azimuth = pi/4,
        elevation = pi/4
    );

    # Earth 3D Plot
    EarthPlot()

    for idx in 1:Ns
        # Plotting Orbits for one satellite in every plane
        if mod(idx,Ns/Np) == 0
            PlotKeplerianOrbit(satellites[idx])
        end

        # Plotting Satellite Initial Positions
        PlotSatellite(satellites[idx])
    end

    return f
end

# Creates animation of Satellite Constellation for specified time
# specify time in seconds (default = 1 hour)
function Animate3DSatCon(satellites; title_str = "3D Satellite Constellation Animation", time = 3600, nump::Int64 = 100, framerate::Int64 = 24)
    ## Initialization
    dt = time/nump

    nump += 1 # accounts for initial state

    # creating observables for plotting
    sats = Observable(Vector{Point{3, Float32}}())

    title_str_update = Observable("")
    azimuth = Observable(0.0)

    ## Position Calculations

    # malloc
    cart = zeros(nump,3,Ns)

    for idx1 in 1:Ns
        u₀ = Orbit2Cart(μ, satellites[idx1])

        tspan = [0., time]

        prob = ODEProblem(OrbitPropCartesian!,u₀,tspan,μ)
        sol = solve(prob,Tsit5(),reltol = 1e-12,abstol = 1e-12,saveat = dt)

        for idx2 in 1:nump
            cart[idx2,:,idx1] = sol.u[idx2][1:3]
        end
    end

    ## Plotting
    # static plotting
    
    f = Figure(resolution=(1920,1080))
    ax = Axis3(
        f[1,1],
        title = title_str_update,
        xlabel = "x [km]",
        xlabeloffset = 80,
        xlabelrotation = 0,
        ylabel = "y [km]",
        ylabeloffset = 80,
        ylabelrotation = 0,
        zlabel = "z [km]",
        zlabeloffset = 80,
        zlabelrotation = 0,
        aspect = :data,
        viewmode = :fit,
        titlegap = -100.0,
        azimuth = azimuth,
        elevation = pi/4
    );

    # Earth 3D Plot
    earth_mesh = EarthPlot()

    for idx in 1:Ns
        # Plotting Orbits for one satellite in every plane
        if mod(idx,Ns/Np) == 0
            PlotKeplerianOrbit(satellites[idx])
        end
    end

    scatter!(sats)

    ## Recording Animation
    record(f, "SatConAnimation.mp4", 1:nump; framerate = framerate) do frame
        temp_sats = Vector{Point{3, Float32}}()
        for idx in 1:Ns
            push!(temp_sats,Point3f(cart[frame,:,idx]))
        end

        sim_time = round((frame-1)*dt/60,digits=2)
        frame_time = round(dt,digits=2)

        RotateEarth(earth_mesh, (frame-1)*dt)

        if rot   
            azimuth[] = ω⨁*(frame-1)*dt
        end
        
        title_str_update[] = title_str*"\nFrametime: $frame_time Seconds\nTime: $sim_time Minutes"
        sats[] = temp_sats
    end
end

# 2D representation of Satellite Constellation
function Plot2DSatCon(constellation; title_str = "2D Satellite Constellation")
    return scatter(rad2deg.(constellation);
        figure = (; resolution = (800,800)),
        axis = (; title = title_str,
            xlabel = "Ω (degrees)",
            ylabel = "M (degrees)",
            limits = (0, 360, 0, 360))
    );
end
