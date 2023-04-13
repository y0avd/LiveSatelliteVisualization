# Initializations
using GLMakie, FileIO, GeometryBasics, Colors, LinearAlgebra, DifferentialEquations

include("OrbitFunctions.jl")
include("PlottingFunctions.jl")

# constants
const μ = 398600.4418 # [km^3/s^2] Earth's graviation parameter
const R⨁ = 6378.137 # [km] Earth's equatorial radius
const ω⨁ = 2pi/(24*3600); # [rad/sec] Earth angular rotation

# config
screen_res = (1920, 1080)
sats_opacity = 1
orbits_opacity = 0.25

function getNp(Ns, percent)
    # for a given Ns, Np can only be factor of Ns
    Np_options = []

    for idx in 1:Integer(floor(Ns/2))
        if mod(Ns,idx) == 0
            push!(Np_options, idx)
        end
    end

    push!(Np_options, Ns)

    # once we have Np options, we interpolate the percentage into one of those choices
    idx = Integer(ceil(percent*length(Np_options)/100))

    if idx < 1; idx = 1
    elseif idx > length(Np_options); idx = length(Np_options)
    end

    return Np_options[idx]
end

function getNc(Np, percent)
    # for a given Np, Nc can be between 0 and Np-1
    Nc_options = 0:1:(Np-1)

    # once we have Np options, we interpolate the percentage into one of those choices
    idx = Integer(ceil(percent*length(Nc_options)/100))

    if idx < 1; idx = 1
    elseif idx > length(Nc_options); idx = length(Nc_options)
    end

    return Nc_options[idx]
end

# Observables from interactive (sets default)
slider_vals_obs = Observable((10, 50., 25., 1.5, 30., 0.))

# Observables from satellite
Ns = Observable(slider_vals_obs[][1])
Np = Observable(getNp(Ns[],slider_vals_obs[][2]))
Nc = Observable(getNc(Np[],slider_vals_obs[][3]))
a = Observable(R⨁*slider_vals_obs[][4])
i = Observable(deg2rad(slider_vals_obs[][5]))
e = Observable(slider_vals_obs[][6])

# Observables from plotting
azimuth = Observable(0.4)
elevation = Observable(0.75)

# Observables for plotting
title_strA = Observable("")
title_strB = Observable("")
sats_color = Observable((:blue, sats_opacity))
orbits_color = Observable((:blue, orbits_opacity))

## PLOTTING ##

f = Figure(resolution = screen_res)
display(f)

function updateFigure()
    # save camera angles
    azimuth[] = f.content[2].azimuth[]
    elevation[] = f.content[2].elevation[]
    empty!(f)
    setUp(slider_vals_obs[])
end

function Sliders(slider_vals)
    # Sliders
    sg = SliderGrid(
        f[2, 2:3],
        (label = "Number of Satellites", range = 1:1:500, format = "{:d}", startvalue = slider_vals[1]),
        (label = "Number of Planes", range = 0:.1:100, format = "{:.1f} %", startvalue = slider_vals[2]),
        (label = "Configuration", range = 0:.1:100, format = "{:.1f} %", startvalue = slider_vals[3]),
        (label = "Semi-major Axis", range = 1:.001:10, format = "{:.3f} R⨁", startvalue = slider_vals[4]),
        (label = "Inclination", range = 0:.01:90, format = "{:.2f} ᵒ", startvalue = slider_vals[5]),
        (label = "Eccentricity", range = 0:.01:0.9, format = "{:.2f} ", startvalue = slider_vals[6]),
        halign = :right
    )

    # update Observables when slider used
    sliderobservables = [s.value for s in sg.sliders]
    lift(sliderobservables...) do slvalues...
        slider_vals_obs[] = slvalues

        Ns[] = slvalues[1]
        Np[] = getNp(Ns[], slvalues[2])
        Nc[] = getNc(Np[], slvalues[3])

        a[] = R⨁*slvalues[4]
        i[] = deg2rad(slvalues[5])
        e[] = slvalues[6]

        title_strB[] = "Ns = $(Ns[]) Np = $(Np[]) Nc = $(Nc[])"
    end
end

function Buttons()
    # Buttons
    f[2, 1] = menusgrid = GridLayout()

    button = menusgrid[1:3, 1] = [
        Button(f, label="Update\n3D Plot", buttoncolor = RGBf(0.8, 0.94, 0.8)),
        Button(f, label="Satellites\nToggle", buttoncolor = RGBf(0.94, 0.94, 0.94)),
        Button(f, label="Orbits\nToggle", buttoncolor = RGBf(0.94, 0.94, 0.94))]
    
    # Button Uses
    on(button[1].clicks) do click
        updateFigure()
    end

    on(button[2].clicks) do click
        if sats_color[][2] == sats_opacity
            sats_color[] = (sats_color[][1],0)
        else
            sats_color[] = (sats_color[][1],sats_opacity)
        end
    end

    on(button[3].clicks) do click
        if orbits_color[][2] == orbits_opacity
            orbits_color[] = (orbits_color[][1],0)
        else
            orbits_color[] = (orbits_color[][1],orbits_opacity)
        end
    end
end

function SatToCart(orbit::Orbit)
    cart = Orbit2Cart(μ, orbit)

    return Point3f(cart[1:3])
end

function setUp(slider_vals)
    # creating observables for plotting
    cartPos = Observable([Vector{Point{3, Float32}}()])

    # Layout Content
    # 2D
    gB = f[1, 3] = GridLayout()
    axB = Axis(gB[1,1], title = title_strB,
        xlabel = "Ω (degrees)", xticks = 0:30:360,
        ylabel = "M (degrees)", yticks = 0:30:360,
        limits = (-5, 360, -5, 360),
        xpanlock = true, xzoomlock = true, xrectzoom = false,
        ypanlock = true, yzoomlock = true, yrectzoom = false,
        aspect = 1
    )

    # 3D Plot
    gA = f[1, 1:2] = GridLayout()
    axA = Axis3(gA[1,1], title = title_strA,
        aspect = :data,
        azimuth = azimuth[],
        elevation = elevation[]
    );

    # Prettify
    hidespines!(axA)
    hidedecorations!(axA)
    hidespines!(axB)

    Sliders(slider_vals)

    # calculate constellation
    constellation = @lift (LatticeFlower2D($Ns, $Np, $Nc))
    con_deg = @lift (rad2deg.($constellation))

    # plot constellation
    scatter!(axB, con_deg, color=:blue);

    # Now we must get all satellite positions and save them as an observable from constellation
    sats = Observable(Vector{Orbit}())

    cartPos = @lift SatToCart.($sats)

    axA
    EarthPlot()

    scatter!(axA, cartPos, color=sats_color)

    θ = zeros(Ns[])

    for idx in 1:Ns[]
        M = constellation[][idx,2]
        θ[idx] = Kepler(M,e[])
    end

    sats[] = eMagOrbit.(a[],e[],i[],constellation[][:,1],θ,0);

    for idx in 1:Ns[]
        # Plotting Orbits for one satellite in every plane
        if mod(idx,Ns[]/Np[]) == 0
            PlotKeplerianOrbit(to_value(sats)[idx], color = orbits_color)
        end
    end

    Buttons()

    reset_limits!(axA)
end

setUp(slider_vals_obs[])