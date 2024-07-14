module MinPathAnimation
using MinPath: MinPath, Coefficients, action, ∇action, K_energy, U, build_path, extend_to_period
using GLMakie, LaTeXStrings, LinearAlgebra, OffsetArrays

"""
    tolatex(x::Float64)::LaTeXString

Convert a Float64 number to a LaTeXString using scientific notation.
The number is rounded to 6 significant digits.

# Examples
```jldoctest
julia> tolatex(1.23456789e-10)
L"1.234567\times 10^{-10}"
```
```jldoctest
julia> tolatex(124.0)
L"124.0"
```
"""
function tolatex(x::Float64)
    LaTeXString(replace(string(round(x, sigdigits=6)), r"e(-?\d+)" => s"\\times 10^{\1}"))
end

"""
    register_click_interaction!(ax::Axis, t::Observable)
    
Register the click interaction for the plots: when the user clicks on the plot `ax`,
the Observable `t` will be updated with the time where the user clicked.
"""
function register_click_interaction!(ax::Axis, t::Observable)
    register_interaction!(ax, :click_on_plot) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick
            t[] = round(event.data[1])
        end
    end
end

"""
    path_animation(Γ::Coefficients; [period = 12.0, nsteps = 100])

Create an animation of the path of the N bodies given the coefficients `Γ`. The animation will
show the path of the bodies, the energy of the system, the relative distance between the bodies,
the action, the norm of the gradient of the action, and the dispersion of the energy. The animation
will have a period of `period` seconds and path will be composed by `nsteps` steps.
"""
function MinPath.path_animation(Γ::MinPath.Coefficients; period=12.0, nsteps=100)
    # Calculate the path from the Fourier coefficients
    path = extend_to_period(build_path(Γ, nsteps))
    l = length(path) - 1
    N = length(path[end])
    dim = length(path[end][end])
    # Calculate the fps of the animation
    fps = length(path) / period

    # Compute the needed values for the plots
    action_value = tolatex(action(Γ))
    ∇action_norm = tolatex(norm(∇action(Γ)))
    E = extend_to_period(K_energy(Γ, nsteps) - U(Γ, nsteps))
    meanE = sum(E) / length(E)

    # Set the theme of the window that will contain
    theme = merge(theme_dark(), theme_latexfonts())
    set_theme!(theme)

    # Define the window
    f = Figure(size=(1500, 800))

    # ==== LEFT SIDE OF THE WINDOW ====S
    # On the left, we will have the animation
    # Start by defining the lights and the scene
    al = AmbientLight(RGBf(0.3, 0.3, 0.3))
    dl = DirectionalLight(RGBf(0.9, 0.9, 0.9), Vec3f(-1, 1, -1))
    scene = LScene(f[1:5, 1], show_axis=false, scenekw=(lights=[al, dl],
        backgroundcolor=:black, clear=true))

    # Plot the path
    for i ∈ 1:N
        lines!(scene,
            [[path[j][i][h] for j ∈ 0:l] for h ∈ 1:dim]...,
            color=i, colormap=:lightrainbow, colorrange=(1, N),
            label="Body $i", transparency=true, ssao=true, fxaa=true, linewidth=2)
    end

    # Define the Observables that will contain the position of the N bodies at a given time
    # If the dimension is 2, we will have 2D points, if it is 3, we will have 3D points
    # otherwise we will throw an error
    if dim == 2
        bodies = Observable([Point2(zeros(2)) for _ ∈ 1:N])
    elseif dim == 3
        bodies = Observable([Point3(zeros(3)) for _ ∈ 1:N])
    else
        error("Only 2D and 3D animations are supported")
    end

    # Plot the bodies. Since we are plotting Observables, this plot will be automatically updated
    # when the value of the Observable changes
    meshscatter!(scene, bodies, color=1:N, colormap=:lightrainbow, colorrange=(1, N), markersize=0.05)

    # Add the legend to the plot
    axislegend(scene, labelcolor=:white, labelsize=18, padding=(6.0f0, 12.0f0, 6.0f0, 6.0f0))


    # ==== RIGHT SIDE OF THE WINDOW ====
    # On the top right, we will have a grid containing the value of the action,
    # the norm of the gradient of the action and the disperion of the energy
    grid = GridLayout()
    f.layout[1, 2] = grid
    colsize!(f.layout, 2, Relative(1 / 3))
    padding = (8.0f0, 8.0f0, 4.0f0, 12.0f0)
    Label(grid[1, 1], "Action:", fontsize=24, color=:white, halign=:left, padding=padding)
    Label(grid[1, 2], action_value, halign=:left, fontsize=24, color=:white, padding=padding)
    Label(grid[2, 1], L"||\nabla \text{Action}||:", fontsize=24, halign=:left, color=:white, padding=padding)
    Label(grid[2, 2], ∇action_norm, halign=:left, fontsize=24, color=:white, padding=padding)
    Label(grid[3, 1], L"\Delta E:", fontsize=24, halign=:left, color=:white, padding=padding)
    Label(grid[3, 2], tolatex(maximum(E) - minimum(E)), halign=:left, fontsize=24, color=:white, padding=padding)


    # Below the grid, we will have a slider to control the speed of the animation
    # The scale of the slider is logarithmic, so we can have a wide range of speeds
    sliders = SliderGrid(f[2, 2],
        (label="Speed",
            range=log10(1 / 5):0.1:log10(5.1),
            format=x -> L"\times %$(round(10^x, digits=1))",
            startvalue=0,
            color_active_dimmed=RGBf(0.34, 0.70, 1),
            color_active=RGBf(0.15, 0.40, 0.8),
            color_inactive=RGBf(0.34, 0.34, 0.34)
        ),
        height=Relative(1 / 2)
    )
    sliders.labels[1].color[] = RGBf(1, 1, 1)
    sliders.valuelabels[1].color[] = RGBf(1, 1, 1)
    sliders.labels[1].fontsize[] = 20
    sliders.valuelabels[1].fontsize[] = 20

    # On the bottom right, we will have two plots. One for the energy of the system,
    # and the other for the distance between the bodies.

    # Define the axes
    ax_energy = Axis(f[3, 2], title="Energy", titlecolor=:white, xticklabelcolor=:white,
        yticklabelcolor=:white, ylabel=L"E", ylabelcolor=:white, titlesize=22,
        xzoomlock=true, yzoomlock=true, yrectzoom=false, xrectzoom=false)
    ax_rij = Axis(f[4, 2], title="Relative distance", ylabel=L"r_{ij}", xticklabelcolor=:white,
        yticklabelcolor=:white, ylabelcolor=:white, titlecolor=:white, titlesize=22,
        xzoomlock=true, yzoomlock=true, yrectzoom=false, xrectzoom=false)

    # Plot the energy
    lines!(ax_energy, 0:l, OffsetArrays.no_offset_view(E), label="K Energy", linewidth=2)
    # Make sure thay y-scale isn't too narrow
    ylims!(ax_energy, min(E..., 0.8 * meanE, 1.2 * meanE), max(E..., 1.2 * meanE, 0.8 * meanE))

    # Plot the relative distance
    for i ∈ 1:N-1, j ∈ i+1:N
        r_ij = [norm(path[k][i] - path[k][j]) for k ∈ 0:l]
        lines!(ax_rij, 0:l, r_ij, label=L"r_{%$i%$j}")
    end

    # Add the legend
    Legend(f[5, 2], ax_rij, orientation=:horizontal, labelcolor=:white, labelsize=24)

    # Define the time variable that is an observable, so that the plots can be updated
    t = Observable(0)

    # Register the click interaction for the plots: when the user clicks on the plot,
    # the animation will jump to the time where the user clicked
    register_click_interaction!(ax_energy, t)
    register_click_interaction!(ax_rij, t)

    # Draw a vertical line on the plots to indicate the current time
    vlines!(ax_energy, t, color=:red)
    vlines!(ax_rij, t, color=:red)

    # Display the figure
    display(f, title="Orbit Animation")

    # Start the animation
    while true
        # If the user closed the window, stop the animation
        if !events(scene).window_open[]
            return
        end

        # If the end of the path is reached, start again
        if (t[] == length(path))
            t[] = 0
        end

        # Assign the position of the bodies at the time t to the Observable
        bodies[] = path[t[]]

        # Update the time
        t[] += 1

        # Update the fps according to the speed slider
        fps = length(path) / period * 10^(sliders.sliders[1].value[])

        # Wait for the next frame
        sleep(1.0 / fps)
    end

end
end