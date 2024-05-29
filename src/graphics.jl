
tolatex(x::Float64) = LaTeXString(replace(string(round(x, sigdigits=6)), r"e(-?\d+)" => s"\\times 10^{\1}") )


function show(Γ::Coefficients; period = 12.0, nsteps = 100)
    
    path = reconstruct_path(build_path(Γ, nsteps))
    
    theme = merge(theme_dark(), theme_latexfonts())
    set_theme!(theme)

    action_value = tolatex(action(flatten(Γ)))
    ∇action_norm = tolatex(norm(∇action(flatten(Γ))))

    f = Figure(size=(1500, 800), color=:black, title = "My Custom Window Title")
    al = AmbientLight(RGBf(0.3, 0.3, 0.3))
    dl = DirectionalLight(RGBf(0.9, 0.9, 0.9), Vec3f(-1, 1, -1))
    scene = LScene(f[1:5, 1], show_axis=false, scenekw = (lights = [al, dl], backgroundcolor=:black, clear=true))

    Box(f[1,2], color=:transparent, strokevisible=false)
    
    grid = GridLayout()
    f.layout[1,2] = grid

    padding_left = (8.0f0, 8.0f0, 4.0f0, 12.0f0)
    padding_right = (8.0f0, 8.0f0, 4.0f0, 12.0f0)

    Label(grid[1, 1], "Action:", fontsize=24, color=:white,halign=:left, padding=padding_left)
    Label(grid[1, 2],  action_value, halign=:left, fontsize=24, color=:white, padding=padding_right)


    Label(grid[2, 1], L"||\nabla \text{Action}||:", fontsize=24, halign=:left, color=:white, padding=padding_left)
    Label(grid[2, 2],  ∇action_norm, halign=:left, fontsize=24, color=:white, padding=padding_right)

    colsize!(f.layout, 2, Relative(1/3))
    ax = Axis(f[4,2], title="Distance between bodies", titlecolor=:white, titlesize=22, xzoomlock=true, yzoomlock=true, yrectzoom=false,xrectzoom=false)

    l = length(path)-1

    for i ∈ 1:N
        lines!(scene, 
                [[path[j][i][h] for j ∈ 0:l] for h ∈ 1:dim]..., 
                color=i, colormap=:lightrainbow, colorrange = (1,N), label="Body $i", transparency=true, ssao=true, fxaa=true, linewidth=2);
    end

    for i ∈ 1:N-1, j ∈ i+1:N
        r_ij = [norm(path[k][i] - path[k][j]) for k ∈ 0:l]
        lines!(ax, 0:(length(r_ij)-1), r_ij, label=L"r_{%$i%$j}")
    end

    Legend(f[5,2], ax,  orientation = :horizontal, labelcolor=:white, labelsize=24)
    axislegend(scene,  labelcolor=:white, labelsize=18, padding=(6.0f0, 12.0f0, 6.0f0, 6.0f0))
    
    fps = length(path)/period

    if dim == 2
        bodies = Observable([Point2(zeros(2)) for _ ∈ 1:N])
    elseif dim == 3
        bodies = Observable([Point3(zeros(3)) for _ ∈ 1:N])
    else 
        return
    end


    meshscatter!(scene, bodies, color=1:N, colormap=:lightrainbow, colorrange=(1,N))
    t = 0
    t = Observable(0)
    vlines!(ax, t, color=:red)


    sliders = SliderGrid(f[2, 2], 
        (label = "Speed", 
         range = log10(1/5):0.1:log10(5.1), 
         format = x -> L"\times %$(round(10^x, digits=1))", 
         startvalue = 0, 
         color_active_dimmed = RGBf(0.34, 0.70, 1), 
         color_active=RGBf(0.15, 0.40, 0.8), 
         color_inactive=RGBf(0.34, 0.34, 0.34)
        ),
        height=Relative(1/2)
    )
    
    sliders.labels[1].color[] = RGBf(1, 1, 1)
    sliders.valuelabels[1].color[] = RGBf(1, 1, 1)
    sliders.labels[1].fontsize[] = 20
    sliders.valuelabels[1].fontsize[] = 20
       
    display(f,  title = "Orbit Animation")

    i = Observable(0)
    register_interaction!(ax, :my_interaction) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick
            t[] = round(event.data[1])
        end
    end
    while true
        if ! events(scene).window_open[] return end
        if(t[] == length(path)) t[] = 0 end
        bodies[] = path[t[]]
        t[] += 1
        fps = length(path)/period*10^(sliders.sliders[1].value[])
        sleep(1.0/fps)
    end


end

