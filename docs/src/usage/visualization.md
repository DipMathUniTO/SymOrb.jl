# Loading and visualizing paths
## Loading paths
To load paths saved by the `find_orbits` function, use the `read_path_from_file` function. This function receives a single argument, `filename`, which is a `String` that indicates the path to the file where the path has been saved by the `find_orbits` function.

Note that it outputs a `Problem` object and a `Path` object.

### Example

```@repl
using SymOrb # hide
P, Γ = read_path_from_file("7.4797.toml") 
```

## Orbits visualization

To visualize orbits, use the `path_animation` function. 

!!! note
    The `path_animation` function uses the `GLMakie.jl` package to create the animation. Therefore, you need to install it by typing the following command in the Package Manager:
    ```julia
    ] add GLMakie
    ```
    Note that this package must be explicitly loaded before using the `path_animation` function.

The `path_animation` function receives the following arguments:
-`P`: A `Problem` object.
- `Γ`: A `Path` object.
or 
-`filename`: A `String` that indicates the path to the file where the path has been saved by the `find_orbits` function.

### Example

```julia
path_animation(P, Γ)
```

Or, if you want to load the path directly from the file:

```julia
path_animation("7.4797.toml")
```