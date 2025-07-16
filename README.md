# SymOrb.jl

Symmetric periodic orbits in the n-body problem in Julia

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://dipmathunito.github.io/SymOrb.jl/dev/)

## What

`SymOrb.jl` is a Julia package that provides tools for finding periodic and symmetric solutions to the gravitational $n$-body problem using a variational approach. It is a modular refactoring in a unified Julia environment of the original [symorb](https://github.com/dlfer/symorb) package, which was written in Fortran and Python.

## Documentation 

The online documentation is available at [https://dipmathunito.github.io/SymOrb.jl/dev/](https://dipmathunito.github.io/SymOrb.jl/dev/)

## Installation

A full guide through the installation is available at [https://dipmathunito.github.io/SymOrb.jl/dev/installation/](https://dipmathunito.github.io/SymOrb.jl/dev/installation/)

## Quick start
- Install `SymOrb` follwing the guide linked above
- Write a configuration file `config.toml` 
```toml
symmetry_name = "2d_cyclic_2"
NOB = 3
dim = 2
m = [1, 1, 1]

# Group generators
kern = "TrivialKerTau(2)"
action_type = 0
rotV = "[[-1, 0], [0, -1] ]"
rotS = "(2,3)"

# Other configs
F = 24
Omega = [
    [0, 0],
    [0, 0]
]
```
- Initialize the problem
```julia
using SymOrb
P = initialize("config.toml")
```
- Find a periodic orbit
```julia
orbits = find_orbits(P)
```
- Show the orbit
```julia
using GLMakie
path_animation(P, orbits[1].fourier_coeff)
```