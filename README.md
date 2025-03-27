# SymOrb.jl

Symmetric periodic orbits in the n-body problem in Julia

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://dipmathunito.github.io/SymOrb.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dipmathunito.github.io/SymOrb.jl/dev/)

## What

`SymOrb.jl` is a Julia package provides tools for finding periodic and symmetric solutions to the gravitational $n$-body problem using a variational approach. It is a reimplementation in a unified Julia environment of the original [symorb](https://github.com/dlfer/symorb) package, which was written in Fortran and Python.

## Installation of the dev version
> [!IMPORTANT]
> `SymOrb.jl` uses a package, `GAP.jl`, that is working only on Linux and Mac OS. If you are using Windows, you can still use `SymOrb.jl` in WSL. See Step 0 of the [public guide](https://dipmathunito.github.io/SymOrb.jl/dev/installation/).

### On your Computer
1) Install `git` and log in with your account 
4) Install Julia (see Step 1 of the [public guide](https://dipmathunito.github.io/SymOrb.jl/dev/installation/).)
3) Open a shell. Create a `symorb` directory in your home with the command `mkdir symorb` and enter it with `cd symorb`.
2) Clone this repository doing either:
    - `git clone https://github.com/DipMathUniTO/SymOrb.jl-dev.git`
    - Download the zip file from the "Code" green button and unzip it in the `symorb` directory.
3) To access the `dev` version of the package, run the `git checkout devs` command.

6) Open Julia with the `julia` command.
7) In the Julia shell, type `]` to open the package manager. You will see that the green `julia>` prompt changes to a blue one. Note that the text in round brackets is the Julia environment you are using. 
8) (optional) If you don't want to work in your global Julia environment, create a new one with the command `generate symorb_env` in the  Julia Package Manager. Every other time you want to run Julia using this environment, you have to add a `--project=symorb_env` flag after the `julia` command.
9) Add `SymOrb` to your Julia environment (either the global one or your `symorb_env`) by entering the Package Manager with  `]` (see step 7) and the command `dev SymOrb.jl-dev`. Julia will start downloading all the dependences and pre-compiling them. This can take a while.
10) Install other required packages with `add GLMakie`. This will also take a while.
11) Exit the package manager with the <kbd>Backspace</kbd>  key.
12) To exit Julia, use <kbd>Ctrl</kbd> + <kbd>D</kbd>.

 
### On HPC:
Repeat the steps 2-9 of the previous guide. DO not install `GLMakie` unless you have enabled `Xorg` `ssh` forwarding.

## Quick start
> [!NOTE]
> We are assuming that you installed `SymOrb.jl` using the previous guide. If you installed it to a different location, modify the paths accordingly.

1) Open a shell and navigate to `~/symorb` with `cd symorb`.
2) (JUST ONCE) If you want some example symmetries, they are located in the `exaamples` branch. To open this branch to a different directory, navigate to `SymOrb.jl-dev` and run `git worktree add ../_examples examples`. Now you should have an additional `_examples` directory in your `symorb` one. To navigate back to `symorb`, use `cd ..`.
3) Open julia with `julia` (if you created a separate environment, use `julia --project=symorb_env`).
4) Import `SymOrb` with `using SymOrb`. If it's the first time, it will take some time.
5) If there are some import errors, it might be because you changed git branch or updated it. To solve that, go to `SymOrb.jl-dev`, open `julia`, run `] resolve` and `] instantiate`, go back to `symorb` and do the same. 

## Find new orbits
1) Follow the Quick start guide to have a woking `SymOrb` environment.
2) Choose a symmetry. You can either create your own `.toml` file or choose one example file. 
3) Load the symmetry with 
```julia
P = initialize("_examples/orbits/2d_cyclic.toml")
```
4) Choose an optimization method. There are 4 different "atomic" methods:
    - `GradientDescent(`*`options`*`)` (1st order steepest descent)
    - `ConjugateGradient(`*`options`*`)` (1st order steepest descent)
    - `BFGS(`*`options`*`)` (1st order steepest descent)
    - `Newton(`*`options`*`)` (2nd order)
    
    You can choose to use only a single method or to combine them.
    - To use only one method, for instance BFGS, do 
    ```julia 
    method = OneMethod(BFGS())
    ```
    - You can combine an initialization method and another method using, for instance, 
    ```julia
    method = InitMinimizeMethod(GradientDescent(10), Newton(200))
    ```
    - You can also combine up to three methods. The first one is an initialization method that runs only once. Then the other two methods are executed by alternating between them. For instance, you can do 
    ```julia 
    method = CompoundMethod(GradientDescent(10), BFGS(5), Newton(200))
    ``` 
    
- Install `SymOrb` by typing `] add https://github.com/DipMathUniTO/SymOrb.jl`
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