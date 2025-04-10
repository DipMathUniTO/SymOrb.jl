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
2) Choose a symmetry. You can either create your own `.toml` file or choose one example file. For more information about how to define symmetries, see [here](https://dipmathunito.github.io/SymOrb.jl/dev/usage/problem_definition/) 

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
    method =  CompoundMethod(GradientDescent(10), BFGS(5), Newton(200))
    ``` 
5) Start the optimization process using the function `find_orbits`. You have to provide a Problem `P` and a method `method`.
    ```julia
    orbits = find_orbits(P, method)
    ```
    By default, it will start from a random initial condition and will look for a single critical point. It is possible to customize these choices and many others playing with the parameters, as exposed [here](https://dipmathunito.github.io/SymOrb.jl/dev/usage/new_orbit/)
6) Show an animation of the found orbit
    ```julia
    using GLMakie
    path_animation(P, orbits[1].fourier_coeff)
    ```

## Using the database
1. Install MongoDB following the instructions that can be found [here](https://www.mongodb.com/docs/manual/administration/install-community/)
2. Install a database explorer like Compass [here](https://www.mongodb.com/docs/compass/current/install/).
3.  Make sure that the MongoDB server is running (the installation guide has got a section about how to start the MongoDB server locally).
4. In the Julia REPL, load the SymOrb package and connect to the db
    ```julia
    db = connect_to_database("<dbname>")
    ```
    Note that if the database doesn't exist yet, it will be created.
### Saving orbits to the database
Once a database is loaded into `db`, you can instruct `SymOrb.jl` to save new orbits  to the database by setting the `db` parameter of `find_orbits`.
```julia
find_orbits(P, method, db=db)
```

### Saving symmetries to the database
The database not only will contain the trajectories, but also the symmetry used.
There are 3 different ways of saving symmetries to the database:
1. **Automatically**: when a new orbit is saved to the database (see previous step), also the symmetry is automatically saved in the collection "symmetries".
2. **Manually**: it is possible to save a `Problem` object `P` to the database using the function 
    ```julia
        P = save_symmetry_to_database(db, P)
    ```
    It is important to update `P` with the result of the saving operation because it will contain the id of the symmetry in the database, thus allowing to link new orbits to that specific entry, without creating duplicates.
3. **Via Compass**. It is possible to add symmetries to the database using Compass importing a JSON file with many symmetries. Be careful to rename the collection to "symmetries", otherwise it won't be read by SymOrb.

### Loading orbits and symmetries from the database
It is possible to load symmetries from the database using the same `initialize` function used for the TOML files. The available methods are the following:
- `initialize(db, id)`: Loads the symmetry with given id, that is an integer automatically computed when saving symmetries.
- `initialize(db, symmetry_name)`: Loads the symmetry with given name. The name is a string that can be provided in the TOML file.
- `initialize(db, dim=dim, NOB=NOB, action_type=action_type)`: Shows a menu with all the available symmetries with respect to the given dimension, number of bodies and action type. You can navigate through the menu with the arrow keys of your keyboard and select a symmetry with the enter key. Note that `dim, NOB, action_type` are optional. You can provide all of them, only some or none. The order is not important.

To load orbits from the db, you can use 
```julia
load_path(db, dim=dim, NOB=NOB, action_type=action_type)
```
You will be asked first to choose a symmetry (as above) and then all the orbits with such symmetry will be listed. Then you can choose one of these orbits, if any.
