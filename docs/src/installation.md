# Installation

## Requirements
- Julia `1.10` or later
- OS: Linux or macOS. To run `SymOrb.jl` on Windows, WSL is required (see Step 0 below) 
- [Optional] An OpenGL compatible GPU (for animations). All latest GPUs should work out of the box.
- [Recommended] 8 GB of RAM memory


!!! warning 
    `SymOrb` uses, to interface to the GAP package `symorb`, the package `GAP.jl` that is not supported on Windows. However, you can still use `SymOrb` on Windows by using Windows Subsystem for Linux (WSL). [^1]

### Step 0: Install and configure WSL (Windows only)
!!! info
    A complete guide WSL and its installation can be found [here](https://learn.microsoft.com/en-us/windows/wsl/install). 

Type the  following command in the `Powershell` to install WSL and Ubuntu.
```bash
$ wsl --install
```

Then, Open a `Bash` shell (for instance, opening the `Ubuntu` app from the Start menu). Follow the instructions to set up your user account by choosing a username and password.

Finally, to enable orbits visualisation, you need to install OpenGL [^2]

If you run Windows 11, it is sufficient to install `mesa-utils` by typing in the `Bash` shell:

```bash
$ sudo apt install mesa-utils
```

### Step 1: Install Julia (all platforms)
According to the Julia official documentation, it is strongly advised to install Julia using `juliaup` instead of using the distribution's package manager and binaries (see [here](https://github.com/JuliaLang/juliaup))

Open a shell and type the following command to install `juliaup`:
```bash
$ curl -fsSL https://install.julialang.org | sh -s -- -y
``` 

### Step 2: Install `SymOrb` (all platforms)
!!! info
    `SymOrb.jl` is not yet registered in the Julia General registry. Therefore, you need to install it directly from the GitHub repository.

In the following, it is assumed that you want to create a new project where you will use `SymOrb`. If you want to install `SymOrb` in an existing project, you can skip the first two steps.


Julia projects are simply directories that contain a `Project.toml` file. 

Hence, first of all, create a new directory. If you are using Windows, make sure to create it in the WSL Linux file system. (The Linux file system can be accessed through `File Explorer` as a directory listed in the left sidebar. It is advised to create your documents in the `home` directory located at `/home/your_username`.)

Then, to generate the required `Project.toml` and `Manifest.toml` files, open a `Bash` shell in that directory and launch the Julia REPL by typing `julia --project` and pressing `Enter`.

To install `SymOrb`, you need to access the Julia package manager by typing `]`. You will notice that the prompt changes to `(v1.xx) pkg>`. Then, type the following command and press `Enter`:
```julia
(v1.xx) pkg> add https://github.com/susorb/SymOrb.jl
```

To exit the package manager, simply press `Backspace` to return to the Julia REPL.

To check that `SymOrb` has been correctly installed, type the following command in the Julia REPL:
```@repl
using SymOrb
```

If no error message is displayed, then `SymOrb` has been correctly installed.


[^1]: [https://learn.microsoft.com/en-us/windows/wsl/install]
[^2]: https://gist.github.com/Mluckydwyer/8df7782b1a6a040e5d01305222149f3c
