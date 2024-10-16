# Definition of the problem

The problem is defined in a TOML file. The file must contain the following fields:
 - `symmetry_name`: Name of the symmetry group.
 - `NOB`: Number of bodies.
 - `m`: An array of masses of the bodies.
 - `dim`: The dimension of the space.
 - *Group generators*
    - `action_type`: The type of the action of the group. It can be 
        - 0: Circular action
        - 1: Dihedral action
        - 2: Brake action
    - `kern`: A string containing a GAP command to generate the kernel of the $\tau$ representation. The following helper functions ara available: 
        - `TrivialKerTau(2)` for the trivial kernel of the $\tau$ representation when $d=2$
        - `TrivialKerTau(3)` for the trivial kernel of the $\tau$ representation when $d=3$
    - `rotV`: A string containing a GAP command to generate the $O(d)$ part of the first generator of the group. It must be an orthogonal matrix of size $d \times d$
    - `rotS`: A string containing a GAP command to generate the $\Sigma_n$ part of the first generator of the group. It must be a permutation.
    - If the action is dihedral or brake:
        - `refV`: A string containing a GAP command to generate the $O(d)$ part of the second generator of the group. It must be an orthogonal matrix of size $d \times d$
        - `refS`: A string containing a GAP command to generate the $\Sigma_n$ part of the second generator of the group. It must be a permutation.

The file must also contain the following optional fields:
 - `F`: The number of fourier coefficients (default: `24`)
 - `steps`: The number of steps to be used for the point representation of the path (default: `2 * F`)
 - `Omega`: The infinitesimal generator of the reference frame rotation (default: `0`)
 - `denominator`: The denominator of the potential (default: `x`)


The following is an example of a TOML file defining a problem:
```TOML
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

## Load the problem 
A problem can be loaded using the `initialize` function. The function receives the path to the TOML file and returns a `Problem` object.

```@repl 
using SymOrb # hide
P = initialize("example.toml")
```




