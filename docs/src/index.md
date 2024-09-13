# Index

## Example
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


```@repl
using SymPath
P = initialize("example.json")
result = find_orbits(P, OneMethod(BFGS()), number_of_orbits=1)
```

