# Index

## Example

```@repl
using SymPath
P = initialize("example.toml")
result = find_orbits(P, OneMethod(BFGS()), number_of_orbits=1)
```

