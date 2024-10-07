# Find new orbits

To find new orbits, first initialize a problem using the `initialize` function (see [Problem definition](problem_definition.md)). Then, use the `find_orbits` function to find new orbits. 

The `find_orbits` function receives the following arguments:
- `P`: A `Problem` object.
- `method`: An optimization method. It must be a subtype of `AbstractMethod`. See [below](#methods) for more information.
- `number_of_orbits`: The number of orbits to be found. Default is `1`. If `number_of_orbits = Inf`, the function will continue to search for orbits until it is stopped by the user.
-  `starting_path_type`: a `Symbol` that indicates the type of the starting path. Default is `:random`. See [below](#starting-path-type) for more information.
- `starting_path`: a `Path` object that indicates the starting path. Default is `nothing`. If `starting_path` is provided, `starting_path_type` is ignored.
- `print_path`: a `Bool` that indicates if the path should be printed. Default is `true`. The filename is the value of the action.
- `print_path_folder`: a `String` that indicates the folder where the path should be printed. Default is `"./"` (the current folder).

Any additional keyword arguments are passed to the optimization method.
For example, when using a `CompoundMethod`, you can pass the keyword `max_iterations` in order to specify the maximum number of times the first and second methods are executed.

## Methods

To specify an optimization method, you can combine multiple `AtomicMethod`s. The available atomic methods are:
- **First Order Methods**: 
    - `BFGS()`: BFGS method.
    - `ConjugateGradient()`: Conjugate gradient method.
    - `GradientDescent()`: Gradient descent method.
- **Second Order Methods**:
    - `Newton()`: Newton method.    

Every atomic method accepts the following options:
- `max_iter`: Maximum number of iterations. Default is `200`.
- `tolerance`: Tolerance on the gradient norm for the stopping criterion. Default is `1e-8`.


You can combine them using the following constructors:
- `OneMethod(method)`: Use a single method.
- `InitMinimizeMethod(init, method)`: Use an initialization method and then a minimization method.
- `CompoundMethod(init, first, second)`: Use a compound method. The `init` method is executed first. Then the `first` and `second` methods are executed repetedly until the stopping criterion is met.

The following are valid examples of method combinations:
- `OneMethod(BFGS())`: Use the BFGS method with 200 iterations and a tolerance of `1e-8` (default options)
- `OneMethod(Newton(400))`: Use the Newton method with 400 iterations (custom) and a tolerance of `1e-8` (default)
- `InitMinimizeMethod(ConjugateGradient(10), BFGS())`: Use the Conjugate Gradient method with 10 iterations (custom) and then the BFGS method with 200 iterations and a tolerance of `1e-8` (default)
- `CompoundMethod(BFGS(5), Newton(300), BFGS(10))`: Use the BFGS method with 5 iterations (custom) to initialize, then the Newton method with 300 iterations (custom) and finally the BFGS method with 10 iterations and a tolerance of `1e-8` (default). The last two methods are executed repeatedly until the stopping criterion is met.

## Starting Path Type

The `starting_path_type` argument can be one of the following:
- `:random`: A random path.
- `:circular`: A circular path.
- `:perturbed_circular`: A circular path with a small perturbation.

## Example

```@repl
using SymPath; P = initialize("example.toml") # hide
result = find_orbits(P, OneMethod(BFGS()), number_of_orbits=1)
```

