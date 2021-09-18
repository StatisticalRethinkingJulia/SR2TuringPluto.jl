
### Tag version notes

1. git commit -m "Release v4.0.3"
2. git tag v4.0.3
3. git push origin master --tags

### Cloning the repository

```
# Cd to where you would like to clone to
$ git clone https://github.com/StatisticalRethinkingJulia/StatisticalRethinkingTuring.jl
$ cd StatisticalRethinkingTuring.jl
$ julia
```
and in the Julia REPL:

```
julia> ]                                        # Actvate Pkg mode
(@v1.5) pkg> activate .                         # Activate pkg in .
(StatisticalRethinkingTuring) pkg> instantiate    # Install in pkg environment
(StatisticalRethinkingTuring) pkg> <delete>       # Exit package mode
julia>
```

If above procedure fails, if present, try to delete the Manifest.toml file and repeat above steps. As mentioned above, these steps are only needed the first time.

If you want to use a specific tagged version, use:
```
# cd to cloned directory, switch to e.g. v2.0.0:
$ git checkout v2.0.0
```
