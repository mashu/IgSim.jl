# API

Parameter construction is documented under [Parameters](@ref).

```@autodocs
Modules = [IgSim]
Private = false
Order = [:module, :type, :function, :macro, :constant]
Filter = t -> !(nameof(t) in (:train_params, :SimParams, :DomainFlankMix))
```
