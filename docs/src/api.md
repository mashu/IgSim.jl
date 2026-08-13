# API

Parameter construction is documented under [Parameters](@ref).

```@autodocs
Modules = [IgSim]
Private = false
Order = [:module, :type, :function, :macro, :constant]
Filter = t -> !(nameof(t) in (:train_params, :SimParams, :DomainFlankMix,
                               :GatedRate, :shm, :shm_none, :shm_igm, :shm_igg,
                               :IlluminaError, :illumina_none, :illumina_miseq,
                               :default_flank_5p, :default_flank_3p))
```
