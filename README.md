# Identifying the Smallest Adversarial Load Perturbations that Render DC-OPF Infeasible
<p align="center">
<img width="300" alt="image" src="https://github.com/user-attachments/assets/4d6104b7-bebe-40c4-afc3-b59d95a67239" />
</p>

To set up and initialize this respository locally, clone the folder, launch Julia locally, and then run 
```
julia> ]
(@v1.11) pkg> activate .
(DCAttack) pkg> instantiate
```

* ```run_tests.jl``` will re-run the attack and defense tests, storing all saved data (via ```HDF5```) into the ```data``` folder.
* ```case_study.jl``` will run BaB on the 5-bus network to re-create the callbacks plot
* ```analyze_data.jl``` will re-create the data in the table (once ```run_tests.jl``` is re-run
* please direct questions to Sam Chevalier: ```schevali[@]uvm.edu```

A [Gurobi](https://www.gurobi.com/) license is needed to run BaB.
