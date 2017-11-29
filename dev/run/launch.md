
## on linux and macos workstation

```
bash
source activate petsc
export PYTHONPATH="$PYTHONPATH:path_to_qgsolver"
cd path_to_qgsolver
python setup.py build_ext --inplace
cd dev/run/
mpirun -np 8 python  test_analytical.py -ksp_view -ksp_monitor -ksp_converged_reason
```

## on daparmor
```
cd path_to_qgsolver/dev/run
python run_caparmor.py qg_uniform uniform
```

## PETSc options possible 

```
 -ksp_view -ksp_monitor -ksp_converged_reason
 -ksp_type gmres
 -ksp_type cg
```

