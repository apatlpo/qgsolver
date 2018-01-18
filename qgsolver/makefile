# -*- makefile -*-

# make PYTHON=/Users/aponte/.miniconda3/envs/petsc/bin/python build

MPIEXEC=
#PYTHON := /Users/aponte/.miniconda3/envs/petsc/bin/python
#SHELL := /bin/bash

.PHONY:test
test: run clean

#SCRIPT=run_demo
MODULE=timestepper_adv

.PHONY:build
build: ${MODULE}.so

#.PHONY:run
#run: build
#	${PYTHON} ${SCRIPT}.py

${MODULE}.so: ${MODULE}.pyf ${MODULE}.F90
	@echo $$PATH
	@echo $(PYTHON)
	@echo $(SHELL)
	env \
        F77="${FC}" F77FLAGS="${FC_FLAGS}" \
        F90="${FC}" F90FLAGS="${FC_FLAGS}" \
        LDSHARED="${FC_LINKER}" \
	${PYTHON} setup.py -q build_ext --inplace
	${RM} -r build ${MODULE}module.c

.PHONY:clean
clean::
	${RM} ${MODULE}*.so
	${RM} *.py[co]
	${RM} -r __pycache__

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
MPIEXEC=
