#!/usr/bin/python
# -*- encoding: utf8 -*-

import os

script='test_basic'
#os.system('mpirun -np 8  python  '+script+'.py  -mf -nx 100 -ny 100 -nz 4 -kplt 2')
#os.system('mpirun-openmpi-gcc49 -np 8  python2.7  '+script+'.py  -mf -nx 10 -ny 10 -nz 4')
print 'export PYTHONPATH=\"$PYTHONPATH:/Users/aponte/Current_projects/copernicus/dimup/work/qgsolver/\";'
#print('python2.7 setup.py')
print('python2.7 setup.py build_ext --inplace')
print('mpirun-openmpi-gcc49 -np 4  python2.7  ./dev/'+script+'.py  -mf -ksp_view -ksp_monitor -ksp_converged_reason')

print 'All done'


# PETSc options possible 
# -ksp_view -ksp_monitor -ksp_converged_reason
# -ksp_type gmres
# -ksp_type cg

# make movie
#os.system("ffmpeg -y -r 10 -i figs/psiq_%04d.jpg  figs/psiq.mp4")
