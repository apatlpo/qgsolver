/*  Fast C/OpenMP interpolator for tvs_to_s_fast()
    Runs about 21x faster than single processor tvs_to_s
    Usage: vt = fast_interp(zt,ssh,vs,H)
    Where: zt  is a 1D numpy array, size NZT      , type float64
           z   is a 3D numpy array, size NZxNYxNX, type float64
           vs  is a 3D numpy array, size NZxNYxNX, type float32 or type float64
    It will almost certainly not work without modification
    if the arguments are not exactly like above.
    See tvs_to_s_fast() in utils.py for it's intended use.      */

#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include <stdio.h>
#include <omp.h>

// Linear interpolate/extrapolate function
static inline double lin_int(double x0, double y0, double x1, double y1, double x)
{
   return y0 + (y1-y0)*(x-x0)/(x1-x0);
}

// Linear interpolate a 1D vector, assumes everything is sorted increasing
static inline void myinterp(double *zsk, double *vsk, double *ztk, double *vtk, npy_uint NZ, npy_uint NZT)
{
    npy_uint k, k0, k1 = 1; // initial guess for k1
    for (k=0;k<NZT;k++)
    {
      //if (ztk[k] <= zsk[0])         vtk[k] = lin_int(zsk[0]   ,vsk[0]   ,zsk[1]   ,vsk[1]   ,ztk[k]);
      //else if (ztk[k] >= zsk[NZ-1]) vtk[k] = lin_int(zsk[NZ-2],vsk[NZ-2],zsk[NZ-1],vsk[NZ-1],ztk[k]);
      //if (ztk[k] <= zsk[0])         vtk[k] = vsk[0];
      if (ztk[k] <= zsk[0])         vtk[k] = NPY_NAN;
      else if (ztk[k] >= zsk[NZ-1]) vtk[k] = vsk[NZ-1];
      else {
          while (ztk[k] > zsk[k1]) k1++;
          k0 = k1-1;
          vtk[k] = lin_int(zsk[k0],vsk[k0],zsk[k1],vsk[k1],ztk[k]);
      }
    }
}
// Main fast_interp function
static PyObject* fast_interp3D(PyObject* self, PyObject* args)
{
    PyArrayObject *zt, *z, *vs;
    PyArrayObject *vt=NULL;
    double        *zsk, *vsk, *ztk, *vtk;
    npy_uint      NX,NY,NZ,NZT,nd,i,j,k,typ;

    // Parse numpy ndarray arguments
    if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &zt,
        &PyArray_Type, &z, &PyArray_Type, &vs ))
        return NULL;

    // Check meta data
    nd=PyArray_NDIM(vs);
    if (nd!=3) {
        printf("Input array vs must have three dimensions\n");
        return NULL;
    }
    typ = PyArray_TYPE(vs);

    // Get input dimensions
    NZT = PyArray_DIM(zt,0);
    NZ = PyArray_DIM(vs,0);
    NY = PyArray_DIM(vs,1);
    NX = PyArray_DIM(vs,2);

    // Output array size
    npy_intp dims[] = {NZT,NY,NX};
    vt = (PyArrayObject*)PyArray_SimpleNew(nd, dims, NPY_FLOAT64);
    if (vt == NULL) return NULL;

// These macros make the loops below much cleaner
#define  ZT(k)        ( *((double*)PyArray_GETPTR1(zt,k    )) )
#define  Z(k, j, i)   ( *((double*)PyArray_GETPTR3(z,k,j,i) ) )
#define  VT(k, j, i)  ( *((double*)PyArray_GETPTR3(vt,k,j,i)) )
#define  DVS(k, j, i) ( *((double*)PyArray_GETPTR3(vs,k,j,i)) )
#define  FVS(k, j, i) ( *(( float*)PyArray_GETPTR3(vs,k,j,i)) )

#pragma omp parallel shared(zt,z,vs,vt,NX,NY,NZ,NZT,nd,typ) private(i,j,k,zsk,vsk,ztk,vtk)
    {
      // Thread-local constant and work arrays
      //invH = 1.0/DEPTH(0);
      zsk = malloc(sizeof(double)*NZ);
      vsk = malloc(sizeof(double)*NZ);
      ztk = malloc(sizeof(double)*NZT); for (k=0;k<NZT;k++) ztk[k]=ZT(k);
      vtk = malloc(sizeof(double)*NZT);
#pragma omp for
      for (j=0;j<NY;j++) {
        for (i=0;i<NX;i++) {
          if (typ==NPY_FLOAT32)      for (k=0;k<NZ;k++) vsk[k] = (double)FVS(k,j,i); // data
          else if (typ==NPY_FLOAT64) for (k=0;k<NZ;k++) vsk[k] = DVS(k,j,i); // data
          for (k=0;k<NZ;k++) zsk[k] = Z(k,j,i); // z values
          myinterp(zsk,vsk,ztk,vtk,NZ,NZT); // interpolate this column
          for (k=0;k<NZT;k++) VT(k,j,i) = vtk[k]; // copy to output
        } // end for i
      } // end for j
      free(zsk); free(vsk); free(ztk); free(vtk); // free work arrays
    } // end pragma
    return PyArray_Return(vt);
}

// Define functions in module
static PyMethodDef FastInterp3D[] =
{
    {"interp", fast_interp3D, METH_VARARGS, "fast interp"}, {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION < 3
// Module initialization for Py 2
PyMODINIT_FUNC initfast_interp3D(void)
{
    (void) Py_InitModule("fast_interp3D", FastInterp3D);
    import_array();
}
#else
// Define module for Python 3
static struct PyModuleDef fast_interp_definition = {
    PyModuleDef_HEAD_INIT,
    "fast_interp3D",
    "C implementation of interp.",
    -1,
    FastInterp3D
};

// Module initialization for Python 3
PyMODINIT_FUNC PyInit_fast_interp3D(void)
{
    Py_Initialize();
    PyObject *m = PyModule_Create(&fast_interp_definition);
    import_array();
    return m;
}
#endif


