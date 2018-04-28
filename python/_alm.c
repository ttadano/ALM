#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "alm_wrapper.h"

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * py_alm_new(PyObject *self, PyObject *args);
static PyObject * py_alm_delete(PyObject *self, PyObject *args);
static PyObject * py_run_suggest(PyObject *self, PyObject *args);
static PyObject * py_optimize(PyObject *self, PyObject *args);
static PyObject * py_set_cell(PyObject *self, PyObject *args);
static PyObject * py_set_displacement_and_force(PyObject *self, PyObject *args);
static PyObject * py_get_ndata_used(PyObject *self, PyObject *args);
static PyObject * py_set_constraint_type(PyObject *self, PyObject *args);
static PyObject * py_set_norder(PyObject *self, PyObject *args);
static PyObject * py_set_cutoff_radii(PyObject *self, PyObject *args);
static PyObject * py_generate_force_constant(PyObject *self, PyObject *args);
static PyObject * py_get_atom_mapping_by_pure_translations
(PyObject *self, PyObject *args);
static PyObject * py_get_number_of_displacement_patterns
(PyObject *self, PyObject *args);
static PyObject * py_get_numbers_of_displacements(PyObject *self, PyObject *args);
static PyObject * py_get_displacement_patterns(PyObject *self, PyObject *args);
static PyObject * py_get_number_of_fc_elements(PyObject *self, PyObject *args);
static PyObject * py_get_number_of_irred_fc_elements(PyObject *self, PyObject *args);
static PyObject * py_get_fc_origin(PyObject *self, PyObject *args);
static PyObject * py_get_fc_irreducible(PyObject *self, PyObject *args);
static PyObject * py_get_fc_all(PyObject *self, PyObject *args);
static PyObject * py_set_fc(PyObject *self, PyObject *args);
static PyObject * py_get_matrix_elements(PyObject *self, PyObject *args);

struct module_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
  struct module_state *st = GETSTATE(m);
  PyErr_SetString(st->error, "something bad happened");
  return NULL;
}

static PyMethodDef _alm_methods[] = {
  {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
  {"alm_new", py_alm_new, METH_VARARGS, ""},
  {"alm_delete", py_alm_delete, METH_VARARGS, ""},
  {"run_suggest", py_run_suggest, METH_VARARGS, ""},
  {"optimize", py_optimize, METH_VARARGS, ""},
  {"set_cell", py_set_cell, METH_VARARGS, ""},
  {"set_displacement_and_force", py_set_displacement_and_force, METH_VARARGS, ""},
  {"set_constraint_type", py_set_constraint_type, METH_VARARGS, ""},
  {"set_norder", py_set_norder, METH_VARARGS, ""},
  {"set_cutoff_radii", py_set_cutoff_radii, METH_VARARGS, ""},
  {"generate_force_constant", py_generate_force_constant, METH_VARARGS, ""},
  {"get_atom_mapping_by_pure_translations", py_get_atom_mapping_by_pure_translations,
   METH_VARARGS, ""},
  {"get_number_of_displacement_patterns", py_get_number_of_displacement_patterns,
   METH_VARARGS, ""},
  {"get_numbers_of_displacements", py_get_numbers_of_displacements,
   METH_VARARGS, ""},
  {"get_ndata_used", py_get_ndata_used, METH_VARARGS, ""},
  {"get_displacement_patterns", py_get_displacement_patterns, METH_VARARGS, ""},
  {"get_number_of_fc_elements", py_get_number_of_fc_elements, METH_VARARGS, ""},
  {"get_number_of_irred_fc_elements", py_get_number_of_irred_fc_elements, METH_VARARGS, ""},
  {"get_fc_origin", py_get_fc_origin, METH_VARARGS, ""},
  {"get_fc_irreducible", py_get_fc_irreducible, METH_VARARGS, ""},
  {"get_fc_all", py_get_fc_all, METH_VARARGS, ""},
  {"set_fc", py_set_fc, METH_VARARGS, ""},
  {"get_matrix_elements", py_get_matrix_elements, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int _alm_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int _alm_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_alm",
  NULL,
  sizeof(struct module_state),
  _alm_methods,
  NULL,
  _alm_traverse,
  _alm_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__alm(void)

#else
#define INITERROR return

  void
  init_alm(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  PyObject *module = Py_InitModule("_alm", _alm_methods);
#endif

  if (module == NULL)
    INITERROR;
  struct module_state *st = GETSTATE(module);

  st->error = PyErr_NewException("_alm.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}

static PyObject * py_alm_new(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  const int id = alm_new();

  return PyLong_FromLong((long) id);
}

static PyObject * py_alm_delete(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i",
                              &id)) {
    return NULL;
  }

  alm_delete(id);

  Py_RETURN_NONE;
}

static PyObject * py_run_suggest(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i",
                              &id)) {
    return NULL;
  }

  alm_run_suggest(id);

  Py_RETURN_NONE;
}

static PyObject * py_optimize(PyObject *self, PyObject *args)
{
  int id, info;
  if (!PyArg_ParseTuple(args, "i",
                              &id)) {
    return NULL;
  }

  info = alm_optimize(id);

  return PyLong_FromLong((long) info);
}

static PyObject * py_set_cell(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_lavec;
  PyArrayObject* py_xcoord;
  PyArrayObject* py_kd;
  if (!PyArg_ParseTuple(args, "iOOO",
                              &id,
                              &py_lavec,
                              &py_xcoord,
                              &py_kd)) {
    return NULL;
  }

  double (*lavec)[3] = (double(*)[3])PyArray_DATA(py_lavec);
  double (*xcoord)[3] = (double(*)[3])PyArray_DATA(py_xcoord);
  const int* kd = (int*)PyArray_DATA(py_kd);
  const int nat = PyArray_DIMS(py_kd)[0];

  alm_set_cell(id, nat, lavec, xcoord, kd);

  Py_RETURN_NONE;
}

static PyObject * py_set_displacement_and_force(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_u;
  PyArrayObject* py_f;
  if (!PyArg_ParseTuple(args, "iOO",
                              &id,
                              &py_u,
                              &py_f)) {
    return NULL;
  }

  const double* u = (double*)PyArray_DATA(py_u);
  const double* f = (double*)PyArray_DATA(py_f);

  const int ndata_used = PyArray_DIMS(py_f)[0];
  const int nat = PyArray_DIMS(py_f)[1];
  alm_set_displacement_and_force(id, u, f, nat, ndata_used);

  Py_RETURN_NONE;
}

static PyObject * py_set_constraint_type(PyObject *self, PyObject *args)
{
  int id, iconst;
  if (!PyArg_ParseTuple(args, "ii",
                              &id,
                              &iconst)) {
    return NULL;
  }

  alm_set_constraint_type(id, iconst);

  Py_RETURN_NONE;
}

static PyObject * py_set_norder(PyObject *self, PyObject *args)
{
  int id, norder;
  if (!PyArg_ParseTuple(args, "ii",
                              &id,
                              &norder)) {
    return NULL;
  }

  alm_set_norder(id, norder);

  Py_RETURN_NONE;
}

static PyObject * py_set_cutoff_radii(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_rcs;
  if (!PyArg_ParseTuple(args, "iO",
                              &id,
                              &py_rcs)) {
    return NULL;
  }

  double *rcs = (double(*))PyArray_DATA(py_rcs);

  alm_set_cutoff_radii(id, rcs);

  Py_RETURN_NONE;
}

static PyObject * py_generate_force_constant(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i",
                              &id)) {
    return NULL;
  }

  alm_generate_force_constant(id);

  Py_RETURN_NONE;
}

static PyObject * py_get_atom_mapping_by_pure_translations
(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_map_p2s;

  if (!PyArg_ParseTuple(args, "iO",
                              &id,
                              &py_map_p2s)) {
    return NULL;
  }

  int *map_p2s = (int(*))PyArray_DATA(py_map_p2s);
  const int num_trans = alm_get_atom_mapping_by_pure_translations(id, map_p2s);

  return PyLong_FromLong((long) num_trans);
}

static PyObject * py_get_number_of_displacement_patterns
(PyObject *self, PyObject *args)
{
  int id;
  int fc_order;
  int num_patterns;

  if (!PyArg_ParseTuple(args, "ii",
                              &id,
                              &fc_order)) {
    return NULL;
  }

  num_patterns = alm_get_number_of_displacement_patterns(id, fc_order);

  return PyLong_FromLong((long) num_patterns);
}

static PyObject * py_get_numbers_of_displacements(PyObject *self, PyObject *args)
{
  int id, fc_order;
  PyArrayObject* py_numbers;

  if (!PyArg_ParseTuple(args, "iOi",
                              &id,
                              &py_numbers,
                              &fc_order)) {
    return NULL;
  }

  int (*numbers) = (int(*))PyArray_DATA(py_numbers);

  alm_get_numbers_of_displacements(id, numbers, fc_order);

  Py_RETURN_NONE;
}

static PyObject * py_get_displacement_patterns(PyObject *self, PyObject *args)
{
  int id, fc_order;
  PyArrayObject* py_atom_indices;
  PyArrayObject* py_disp_patterns;

  if (!PyArg_ParseTuple(args, "iOOi",
                              &id,
                              &py_atom_indices,
                              &py_disp_patterns,
                              &fc_order)) {
    return NULL;
  }

  int (*atom_indices) = (int(*))PyArray_DATA(py_atom_indices);
  double (*disp_patterns) = (double(*))PyArray_DATA(py_disp_patterns);

  const int disp_basis = alm_get_displacement_patterns(id,
                                                       atom_indices,
                                                       disp_patterns,
                                                       fc_order);

  return PyLong_FromLong((long) disp_basis);
}

static PyObject * py_get_number_of_fc_elements(PyObject *self, PyObject *args)
{
  int id, fc_order;

  if (!PyArg_ParseTuple(args, "ii",
                              &id,
                              &fc_order)) {
    return NULL;
  }

  const int num_fc_elems = alm_get_number_of_fc_elements(id, fc_order);

  return PyLong_FromLong((long) num_fc_elems);
}

static PyObject * py_get_number_of_irred_fc_elements(PyObject *self, PyObject *args)
{
  int id, fc_order;

  if (!PyArg_ParseTuple(args, "ii",
                              &id,
                              &fc_order)) {
    return NULL;
  }

  const int num_fc_elems = alm_get_number_of_irred_fc_elements(id, fc_order);

  return PyLong_FromLong((long) num_fc_elems);
}

static PyObject * py_get_ndata_used(PyObject *self, PyObject *args)
{
  int id;

  if (!PyArg_ParseTuple(args, "i",
                              &id)) {
    return NULL;
  }
  const int ndata_used = alm_get_ndata_used(id);

  return PyLong_FromLong((long) ndata_used);
}


static PyObject * py_get_fc_origin(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_fc_values;
  PyArrayObject* py_elem_indices;

  if (!PyArg_ParseTuple(args, "iOO",
                              &id,
                              &py_fc_values,
                              &py_elem_indices)) {
    return NULL;
  }

  double (*fc_values) = (double(*))PyArray_DATA(py_fc_values);
  int (*elem_indices) = (int(*))PyArray_DATA(py_elem_indices);
  const int fc_order = PyArray_DIMS(py_elem_indices)[1] - 1;

  alm_get_fc_origin(id, fc_values, elem_indices, fc_order);

  Py_RETURN_NONE;
}

static PyObject * py_get_fc_irreducible(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_fc_values;
  PyArrayObject* py_elem_indices;

  if (!PyArg_ParseTuple(args, "iOO",
                              &id,
                              &py_fc_values,
                              &py_elem_indices)) {
    return NULL;
  }

  double (*fc_values) = (double(*))PyArray_DATA(py_fc_values);
  int (*elem_indices) = (int(*))PyArray_DATA(py_elem_indices);
  const int fc_order = PyArray_DIMS(py_elem_indices)[1] - 1;

  alm_get_fc_irreducible(id, fc_values, elem_indices, fc_order);

  Py_RETURN_NONE;
}

static PyObject * py_get_fc_all(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_fc_values;
  PyArrayObject* py_elem_indices;

  if (!PyArg_ParseTuple(args, "iOO",
                              &id,
                              &py_fc_values,
                              &py_elem_indices)) {
    return NULL;
  }

  double (*fc_values) = (double(*))PyArray_DATA(py_fc_values);
  int (*elem_indices) = (int(*))PyArray_DATA(py_elem_indices);
  const int fc_order = PyArray_DIMS(py_elem_indices)[1] - 1;

  alm_get_fc_all(id, fc_values, elem_indices, fc_order);

  Py_RETURN_NONE;
}

static PyObject * py_set_fc(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_fc_in;

  if (!PyArg_ParseTuple(args, "iO",
                              &id,
                              &py_fc_in)) {
    return NULL;
  }

  double *fc_in = (double(*))PyArray_DATA(py_fc_in);

  alm_set_fc(id, fc_in);

  Py_RETURN_NONE;
}

static PyObject * py_get_matrix_elements(PyObject *self, PyObject *args)
{
  int id;
  int nat;
  int ndata_used;
  PyArrayObject* py_amat;
  PyArrayObject* py_bvec;

  if (!PyArg_ParseTuple(args, "iiiOO",
                              &id,
                              &nat,
                              &ndata_used,
                              &py_amat,
                              &py_bvec)) {
    return NULL;
  }

  double (*amat) = (double(*))PyArray_DATA(py_amat);
  double (*bvec) = (double(*))PyArray_DATA(py_bvec);

  alm_get_matrix_elements(id, nat, ndata_used, amat, bvec);

  Py_RETURN_NONE;
}
