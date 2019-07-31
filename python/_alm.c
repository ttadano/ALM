#include <Python.h>
#include <assert.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "alm_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * py_alm_new(PyObject *self, PyObject *args);
static PyObject * py_alm_delete(PyObject *self, PyObject *args);
static PyObject * py_set_output_filename_prefix(PyObject *self, PyObject *args);
static PyObject * py_suggest(PyObject *self, PyObject *args);
static PyObject * py_optimize(PyObject *self, PyObject *args);
static PyObject * py_get_optimizer_control(PyObject *self, PyObject *args);
static PyObject * py_set_optimizer_control(PyObject *self, PyObject *args);
static PyObject * py_set_cell(PyObject *self, PyObject *args);
static PyObject * py_set_verbosity(PyObject *self, PyObject *args);
static PyObject * py_set_training_data(PyObject *self, PyObject *args);
static PyObject * py_get_nrows_amat(PyObject *self, PyObject *args);
static PyObject * py_get_cv_l1_alpha(PyObject *self, PyObject *args);
static PyObject * py_set_constraint_type(PyObject *self, PyObject *args);
static PyObject * py_define(PyObject *self, PyObject *args);
static PyObject * py_generate_force_constant(PyObject *self, PyObject *args);
static PyObject * py_get_atom_mapping_by_pure_translations
(PyObject *self, PyObject *args);
static PyObject * py_get_number_of_displacement_patterns
(PyObject *self, PyObject *args);
static PyObject * py_get_number_of_displaced_atoms(PyObject *self, PyObject *args);
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
  {"set_output_filename_prefix", py_set_output_filename_prefix, METH_VARARGS, ""},
  {"suggest", py_suggest, METH_VARARGS, ""},
  {"optimize", py_optimize, METH_VARARGS, ""},
  {"get_optimizer_control", py_get_optimizer_control, METH_VARARGS, ""},
  {"set_optimizer_control", py_set_optimizer_control, METH_VARARGS, ""},
  {"set_cell", py_set_cell, METH_VARARGS, ""},
  {"set_verbosity", py_set_verbosity, METH_VARARGS, ""},
  {"set_training_data", py_set_training_data, METH_VARARGS, ""},
  {"set_constraint_type", py_set_constraint_type, METH_VARARGS, ""},
  {"define", py_define, METH_VARARGS, ""},
  {"generate_force_constant", py_generate_force_constant, METH_VARARGS, ""},
  {"get_atom_mapping_by_pure_translations", py_get_atom_mapping_by_pure_translations,
   METH_VARARGS, ""},
  {"get_number_of_displacement_patterns", py_get_number_of_displacement_patterns,
   METH_VARARGS, ""},
  {"get_number_of_displaced_atoms", py_get_number_of_displaced_atoms,
   METH_VARARGS, ""},
  {"get_nrows_amat", py_get_nrows_amat, METH_VARARGS, ""},
  {"get_cv_l1_alpha", py_get_cv_l1_alpha, METH_VARARGS, ""},
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
  if (!PyArg_ParseTuple(args, "i", &id)) {
    return NULL;
  }

  alm_delete(id);

  Py_RETURN_NONE;
}

static PyObject * py_suggest(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id)) {
    return NULL;
  }

  alm_suggest(id);

  Py_RETURN_NONE;
}

static PyObject * py_set_output_filename_prefix(PyObject *self, PyObject *args)
{
  int id;
  const char *prefix;
  if (!PyArg_ParseTuple(args, "is",
                        &id,
                        &prefix)) {
    return NULL;
  }

  alm_set_output_filename_prefix(id, prefix);

  Py_RETURN_NONE;
}

static PyObject * py_optimize(PyObject *self, PyObject *args)
{
  int id, info;
  const char *solver;
  if (!PyArg_ParseTuple(args, "is",
                        &id,
                        &solver)) {
    return NULL;
  }

  info = alm_optimize(id, solver);

  return PyLong_FromLong((long) info);
}

static PyObject * py_get_optimizer_control(PyObject *self, PyObject *args)
{
  int id;
  PyObject* py_optcontrol;
  PyObject* py_value;

  PyObject *array;
  struct optimizer_control optcontrol;
  int len_list, n;

  if (!PyArg_ParseTuple(args, "iO",
                        &id,
                        &py_optcontrol)) {
    return NULL;
  }

  len_list = 15;
  array = PyList_New(len_list);
  n = 0;

  optcontrol = alm_get_optimizer_control(id);

  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.linear_model));
  n++;  /* 0 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.use_sparse_solver));
  n++;  /* 1 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.maxnum_iteration));
  n++;  /* 2 */
  PyList_SetItem(array, n, PyFloat_FromDouble(optcontrol.tolerance_iteration));
  n++;  /* 3 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.output_frequency));
  n++;  /* 4 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.standardize));
  n++;  /* 5 */
  PyList_SetItem(array, n, PyFloat_FromDouble(optcontrol.displacement_normalization_factor));
  n++;  /* 6 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.debiase_after_l1opt));
  n++;  /* 7 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.cross_validation));
  n++;  /* 8 */
  PyList_SetItem(array, n, PyFloat_FromDouble(optcontrol.l1_alpha));
  n++;  /* 9 */
  PyList_SetItem(array, n, PyFloat_FromDouble(optcontrol.l1_alpha_min));
  n++;  /* 10 */
  PyList_SetItem(array, n, PyFloat_FromDouble(optcontrol.l1_alpha_max));
  n++;  /* 11 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.num_l1_alpha));
  n++;  /* 12 */
  PyList_SetItem(array, n, PyFloat_FromDouble(optcontrol.l1_ratio));
  n++;  /* 13 */
  PyList_SetItem(array, n, PyLong_FromLong((long) optcontrol.save_solution_path));
  n++;  /* 14 */

  assert(n == len_list);

  return array;
}

static PyObject * py_set_optimizer_control(PyObject *self, PyObject *args)
{
  int id;
  PyObject* py_optcontrol;
  PyObject* py_value;
  struct optimizer_control optcontrol;
  int i;
  int updated[15];

  if (!PyArg_ParseTuple(args, "iO",
                        &id,
                        &py_optcontrol)) {
    return NULL;
  }

  for (i = 0; i < 15; i++) {
    updated[i] = 0;
  }

  if ((py_value = PyList_GetItem(py_optcontrol, 0)) != Py_None) {
    optcontrol.linear_model = (int)PyLong_AsLong(py_value);
    updated[0] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 1)) != Py_None) {
    optcontrol.use_sparse_solver = (int)PyLong_AsLong(py_value);
    updated[1] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 2)) != Py_None) {
    optcontrol.maxnum_iteration = (int)PyLong_AsLong(py_value);
    updated[2] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 3)) != Py_None) {
    optcontrol.tolerance_iteration = PyFloat_AsDouble(py_value);
    updated[3] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 4)) != Py_None) {
    optcontrol.output_frequency = (int)PyLong_AsLong(py_value);
    updated[4] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 5)) != Py_None) {
    optcontrol.standardize = (int)PyLong_AsLong(py_value);
    updated[5] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 6)) != Py_None) {
    optcontrol.displacement_normalization_factor = PyFloat_AsDouble(py_value);
    updated[6] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 7)) != Py_None) {
    optcontrol.debiase_after_l1opt = (int)PyLong_AsLong(py_value);
    updated[7] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 8)) != Py_None) {
    optcontrol.cross_validation = (int)PyLong_AsLong(py_value);
    updated[8] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 9)) != Py_None) {
    optcontrol.l1_alpha = PyFloat_AsDouble(py_value);
    updated[9] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 10)) != Py_None) {
    optcontrol.l1_alpha_min = PyFloat_AsDouble(py_value);
    updated[10] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 11)) != Py_None) {
    optcontrol.l1_alpha_max = PyFloat_AsDouble(py_value);
    updated[11] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 12)) != Py_None) {
    optcontrol.num_l1_alpha = (int)PyLong_AsLong(py_value);
    updated[12] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 13)) != Py_None) {
    optcontrol.l1_ratio = PyFloat_AsDouble(py_value);
    updated[13] = 1;
  }
  if ((py_value = PyList_GetItem(py_optcontrol, 14)) != Py_None) {
    optcontrol.save_solution_path = (int)PyLong_AsLong(py_value);
    updated[14] = 1;
  }

  py_value = NULL;

  assert(PyList_Size(py_optcontrol) == 15);

  alm_set_optimizer_control(id, optcontrol, updated);

  Py_RETURN_NONE;
}

static PyObject * py_set_cell(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_lavec;
  PyArrayObject* py_xcoord;
  PyArrayObject* py_kind_in;
  PyArrayObject* py_kind_indices;

  if (!PyArg_ParseTuple(args, "iOOOO",
                        &id,
                        &py_lavec,
                        &py_xcoord,
                        &py_kind_in,
                        &py_kind_indices)) {
    return NULL;
  }

  double (*lavec)[3] = (double(*)[3])PyArray_DATA(py_lavec);
  double (*xcoord)[3] = (double(*)[3])PyArray_DATA(py_xcoord);
  const int* kind_in = (int*)PyArray_DATA(py_kind_in);
  const size_t nat = (size_t)PyArray_DIMS(py_kind_in)[0];
  int* kind_indices = (int*)PyArray_DATA(py_kind_indices);

  alm_set_cell(id, nat, lavec, xcoord, kind_in, kind_indices);

  Py_RETURN_NONE;
}

static PyObject * py_set_verbosity(PyObject *self, PyObject *args)
{
  int id, verbosity;
  if (!PyArg_ParseTuple(args, "ii",
                        &id,
                        &verbosity)) {
    return NULL;
  }

  alm_set_verbosity(id, verbosity);

  Py_RETURN_NONE;
}

static PyObject * py_set_training_data(PyObject *self, PyObject *args)
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

  const size_t ndata_used = (size_t)PyArray_DIMS(py_f)[0];
  const size_t nat = (size_t)PyArray_DIMS(py_f)[1];

  alm_set_training_data(id, u, f, nat, ndata_used);

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

static PyObject * py_define(PyObject *self, PyObject *args)
{
  int id;
  int maxorder;

  PyArrayObject* py_nbody_include;
  PyArrayObject* py_cutoff_radii;

  size_t nkd;
  double *cutoff_radii;

  if (!PyArg_ParseTuple(args, "iiOO",
                        &id,
                        &maxorder,
                        &py_nbody_include,
                        &py_cutoff_radii)) {
    return NULL;
  }


  if ((PyObject*)py_cutoff_radii == Py_None) {
    cutoff_radii = NULL;
    nkd = 0;
  } else {
    nkd = (size_t)PyArray_DIM(py_cutoff_radii, 1);
    cutoff_radii = (double*)PyArray_DATA(py_cutoff_radii);
  }
  const int *nbody_include = (int*)PyArray_DATA(py_nbody_include);

  alm_define(id, maxorder, nkd, nbody_include, cutoff_radii);

  Py_RETURN_NONE;
}

static PyObject * py_generate_force_constant(PyObject *self, PyObject *args)
{
  int id;
  if (!PyArg_ParseTuple(args, "i", &id)) {
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

  int *map_p2s = (int*)PyArray_DATA(py_map_p2s);
  const int num_trans = alm_get_atom_mapping_by_pure_translations(id, map_p2s);

  return PyLong_FromLong((long) num_trans);
}

static PyObject * py_get_number_of_displacement_patterns
(PyObject *self, PyObject *args)
{
  int id;
  int fc_order;
  size_t num_patterns;

  if (!PyArg_ParseTuple(args, "ii",
                        &id,
                        &fc_order)) {
    return NULL;
  }

  num_patterns = alm_get_number_of_displacement_patterns(id, fc_order);

  return PyLong_FromSize_t(num_patterns);
}

static PyObject * py_get_number_of_displaced_atoms(PyObject *self, PyObject *args)
{
  int id, fc_order;
  PyArrayObject* py_numbers;

  if (!PyArg_ParseTuple(args, "iOi",
                        &id,
                        &py_numbers,
                        &fc_order)) {
    return NULL;
  }

  int (*numbers) = (int*)PyArray_DATA(py_numbers);

  alm_get_number_of_displaced_atoms(id, numbers, fc_order);

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

  int (*atom_indices) = (int*)PyArray_DATA(py_atom_indices);
  double (*disp_patterns) = (double*)PyArray_DATA(py_disp_patterns);

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

  const size_t num_fc_elems = alm_get_number_of_fc_elements(id, fc_order);

  return PyLong_FromSize_t(num_fc_elems);
}

static PyObject * py_get_number_of_irred_fc_elements(PyObject *self, PyObject *args)
{
  int id, fc_order;

  if (!PyArg_ParseTuple(args, "ii",
                        &id,
                        &fc_order)) {
    return NULL;
  }

  const size_t num_fc_elems = alm_get_number_of_irred_fc_elements(id, fc_order);

  return PyLong_FromSize_t(num_fc_elems);
}


static PyObject * py_get_nrows_amat(PyObject *self, PyObject *args)
{
    int id;

    if (!PyArg_ParseTuple(args, "i", &id)) {
        return NULL;
    }
    size_t nrows = alm_get_nrows_sensing_matrix(id);

    return PyLong_FromSize_t(nrows);
}

static PyObject * py_get_cv_l1_alpha(PyObject *self, PyObject *args)
{
    int id;

    if (!PyArg_ParseTuple(args, "i", &id)) {
        return NULL;
    }
    double cv_l1_alpha = alm_get_cv_l1_alpha(id);

    return PyFloat_FromDouble(cv_l1_alpha);
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

  double (*fc_values) = (double*)PyArray_DATA(py_fc_values);
  int (*elem_indices) = (int*)PyArray_DATA(py_elem_indices);
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

  double (*fc_values) = (double*)PyArray_DATA(py_fc_values);
  int (*elem_indices) = (int*)PyArray_DATA(py_elem_indices);
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

  double (*fc_values) = (double*)PyArray_DATA(py_fc_values);
  int (*elem_indices) = (int*)PyArray_DATA(py_elem_indices);
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

  double *fc_in = (double*)PyArray_DATA(py_fc_in);

  alm_set_fc(id, fc_in);

  Py_RETURN_NONE;
}

static PyObject * py_get_matrix_elements(PyObject *self, PyObject *args)
{
  int id;
  PyArrayObject* py_amat;
  PyArrayObject* py_bvec;

  if (!PyArg_ParseTuple(args, "iOO",
                        &id,
                        &py_amat,
                        &py_bvec)) {
    return NULL;
  }

  double (*amat) = (double*)PyArray_DATA(py_amat);
  double (*bvec) = (double*)PyArray_DATA(py_bvec);

  alm_get_matrix_elements(id, amat, bvec);

  Py_RETURN_NONE;
}

#ifdef __cplusplus
}
#endif
