#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "alm_wrapper.h"

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * py_get_version(PyObject *self, PyObject *args);
static PyObject * py_get_fc(PyObject *self, PyObject *args);

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
  {"version", py_get_version, METH_VARARGS, "Alm version"},
  {"get_fc", py_get_fc, METH_VARARGS, ""},
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

static PyObject * py_get_fc(PyObject *self, PyObject *args)
{
  int fc_order;
  PyArrayObject* py_fc_value;
  PyArrayObject* py_elem_indices;

  if (!PyArg_ParseTuple(args, "OOi",
			&py_fc_value,
			&py_elem_indices,
			&fc_order)) {
    return NULL;
  }

  double (*fc_value) = (double(*))PyArray_DATA(fc_value);
  int (*elem_indices) = (int(*))PyArray_DATA(elem_indices);

  alm_get_fc(fc_value, elem_indices, fc_order);

  Py_RETURN_NONE;
}

static PyObject * py_get_fc_length(PyObject *self, PyObject *args)
{
  int fc_order;

  if (!PyArg_ParseTuple(args, "i",
			&fc_order)) {
    return NULL;
  }

  alm_get_fc_length(fc_order);

  return PyLong_FromLong((long) fc_order);
}

static PyObject * py_get_fc_length(PyObject *self, PyObject *args)
{
  PyArrayObject* py_lavec;
  PyArrayObject* py_xcoord;
  PyArrayObject* py_kd;
  if (!PyArg_ParseTuple(args, "OO",
			&py_lavec,
			&py_xcoord)) {
    return NULL;
  }
}
