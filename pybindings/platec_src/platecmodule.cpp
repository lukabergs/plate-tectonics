#include "platecapi.hpp"
#ifdef __MINGW32__ // this is to avoid a problem with the hypot function which is messed up by Python...
#undef __STRICT_ANSI__
#endif
#include <cmath>
#include <Python.h>
#include <vector>

#define __STDC_CONSTANT_MACROS
#if _WIN32 || _WIN64
#include <Windows.h>
typedef UINT32 uint32_t;
typedef INT32 int32_t;
#else
#include <stdint.h>
#endif

static PyObject * platec_create(PyObject *self, PyObject *args, PyObject *kwargs)
{
    unsigned int seed;
    unsigned int width;
    unsigned int height;
    float sea_level;
    unsigned int erosion_period;
    float folding_ratio;
    unsigned int aggr_overlap_abs;
    float aggr_overlap_rel;
    unsigned int cycle_count;
    unsigned int num_plates;

    static char *kwlist[] = {
        (char*)"seed",
        (char*)"width",
        (char*)"height",
        (char*)"sea_level",
        (char*)"erosion_period",
        (char*)"folding_ratio",
        (char*)"aggr_overlap_abs",
        (char*)"aggr_overlap_rel",
        (char*)"cycle_count",
        (char*)"num_plates",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "IIIfIfIfII", kwlist,
                                     &seed, &width, &height, &sea_level, &erosion_period,
                                     &folding_ratio, &aggr_overlap_abs, &aggr_overlap_rel,
                                     &cycle_count, &num_plates))
        return nullptr;
    srand(seed);

    void *litho = platec_api_create(seed, width, height, sea_level, erosion_period,
                                    folding_ratio, aggr_overlap_abs, aggr_overlap_rel,
                                    cycle_count, num_plates);

    Py_ssize_t pointer = (Py_ssize_t)litho;
    return Py_BuildValue("n", pointer);
}

static PyObject * platec_step(PyObject *self, PyObject *args)
{
    void *litho;
    if (!PyArg_ParseTuple(args, "n", &litho))
        return nullptr;
    platec_api_step(litho);
    return Py_BuildValue("i", 0);
}

static PyObject * platec_destroy(PyObject *self, PyObject *args)
{
    void *litho;
    if (!PyArg_ParseTuple(args, "n", &litho))
        return nullptr;
    platec_api_destroy(litho);
    return Py_BuildValue("i", 0);
}

static PyObject * platec_load_heightmap_u16(PyObject *self, PyObject *args)
{
    void *litho;
    PyObject* sequence_obj;
    unsigned int sea_level_m;
    if (!PyArg_ParseTuple(args, "nOI", &litho, &sequence_obj, &sea_level_m))
        return nullptr;
    if (sea_level_m > 65535U) {
        PyErr_SetString(PyExc_ValueError, "sea_level_m must be <= 65535");
        return nullptr;
    }

    PyObject* sequence = PySequence_Fast(sequence_obj, "heightmap must be a sequence");
    if (sequence == nullptr) {
        return nullptr;
    }

    const Py_ssize_t expected =
        static_cast<Py_ssize_t>(lithosphere_getMapWidth(litho)) *
        static_cast<Py_ssize_t>(lithosphere_getMapHeight(litho));
    if (PySequence_Fast_GET_SIZE(sequence) != expected) {
        Py_DECREF(sequence);
        PyErr_SetString(PyExc_ValueError, "heightmap length does not match simulation dimensions");
        return nullptr;
    }

    std::vector<uint16_t> heightmap(static_cast<size_t>(expected));
    PyObject** items = PySequence_Fast_ITEMS(sequence);
    for (Py_ssize_t i = 0; i < expected; ++i) {
        const unsigned long value = PyLong_AsUnsignedLong(items[i]);
        if (PyErr_Occurred() != nullptr) {
            Py_DECREF(sequence);
            return nullptr;
        }
        if (value > 65535UL) {
            Py_DECREF(sequence);
            PyErr_SetString(PyExc_ValueError, "heightmap samples must be in [0, 65535]");
            return nullptr;
        }
        heightmap[static_cast<size_t>(i)] = static_cast<uint16_t>(value);
    }
    Py_DECREF(sequence);

    platec_api_load_heightmap_u16(litho, heightmap.data(), static_cast<uint16_t>(sea_level_m));
    return Py_BuildValue("i", 0);
}

PyObject *makelist(float array[], size_t size) {
    PyObject *l = PyList_New(size);
    for (size_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l, i, Py_BuildValue("f",array[i]));
    }
    return l;
}

PyObject *makelist_int(uint32_t array[], uint32_t size) {
    PyObject *l = PyList_New(size);
    for (uint32_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l, i, Py_BuildValue("i",array[i]));
    }
    return l;
}

static PyObject * platec_get_heightmap(PyObject *self, PyObject *args)
{
    void *litho;
    if (!PyArg_ParseTuple(args, "n", &litho))
        return nullptr;
    float *hm = platec_api_get_heightmap(litho);

    size_t width = lithosphere_getMapWidth(litho);
    size_t height = lithosphere_getMapHeight(litho);

    PyObject* res =  makelist(hm,width*height);
    Py_INCREF(res);
    return res;
}

static PyObject * platec_get_platesmap(PyObject *self, PyObject *args)
{
    void *litho;
    if (!PyArg_ParseTuple(args, "n", &litho))
        return nullptr;
    uint32_t *hm = platec_api_get_platesmap(litho);

    size_t width = lithosphere_getMapWidth(litho);
    size_t height = lithosphere_getMapHeight(litho);

    PyObject* res =  makelist_int(hm,width*height);
    Py_INCREF(res);
    return res;
}

static PyObject * platec_is_finished(PyObject *self, PyObject *args)
{
    void *litho;
    if (!PyArg_ParseTuple(args, "n", &litho))
        return nullptr;
    PyObject* res = Py_BuildValue("b",platec_api_is_finished(litho));
    return res;
}

static PyObject * platec_get_sea_level_m(PyObject *self, PyObject *args)
{
    void *litho;
    if (!PyArg_ParseTuple(args, "n", &litho))
        return nullptr;
    return Py_BuildValue("I", static_cast<unsigned int>(platec_api_get_sea_level_m(litho)));
}

static PyMethodDef PlatecMethods[] = {
    {   "create",  (PyCFunction)platec_create, METH_VARARGS | METH_KEYWORDS,
        "Create initial plates configuration."
    },
    {   "destroy",  platec_destroy, METH_VARARGS,
        "Release the data for the simulation."
    },
    {   "load_heightmap_u16",  platec_load_heightmap_u16, METH_VARARGS,
        "Load a uint16 metric heightmap into an existing simulation."
    },
    {   "get_heightmap",  platec_get_heightmap, METH_VARARGS,
        "Get current heightmap."
    },
    {   "get_sea_level_m",  platec_get_sea_level_m, METH_VARARGS,
        "Get the active metric sea level threshold."
    },
    {   "get_platesmap",  platec_get_platesmap, METH_VARARGS,
        "Get current plates map."
    },
    {   "step", platec_step, METH_VARARGS,
        "Perform next step of the simulation."
    },
    {   "is_finished",  platec_is_finished, METH_VARARGS,
        "Is the simulation finished?"
    },
    {nullptr, nullptr, 0, nullptr} /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "platec",     /* m_name */
    "Plate tectonics simulation",  /* m_doc */
    -1,                  /* m_size */
    PlatecMethods,       /* m_methods */
    nullptr,                /* m_reload */
    nullptr,                /* m_traverse */
    nullptr,                /* m_clear */
    nullptr,                /* m_free */
};

PyMODINIT_FUNC PyInit_platec(void)
{
    return PyModule_Create(&moduledef);
}
