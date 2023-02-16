/**
 * cppy3 -- embed python3 scripting layer into your c++ app in 10 minutes
 *
 * (c) 2018 Dennis Shilko
 *
 * Minimalistic library for embedding python 3 in C++ application.
 * No additional dependencies required.
 * Linux, Windows platforms supported
 *
 * - Convenient rererence-counted holder for PyObject*
 * - Simple wrapper over init/shutdown of interpreter in 1 line of code
 * - Manage GIL with scoped lock/unlock guards
 * - Translate exceptions from Python to to C++ layer
 * - C++ abstractions for list, dict and numpy.ndarray
 *
 */
#pragma once



// std
#include <exception>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include <streambuf>

#include <string>
#include <locale>

#if (defined(_WIN32) || defined(__CYGWIN__))
 #pragma warning(disable : 4996 4244 4180 4800 4661)
 // if we are building the dynamic library (instead of using it)
 #if defined(cppy2_EXPORTS)
   #define LIB_API __declspec(dllexport)
 #else
   #define LIB_API __declspec(dllimport)
 #endif
#else
 #if defined(cppy2_EXPORTS)
  #define LIB_API __attribute__((visibility("default")))
 #else
  #define LIB_API
 #endif
#endif

#ifdef CPPY3_USE_BOOST_CONVERT
    #include <boost/locale/encoding_utf.hpp>
#else
    #include <codecvt>
#endif


#ifndef _NDEBUG
 #define DLOG(MSG) {std::cerr<<MSG<<std::endl;}
 #define DWLOG(MSG) {std::wcerr<<MSG<<std::endl;}
#else
 #define DLOG(MSG) {}
 #define DWLOG(MSG) {}
#endif

// undef will surpress python warnings
#ifdef _POSIX_C_SOURCE
 #undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
 #undef _XOPEN_SOURCE
#endif
#include <Python.h>

// deal with crazy numpy 1.7.x api init procedure
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_API__CPPY3_APP_TOKEN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// python2.7/dist-packages/numpy/core/include
#include <numpy/arrayobject.h>


#include <cassert>

// fill with zeros by default new NDArray objects
#define SLOWER_AND_CLEARNER false

namespace cppy3 {

    std::wstring UTF8ToWide(const std::string& text);
    std::string WideToUTF8(const std::wstring& text);

    // forward decl
    class Var;
    class List;
    class Dict;
    class PyExceptionData;
    class PythonException;
    class Main;

    /** exec python script from text string */
    LIB_API Var exec(const char* pythonScript);
    LIB_API Var exec(const std::wstring& pythonScript);
    LIB_API Var exec(const std::string& pythonScript);
    LIB_API Var eval(const char* pythonScript);

    /** exec python script file */
    LIB_API Var execScriptFile(const std::wstring& path);

    /** @return true if python exception occured */
    LIB_API bool error();

    /** Add to sys.path */
    void appendToSysPath(const std::list<std::wstring>& paths);

    /** Make instance of a class */
    Var createClassInstance(const std::wstring& callable);

    /** Send ctrl-c */
    void interrupt();

    /** Set sys.argv */
    void setArgv(const std::list<std::wstring>& argv);


    /**
    * @returns pointer to object of root python module __main__
    * can be reached with PyImport_AddModule("__main__")
    * or PyDict_GetItemString(PyImport_GetModuleDict(), "__main__")
    * either way is right
    */
    LIB_API PyObject* getMainModule();
    LIB_API PyObject* getMainDict();

    /**
    * writes python error output to object
    * @param clearError - if true function will reset python error status before exit
    * @return python exception info and traceback in object
    */
    LIB_API PyExceptionData getErrorObject(bool clearError = false);

    /** throws c++ exception if python exception occured */
    LIB_API void rethrowPythonException();

    /** import python module into given context */
    LIB_API Var import(const char* moduleName, PyObject* globals = NULL, PyObject* locals = NULL);


    /** call python callable object and return result */
    typedef std::vector<Var> arguments;
    LIB_API PyObject* call(PyObject* callable, const arguments& args = arguments());
    LIB_API PyObject* call(const char* callable, const arguments& args = arguments());

    /** get reference to an object in python's namespace */
    LIB_API Var lookupObject(PyObject* module, const std::wstring& name);
    LIB_API Var lookupCallable(PyObject* module, const std::wstring& name);

    /**
    * Tiny wrapper over CPython interpreter instance
    * to manage init/shutdown and expose api in most simple way
    */
    class PythonVM {

    public:
    PythonVM(const std::string& programName = "");
    ~PythonVM();
    };

    struct PyExceptionData {
    std::wstring type;
    std::wstring reason;
    std::vector<std::wstring> trace;

    explicit PyExceptionData(const std::wstring& reason = std::wstring()) throw() : reason(reason) {}
    explicit PyExceptionData(const std::wstring& type, const std::wstring& reason, const std::vector<std::wstring>& trace) throw() : type(type), reason(reason), trace(trace) {}

    bool isEmpty() const throw() {
        return type.empty() && reason.empty();
    }

    std::wstring toString() const {
        std::wstring traceText;
        for (auto t: trace) {
        traceText += t + L'\n';
        }
        return isEmpty() ? std::wstring() : type + L"\n" + reason + L"\n" + traceText;
    }
    };

    class PythonException : public std::exception
    {
    public:
    PythonException(const PyExceptionData& info_) : info(info_), _what(WideToUTF8(info.toString())) {}
    PythonException(const std::wstring& reason) : info(reason), _what(WideToUTF8(info.toString())) {}
    ~PythonException() throw() {}

    const char* what() const throw() {
        return _what.c_str();
    }
    const PyExceptionData info;
    private:
    const std::string _what;
    };

    void importNumpy();


    NPY_TYPES toNumpyDType(double);
    NPY_TYPES toNumpyDType(int);

    /**
    * Simple wrapper for numpy.ndarray
    */
    template<typename Type>
    class NDArray {
    public:
    NDArray() : _ndarray(NULL) {}


    NDArray(int n) : _ndarray(NULL) {
        create(n);
    }

    NDArray(int n1, int n2) : _ndarray(NULL) {
        create(n1, n2);
    }

    NDArray(const Type* data, int n) : _ndarray(NULL) {
        copy(data, n);
    }

    NDArray(const Type* data, int n1, int n2) : _ndarray(NULL) {
        copy(data, n1, n2);
    }

    ~NDArray() {
        decref();
    }

    /**
    * Create 1d array of given size
    * @param n - size of dimension
    * @param fillZeros - initialize allocated array with zeros
    */
    void create(int n, bool fillZeros = SLOWER_AND_CLEARNER) {

        decref();

        npy_intp dim1[1];
        dim1[0] = n;
        Type impltype = 0;
        if (fillZeros) {
        _ndarray = (PyArrayObject*)PyArray_ZEROS(1, dim1, toNumpyDType(impltype), 0);
        } else {
        _ndarray = (PyArrayObject*)PyArray_SimpleNew(1, dim1, toNumpyDType(impltype));
        }
        assert(_ndarray);
    }

    /**
    * Create 2d array of given size
    * @param n1 - rows size
    * @param n2 - cols size
    * @param fillZeros - initialize allocated array with zeros
    */
    void create(int n1, int n2, bool fillZeros = SLOWER_AND_CLEARNER) {

        decref();

        npy_intp dim2[2];
        dim2[0] = n1;
        dim2[1] = n2;
        Type impltype = 0;
        if (fillZeros) {
        _ndarray = (PyArrayObject*)PyArray_ZEROS(2, dim2, toNumpyDType(impltype), 0);
        } else {
        _ndarray = (PyArrayObject*)PyArray_SimpleNew(2, dim2, toNumpyDType(impltype));
        }
        assert(_ndarray);
    }

    bool isset() {
        return (_ndarray);
    }

    /**
    * Wrap an existing 1d array, pointed to by a single Type* pointer, and wraps it in a Numpy ndarray instance
    * @param data - points to allocated array
    * @param n - number of elements of type Type in array
    */
    void wrap(Type* data, int n) {
        npy_intp dim1[1];
        dim1[0] = n;
        _ndarray = (PyArrayObject*)PyArray_SimpleNewFromData(1, dim1, toNumpyDType(*data), (void*)&data);
    }

    /**
    * Wrap an existing 2d array, pointed to by a single Type* pointer, and wraps it in a Numpy ndarray instance
    * @param data - points to allocated array
    * @param[in] n1 - number of elements of type Type in array's row
    * @param[in] n2 - number of elements of type Type in array's column
    * @param[in] data - array data
    */
    void wrap(Type* data, int n1, int n2) {
        npy_intp dim2[2];
        dim2[0] = n1;
        dim2[1] = n2;
        _ndarray = (PyArrayObject*)PyArray_SimpleNewFromData(2, dim2, toNumpyDType(*data), (void*)data);
    }

    /**
    * Create a Numpy ndarray copy of data
    * @param data - 1d array
    * @param n - size
    */
    void copy(const Type* data, int n) {
        create(n, false);

        for (int i = 0; i < n; i++) {
        *(Type*)(PyArray_GETPTR1(_ndarray, i)) = data[i];
        }
    }

    /**
    * Create a Numpy ndarray copy of data
    * @param[in] data - 2d array data
    * @param[in] n1 - rows count
    * @param[in] n2 - columns count
    */
    void copy(const Type* data, size_t n1, size_t n2) {
        create(n1, n2, false);
        for (size_t r = 0; r < n1; ++r) {
        size_t rowOffset = r * n2;
        for (size_t c = 0; c < n2; ++c) {
            *(Type*)PyArray_GETPTR2(_ndarray, r, c) = data[rowOffset + c];
        }
        }
    }

    Type& operator() (int i) {
        assert(_ndarray);
        assert(PyArray_NDIM(_ndarray) == 1 && i >= 0 && i < PyArray_DIM(_ndarray, 0));
        return *((Type*)PyArray_GETPTR1(_ndarray, i));
    }

    Type& operator() (int i, int j) {
        assert(_ndarray);
        assert(PyArray_NDIM(_ndarray) == 2);
        assert(i >= 0 && i < PyArray_DIM(_ndarray, 0));
        assert(j >= 0 && j < PyArray_DIM(_ndarray, 1));
        return *((Type*)PyArray_GETPTR2(_ndarray, i, j));
    }

    operator PyObject*() {
        assert(_ndarray);
        return (PyObject*)_ndarray;
    }
    operator PyArrayObject*() {
        assert(_ndarray);
        return _ndarray;
    }

    /**
    * @return number of dimensions
    */
    int nd() const {
        assert(_ndarray);
        return PyArray_NDIM(_ndarray);
    }

    /**
    * @return size of dimension n
    */
    int dim(size_t n) const {
        assert(_ndarray);
        assert(PyArray_NDIM(_ndarray) > n);
        return PyArray_DIM(_ndarray, n);
    }

    int dim1() const {
        return dim(1);
    }

    int dim2() const {
        return dim(2);
    }

    /**
    * @return raw pointer to data array
    */
    Type* getData() {
        assert(_ndarray);
        return PyArray_DATA(_ndarray);
    }

    private:
    PyArrayObject* _ndarray;

    void decref() {
        Py_XDECREF(_ndarray);
    }

    };


    // workaround numpy & python3 https://github.com/boostorg/python/issues/214
    static void * wrap_import_array() { import_array(); }

    void importNumpy() {
    static bool imported = false;
    if (!imported) {
        // @todo double-lock-singleton pattern against multithreaded race condition
        imported = true;
        wrap_import_array();
        rethrowPythonException();
    }
    }


    NPY_TYPES toNumpyDType(double) {
    return NPY_DOUBLE;
    }


    NPY_TYPES toNumpyDType(int) {
    return NPY_INT;
    }

    std::wstring UTF8ToWide(const std::string& text) {
    #ifdef CPPY3_USE_BOOST_CONVERT
        using boost::locale::conv::utf_to_utf;
        return utf_to_utf<wchar_t>(text.c_str(), text.c_str() + text.size());
    #else
        std::wstring_convert<std::codecvt_utf8_utf16<wchar_t, 0x10ffff, std::little_endian>> converter;
        return converter.from_bytes(text);
    #endif
    }

    std::string WideToUTF8(const std::wstring& text) {
    #ifdef CPPY3_USE_BOOST_CONVERT
        using boost::locale::conv::utf_to_utf;
        return utf_to_utf<char>(text.c_str(), text.c_str() + text.size());
    #else
        std::wstring_convert<std::codecvt_utf8<wchar_t, 0x10ffff, std::little_endian>, wchar_t> converter;
        return converter.to_bytes(text);
    #endif
    }
    /**
    * Setters / getters for access and manipulation with python vars and namespaces
    */
    LIB_API PyObject* convert(const char* value);
    LIB_API PyObject* convert(const std::wstring& value);
    LIB_API PyObject* convert(const std::wstring& value);
    LIB_API PyObject* convert(const int& value);
    LIB_API PyObject* convert(const double& value);

    template<typename T>
    PyObject* convert(const std::vector<T>& value) {
    PyObject* o = PyList_New(value.size());
    assert(o);

    for (size_t i = 0; i < value.size(); ++i) {
        PyObject* item = convert(value[i]);
        assert(item);
        int r = PyList_SetItem(o, i, item);
        assert(r == 0);
    }
    return o;
    }

    #if HAVE_MATRIX_CONTAINER
    template<typename T>
    PyObject* varCreator(const Matrix<T>& value) {
    PyObject* o = PyList_New(value.rows());
    assert(o);

    for (size_t r = 0; r < value.rows(); ++r) {
        PyObject* row = PyList_New(value.rows());
        int res = PyList_SetItem(o, r, row);
        assert(res == 0);

        for (size_t c = 0; c < value.cols(); ++c) {
        PyObject* item = varCreator(value[r][c]);
        int r = PyList_SetItem(row, c, item);
        assert(r == 0);
        }
    }
    return o;
    }
    #endif

    LIB_API void extract(PyObject* o, std::wstring& value);
    LIB_API void extract(PyObject* o, long& value);
    LIB_API void extract(PyObject* o, double& value);

    template<typename T>
    void extract(PyObject* o, std::vector<T>& value) {
    assert(o);
    value.reset(0);

    if (PyList_Check(o)) {
        value.reset(PyList_Size(o));
        for (int i = 0; i < PyList_Size(o); i++) {
        PyObject* item = PyList_GetItem(o, i);
        T itemValue;
        extract(item, itemValue);
        value[i] = itemValue;
        }
    } else {throw PythonException("variable is not a list");
    }
    }


    #if HAVE_MATRIX_CONTAINER
    template<typename T>
    void varGetter(PyObject* o, Matrix<T>& value) {
    assert(o);
    value.reset(0, 0);

    if (!PyList_Check(o)) {
        throw PyException("getMatrix(): variable is not a list");
    }

    if (PyList_Size(o) == 0) {
        value.reset(0, 0);
        return;
    }

    PyObject* row = PyList_GetItem(o, 0);
    assert(row);
    int rowSize = PyList_Size(o);
    int colSize = PyList_Size(row);
    value.reset(rowSize, colSize);

    for (int r = 0; r < rowSize; r++) {
        row = PyList_GetItem(o, r);
        if (!PyList_Check(row)) {
        throw PyException("getMatrix(): list item is not a list");
        }

        if (colSize != PyList_Size(row)) {
        throw PyException("getMatrix(): matrix rows have different size");
        }

        for (int c = 0; c < colSize; c++) {
        PyObject* item = PyList_GetItem(row, c);
        T itemValue;
        varGetter(item, itemValue);
        value[r][c] = itemValue;
        }
    }
    }
    #endif

    /**
    * In python everything is a variable object instance.
    * Base wrapper for PyObject with reference counter
    */
    class LIB_API Var {
    public:

    /**
    * Python's basic object types
    * Used mostly for internal checks
    */
    enum Type {
        UNKNOWN,
        LONG,
        BOOL,
        FLOAT,
        STRING,
        LIST,
        DICT,
        TUPLE,
        NUMPY_NDARRAY,
        MODULE
    };

    /** Construct empty */
    Var() : _o(NULL) {}

    /** Construct holder for given PyObject */
    explicit Var(PyObject* object) : _o(NULL) {
        reset(object);
    }

    /** Construct holder for given PyObject */
    Var(const Var& other) : _o(NULL) {
        reset(other.data());
    }

    /**
    * Construct holder for object parent[name]
    */
    explicit Var(const char* name, const PyObject* parent) : _o(NULL) {
        assert(name && parent);
        reset(PyDict_GetItemString(const_cast<PyObject*>(parent), name));
        assert(_o);
    }

    /**
    * Construct holder for object parent[name]
    */
    explicit Var(const std::wstring& name, const PyObject* parent) : _o(NULL) {
        assert(parent);
        reset(PyDict_GetItem(const_cast<PyObject*>(parent), Var::from(convert(name))));
        assert(_o);
    }

    ~Var() {
        decref();
    }

    /**
    * Getter for object with @b name in context of @b this
    */
    Var var(const char* name) const {
        assert(_o && "to get child object this must have parent");
        return Var(name, _o);
    }

    Var var(const std::wstring& name) const {
        assert(_o && "to get child object this must have parent");
        return Var(name, _o);
    }


    /**
    * Hold PyObject
    * reference increased
    * Used for borrowed objects
    * @param[in] o - pointer to PyObject
    */
    void reset(PyObject* o) {
        if (o != _o) {
        decref();
        _o = o;
        Py_XINCREF(_o);
        }
    }

    /**
    * Hold PyObject without increase of reference
    * Become owner of given object
    * Used mostly for owned objects
    * @param[in] o - pointer to PyObject
    */
    void newRef(PyObject* o) {
        if (o != _o) {
        decref();
        _o = o;
        }
    }

    /** static alias for newRef() */
    static Var from(PyObject* o) {
        Var v;
        v.newRef(o);
        return v;
    }

    /**
    * Get pointer to contained PyObject
    */
    PyObject* data() const {
        return _o;
    }

    /**
    * Check if contained pointer is NULL or PyObject is None
    */
    bool none() const {
        return (null() || (_o == Py_None));
    }

    /**
    * Check if contained pointer is NULL
    */
    bool null() const {
        return (_o == NULL);
    }

    /**
    * Decrement reference counter
    */
    void decref() {
        assert(!_o || (_o && (Py_REFCNT(_o) > 0)));
        Py_XDECREF(_o);
    }

    /**
    * Get text representation of PyObject type
    */
    const char* typeName() const {
        assert(_o);
        return _o->ob_type->tp_name;
    }

    /**
    * PyObject access operators
    */
    PyObject* operator->() const {
        return data();
    }

    PyObject& operator*() const {
        return *data();
    }

    operator PyObject*() const {
        return data();
    }


    /**
    * Make python object and inject into this
    */
    template<typename T>
    void injectVar(const std::wstring& varName, const T& value) {
        PyObject* o = convert(value);
        inject(varName, o);
    }

    /**
    * Make python object and inject into this
    */
    template<typename T>
    void injectVar(const std::string& varName, const T& value) {
        PyObject* o = convert(value);
        inject(varName, o);
    }

    /**
    * inject object
    */
    void inject(const std::string& varName, PyObject* o) {
        int r = PyDict_SetItemString(*this, varName.c_str(), o);
        r = r; // surpress compiler warning in release build
        assert(r == 0);
    }

    /**
    * inject object
    */
    void inject(const std::wstring& varName, PyObject* o) {
        inject(WideToUTF8(varName), o);
    }



    /**
    * Get PyObject with @b varName in @b this context convert data to @b value
    */
    template<typename T>
    void getVar(const std::wstring& varName, T& value) const {
        PyObject* o = PyDict_GetItemString(*this, WideToUTF8(varName).data());
        assert(o);
        extract(o, value);
    }

    template<typename T>
    void getList(const std::wstring& varName, std::vector<T>& value) const {
        PyObject* o = PyDict_GetItemString(*this, WideToUTF8(varName).data());
        assert(o);
        extract(o, value);
    }

    #if HAVE_MATRIX_CONTAINER
    template<typename T>
    void getMatrix(const std::wstring& varName, Matrix<T>& value) const {
        PyObject* o = PyDict_GetItemString(*this, WideToUTF8(varName).data());
        assert(o);
        varGetter(o, value);
    }
    #endif

    /**
    * Get text representation. Equal to python str() or repr()
    */
    static std::wstring toString(PyObject* val);
    std::wstring toString() const;
    std::string toUTF8String() const {return WideToUTF8(toString());}
    
    /**
    * Get cast to scalar POD types
    */
    long toLong() const;
    double toDouble() const;

    /**
    * Get type
    */
    Type type() const;

    protected:
    PyObject* _o;

    };

    /**
    * Adapter for python list type
    */
    class LIB_API List : public Var {
    public:
    List(const char* name, PyObject* parent) : Var(name, parent) {
        _validate();
    }

    List(PyObject* obj = NULL) : Var(obj) {
        _validate();
    }

    void reset(PyObject* o) {
        Var::reset(o);
        _validate();
    }

    size_t size() {
        const Py_ssize_t size = PyList_Size(_o);
        if (size == -1) {
        throw PythonException(L"invalid python list object");
        }
        return size;
    }

    Var operator[](const size_t i) {
        if (i >= size()) {
        throw PythonException(L"List index of of bounds");
        }
        return Var(PyList_GetItem(_o, i));
    }

    void remove(const size_t i) {
        const int result = PySequence_DelItem(_o, i);
        if (result == -1) {
        throw PythonException(L"PySequence_DelItem error");
        }
    }

    bool contains(PyObject* element) {
        const int result = PySequence_Contains(_o, element);
        if (result == -1) {
        throw PythonException(L"PySequence_Contains failed on list object");
        }
        return result;
    }

    void append(PyObject* element) {
        const int result = PyList_Append(_o, element);
        if (result == -1) {
        throw PythonException(L"PyList_Append failed on list object");
        }
    }

    void insert(const size_t index, PyObject* element) {
        const int result = PyList_Insert(_o, index, element);
        if (result == -1) {
        throw PythonException(L"PyList_Insert failed on list object");
        }
    }

    private:
    void _validate() const {
        if (_o) {
        assert(type() == LIST);
        }
    }

    };


    /**
    * Adapter for python dict type
    */
    class LIB_API Dict : public Var {
    public:
    Dict(const char* name, PyObject* parent) : Var(name, parent) {}

    Dict(PyObject* o) : Var(o) {
        assert(type() == DICT);
    }

    /**
    * Accessors to nested variables
    */
    Dict dict(const char* name) {
        assert(Var(name, _o).type() == DICT);
        return Dict(name, _o);
    }

    Dict dict(const std::wstring& name) { return dict(WideToUTF8(name).c_str()); }


    List list(const char* name) const {
        assert(Var(name, _o).type() == LIST);
        return List(name, _o);
    }

    Dict moduledict(const char* name) {
        assert(Var(name, _o).type() == MODULE);
        return Dict(PyModule_GetDict(Var(name, _o)));
    }

    Dict moduledict(const std::wstring& name) { return moduledict(WideToUTF8(name).c_str()); }


    /** @} */

    bool contains(const char* name) const {
        assert(type() == DICT);
        Var key;
        key.newRef(convert(name));
        return PyDict_Contains(_o, key.data());
    }

    void clear() {
        PyDict_Clear(_o);
    }

    };


    /**
    * Adapter for python root '__main__' namespace dict
    */
    class LIB_API Main : public Dict {
    public:
    Main() : Dict(getMainDict()) {}
    };


    /**
    * GIL state scoped-lock
    * can be used recursively (like recursive mutex)
    */
    class LIB_API GILLocker {
    public:
    GILLocker();
    ~GILLocker();
    private:
    void lock();
    void release();
    bool _locked;
    PyGILState_STATE _pyGILState;
    };


    /**
    * @brief The Scoped GIL unlocker
    */
    class ScopedGILRelease {
    public:
    ScopedGILRelease() {
        _threadState = PyEval_SaveThread();
    }

    ~ScopedGILRelease() {
        PyEval_RestoreThread(_threadState);
    }
    private:
    PyThreadState* _threadState;
    };


    /**
    * @brief The Scoped GIL locker
    */
    class ScopedGILLock {
    public:
    ScopedGILLock() {
        _state = PyGILState_Ensure();
    }

    ~ScopedGILLock() {
        PyGILState_Release(_state);
    }
    private:
    PyGILState_STATE _state;
    };




    PythonVM::PythonVM(const std::string& programName) {

    setenv("PYTHONDONTWRITEBYTECODE", "1", 0);

    #ifdef _WIN32
    // force utf-8 on windows
    setenv("PYTHONIOENCODING", "UTF-8");
    #endif

    Py_SetProgramName(const_cast<wchar_t*>(UTF8ToWide(programName).c_str()));

    // create CPython instance without registering signal handlers
    Py_InitializeEx(0);

    // initialize GIL
    PyEval_InitThreads();
    }

    PythonVM::~PythonVM() {
    if (!PyImport_AddModule("dummy_threading")) {
        PyErr_Clear();
    }
    Py_Finalize();
    }


    void setArgv(const std::list<std::wstring>& argv) {
    GILLocker lock;

    std::vector<const wchar_t*> cargv;
    for (std::list<std::wstring>::const_iterator it = argv.begin(); it != argv.end(); ++it) {
        cargv.push_back(it->data());
    }
    PySys_SetArgvEx(argv.size(), (wchar_t**)(&cargv[0]), 0);
    }


    Var createClassInstance(const std::wstring& callable) {
    GILLocker lock;
    Var instance;

    instance.newRef(call(lookupCallable(getMainModule(), callable)));
    rethrowPythonException();
    if (instance.none()) {
        std::wstringstream ss;
        ss << L"error instantiating '"<< callable <<"': "<< getErrorObject().toString();
        throw PythonException(ss.str());
    }
    return instance;
    }


    void appendToSysPath(const std::vector<std::wstring>& paths) {
    GILLocker lock;

    Var sys = import("sys");
    List sysPath(lookupObject(sys, L"path"));
    for (auto path: paths) {
        Var pyPath(convert(path));
        if (!sysPath.contains(pyPath)) {
        // append into the 'sys.path'
        sysPath.append(pyPath);
        }
    }
    }

    void interrupt() {
    PyErr_SetInterrupt();
    }


    Var exec(const char* pythonScript) {
    GILLocker lock;
    PyObject* mainDict = getMainDict();
    Var result;

    result.newRef(PyRun_String(pythonScript, Py_file_input, mainDict, mainDict));
    if (result.data() == NULL) {
        rethrowPythonException();
    }
    return result;
    }

    LIB_API Var eval(const char* pythonScript) {
    GILLocker lock;
    PyObject* mainDict = getMainDict();
    Var result;

    result.newRef(PyRun_String(pythonScript, Py_eval_input, mainDict, mainDict));
    if (result.data() == NULL) {
        const PyExceptionData excData = getErrorObject(false);
        if (excData.type == L"<class 'SyntaxError'>") {
        // eval() throws SyntaxError when called with expressions
        // use exec() for expressions
        getErrorObject(true);
        return exec(pythonScript);
        } else {
        rethrowPythonException();
        }
    }
    return result;
    }

    Var exec(const std::string& pythonScript) {
    return exec(pythonScript.data());
    }

    Var exec(const std::wstring& pythonScript) {
    ScopedGILLock lock;

    // encode unicode std::wstring to utf8
    std::wstring script = L"# -*- coding: utf-8 -*-\n";
    script += pythonScript;
    return exec(WideToUTF8(script).c_str());
    }

    Var execScriptFile(const std::wstring& path) {
    std::ifstream t("file.txt");

    if (!t.is_open()) {
        throw PythonException(L"canot open file " + path);
    }

    std::string script((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
    return exec(script.c_str());
    }

    bool error() {
    GILLocker lock;
    return Py_IsInitialized() && PyErr_Occurred();
    }

    void rethrowPythonException() {
    if (error()) {
        const PyExceptionData excData = getErrorObject(true);
        throw PythonException(excData);
    }
    }

    std::wstring pyUnicodeToWstring(PyObject* object) {
    std::wstring result;
    if (PyUnicode_Check(object)) {
        PyObject* bytes = PyUnicode_AsEncodedString(object, "UTF-8", "strict"); // Owned reference
        if (bytes != NULL) {
        char* utf8String = PyBytes_AS_STRING(bytes); // Borrowed pointer
        result = UTF8ToWide(std::string(utf8String));
        Py_DECREF(bytes);
        }
    }
    return result;
    }

    PyExceptionData getErrorObject(const bool clearError) {
    GILLocker lock;
    std::wstring exceptionType;
    std::wstring exceptionMessage;
    std::vector<std::wstring> exceptionTrace;
    if (PyErr_Occurred()) {
        // get error context
        PyObject* excType = NULL;
        PyObject* excValue = NULL;
        PyObject* excTraceback = NULL;
        PyErr_Fetch(&excType, &excValue, &excTraceback);
        PyErr_NormalizeException(&excType, &excValue, &excTraceback);

        // get traceback module
        PyObject* name = PyUnicode_FromString("traceback");
        PyObject* tracebackModule = PyImport_Import(name);
        Py_DECREF(name);

        // write text type of exception
        exceptionType = pyUnicodeToWstring(PyObject_Str(excType));

        // write text message of exception
        exceptionMessage = pyUnicodeToWstring(PyObject_Str(excValue));

        if (excTraceback != NULL && tracebackModule != NULL) {
        // get traceback.format_tb() function ptr
        PyObject* tbDict = PyModule_GetDict(tracebackModule);
        PyObject* format_tbFunc = PyDict_GetItemString(tbDict, "format_tb");
        if (format_tbFunc && PyCallable_Check(format_tbFunc)) {
            // build argument
            PyObject* excTbTupleArg = PyTuple_New(1);
            PyTuple_SetItem(excTbTupleArg, 0, excTraceback);
            Py_INCREF(excTraceback); // because PyTuple_SetItem() steals reference
            // call traceback.format_tb(excTraceback)
            PyObject* list = PyObject_CallObject(format_tbFunc, excTbTupleArg);
            if (list != NULL) {
            // parse list and extract traceback text lines
            const int len = PyList_Size(list);
            for (int i = 0; i < len; i++) {
                PyObject* tt = PyList_GetItem(list, i);
                PyObject* t = Py_BuildValue("(O)", tt);
                char* buffer = NULL;
                if (PyArg_ParseTuple(t, "s", &buffer)) {
                exceptionTrace.push_back(UTF8ToWide(buffer));
                }
                Py_XDECREF(t);
            }
            Py_DECREF(list);
            }
            Py_XDECREF(excTbTupleArg);
        }
        }
        Py_XDECREF(tracebackModule);

        if (clearError) {
        Py_XDECREF(excType);
        Py_XDECREF(excValue);
        Py_XDECREF(excTraceback);
        } else {
        PyErr_Restore(excType, excValue, excTraceback);
        }
    }
    return PyExceptionData(exceptionType, exceptionMessage, exceptionTrace);
    }


    PyObject* convert(const int& value) {
    PyObject* o = PyLong_FromLong(value);
    assert(o);
    return o;
    }


    PyObject* convert(const double& value) {
    PyObject* o = PyFloat_FromDouble(value);
    assert(o);
    return o;
    }


    PyObject* convert(const char* value) {
    PyObject* o = PyUnicode_FromString(value);
    return o;
    }

    PyObject* convert(const std::wstring& value) {
    PyObject* o = PyUnicode_FromWideChar(value.data(), value.size());
    return o;
    }

    void extract(PyObject* o, std::wstring& value) {
    Var str(o);
    if (!PyUnicode_Check(o)) {
        // try cast to string
        str.newRef(PyObject_Str(o));
        if (!str.data()) {
        throw PythonException(L"variable has no string representation");
        }
    }
    value = PyUnicode_AsUnicode(str);
    }


    void extract(PyObject* o, double& value) {
    if (PyFloat_Check(o)) {
        value = PyFloat_AsDouble(o);
    } else {
        throw PythonException(L"variable is not a real type");
    }
    }


    void extract(PyObject* o, long& value) {
    if (PyLong_Check(o)) {
        value = PyLong_AsLong(o);
    } else {
        throw PythonException(L"variable is not a long type");
    }
    }


    std::wstring Var::toString() const {
    return toString(_o);
    }


    std::wstring Var::toString(PyObject* val) {
    assert(val);
    std::wstringstream result;
    // try str() operator
    PyObject* str = PyObject_Str(val);
    if (!str) {
        // try repr() operator
        str = PyObject_Repr(val);
    }
    if (str) {
        result << pyUnicodeToWstring(str);
    } else {
        result << "< type='" << val->ob_type->tp_name << L"' has no string representation >";
    }
    return result.str();
    }


    long Var::toLong() const {
    long value = 0;
    extract(_o, value);
    return value;
    }


    double Var::toDouble() const {
    long value = 0;
    extract(_o, value);
    return value;
    }


    Var::Type Var::type() const {
    if (PyLong_Check(_o)) {
        return LONG;
    } else if (PyFloat_Check(_o)) {
        return FLOAT;
    } else if (PyUnicode_Check(_o)) {
        return STRING;
    } else if (PyTuple_Check(_o)) {
        return TUPLE;
    } else if (PyDict_Check(_o)) {
        return DICT;
    } else if (PyList_Check(_o)) {
        return LIST;
    } else if (PyBool_Check(_o)) {
        return BOOL;
    }
    #ifdef NPY_NDARRAYOBJECT_H
    else if (PyArray_Check(_o)) {
        return NUMPY_NDARRAY;
    }
    #endif
    else if (PyModule_Check(_o)) {
        return MODULE;
    } else {
        return UNKNOWN;
    }
    }


    Var import(const char* moduleName, PyObject* globals, PyObject* locals) {
    Var module;
    module.newRef(PyImport_ImportModuleEx(const_cast<char*>(moduleName), globals, locals, NULL));
    if (module.null()) {
        const PyExceptionData excData = getErrorObject();
        throw PythonException(excData);
    }
    return module;
    }


    Var lookupObject(PyObject* module, const std::wstring& name) {
    std::wstring temp;
    std::vector<std::wstring> items;
    std::wstringstream wss(name);
    while(std::getline(wss, temp, L'.'))
        items.push_back(temp);

    Var p(module);
    // Var prev;  // (1) cause refcount bug
    std::string itemName;
    for (auto it = items.begin(); it != items.end() && !p.null(); ++it) {

        // prev = p; // (2) cause refcount bug
        itemName = WideToUTF8(*it);
        if (PyDict_Check(p)) {
        p = Var(PyDict_GetItemString(p, itemName.data()));
        } else {
        // PyObject_GetAttrString returns new reference
        p.newRef(PyObject_GetAttrString(p, itemName.data()));
        }

        if (p.null()) {
        std::wstringstream wss;
        wss << L"lookup " << name << L" failed: no item " << UTF8ToWide(itemName);
        throw PythonException(wss.str());
        }
    }
    return p;
    }


    Var lookupCallable(PyObject* module, const std::wstring& name) {
    Var p = lookupObject(module, name);

    if (!PyCallable_Check(p)) {
        std::wstringstream wss;
        wss << L"PyObject " << name << L" is not callable";
        throw PythonException(wss.str());
    }

    return p;
    }


    PyObject* call(PyObject* callable, const arguments& args) {
    assert(callable);
    if (!PyCallable_Check(callable)) {
        std::wstringstream wss;
        wss << L"PyObject " << callable << L" is not callable";
        throw PythonException(wss.str());
    }

    PyObject* result = NULL;
    Var argsTuple;
    const int argsCount = args.size();
    if (argsCount > 0) {
        argsTuple.newRef(PyTuple_New(argsCount));
        for (int i = 0; i < argsCount; i++) {
        // steals reference
        PyTuple_SetItem(argsTuple, i, args[i]);
        }
    }

    PyErr_Clear();
    result = PyObject_CallObject(callable, argsTuple);
    rethrowPythonException();

    return result;
    }


    PyObject* call(const char* callable, const arguments& args) {
    return call(lookupCallable(getMainModule(), UTF8ToWide(callable)), args);
    }


    GILLocker::GILLocker() :
    _locked(false) {
    // autolock GIL in scoped_lock style
    lock();
    }


    GILLocker::~GILLocker() {
    release();
    }


    void GILLocker::release() {
    if (_locked) {
        assert(Py_IsInitialized());
        PyGILState_Release(_pyGILState);
        _locked = false;
    }
    }


    void GILLocker::lock() {
    if (!_locked) {
        assert(Py_IsInitialized());
        _pyGILState = PyGILState_Ensure();
        _locked = true;
    }
    }


    PyObject* getMainModule() {
    PyObject* mainModule = PyImport_AddModule("__main__");
    assert(mainModule);
    return mainModule;
    }


    PyObject* getMainDict() {
    PyObject* mainDict = PyModule_GetDict(getMainModule());
    assert(mainDict);
    return mainDict;
    }
}