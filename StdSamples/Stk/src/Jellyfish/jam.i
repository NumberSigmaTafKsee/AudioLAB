%module jam
%{
#include "cppy3.hpp"
#include "cppy3_numpy.hpp"
#include <string>
#include <vector>
#include <fstream>

using namespace cppy3;
using namespace std;
%}

%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"

using namespace std;

namespace cppy3 
{

// it is backwards
Var exec(const char* pythonScript);
Var exec(const std::wstring& pythonScript);
Var exec(const std::string& pythonScript);
Var eval(const char* pythonScript);


//Var execScriptFile(const std::wstring& path);



bool error();
//void appendToSysPath(const std::list<std::wstring>& paths);
Var createClassInstance(const std::wstring& callable);
void interrupt();
void setArgv(const std::list<std::wstring>& argv);
PyObject* getMainModule();
PyObject* getMainDict();
PyExceptionData getErrorObject(bool clearError = false);
void rethrowPythonException();
Var import(const char* moduleName, PyObject* globals = NULL, PyObject* locals = NULL);

typedef std::vector<Var> arguments;
PyObject* call(PyObject* callable, const arguments& args = arguments());
PyObject* call(const char* callable, const arguments& args = arguments());

Var lookupObject(PyObject* module, const std::wstring& name);
Var lookupCallable(PyObject* module, const std::wstring& name);

class PythonVM {

public:
  PythonVM(const std::string& programName = "");
  ~PythonVM();
};

struct PyExceptionData {
  std::wstring type;
  std::wstring reason;
  std::vector<std::wstring> trace;

  explicit PyExceptionData(const std::wstring& reason = std::wstring()) throw();
  explicit PyExceptionData(const std::wstring& type, const std::wstring& reason, const std::vector<std::wstring>& trace) throw();

  bool isEmpty() const throw();
  std::wstring toString() const;
};

class PythonException 
{
public:
  PythonException(const PyExceptionData& info_);
  PythonException(const std::wstring& reason);
  ~PythonException() throw();
  const char* what() const throw();
  const PyExceptionData info;
};


//PyObject* convert(const char* value);
//PyObject* convert(const std::wstring& value);
//PyObject* convert(const std::wstring& value);
//PyObject* convert(const int& value);
//PyObject* convert(const double& value);

/*
PyObject* convert(const std::vector<double>& value);
PyObject* convert(const std::vector<std::string>& value);
*/

void extract(PyObject* o, std::wstring& value);
void extract(PyObject* o, long& value);
void extract(PyObject* o, double& value);

//void extract(PyObject* o, std::vector<double>& value);
//void extract(PyObject* o, std::vector<std::string>& value);


class Var {
public:

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


  Var();
  explicit Var(PyObject* object);
  Var(const Var& other);
  explicit Var(const char* name, const PyObject* parent);
  explicit Var(const std::wstring& name, const PyObject* parent);
  ~Var();

  Var var(const char* name) const;
  Var var(const std::wstring& name) const;

  void reset(PyObject* o);
  void newRef(PyObject* o);
  static Var from(PyObject* o);
  PyObject* data() const;
  bool none() const;
  bool null() const;
  void decref();
  const char* typeName() const;
    
  void inject(const std::string& varName, PyObject* o);
  void inject(const std::wstring& varName, PyObject* o);


  static std::wstring toString(PyObject* val);
  std::wstring toString() const;
  std::string toUTF8String() const;
  long toLong() const;
  double toDouble() const;

};

/**
 * Adapter for python list type
 */
class List : public Var {
public:
  List(const char* name, PyObject* parent);
  List(PyObject* obj = NULL) : Var(obj);
  
  void reset(PyObject* o);
  size_t size();
  

  %extend {
      Var __getitem__(const size_t i) { return (*$self)[i]; }
  }
  
  void remove(const size_t i);
  bool contains(PyObject* element);
  void append(PyObject* element);
  void insert(const size_t index, PyObject* element);
};

class Dict : public Var {
public:

  Dict(const char* name, PyObject* parent);
  Dict(PyObject* o);

  Dict dict(const char* name);
  Dict dict(const std::wstring& name);
  List list(const char* name) const;
  Dict moduledict(const char* name);
  Dict moduledict(const std::wstring& name);
  bool contains(const char* name) const;
  void clear();
};


class Main : public Dict {
public:
  Main();
};


class GILLocker {
public:
  GILLocker();
  ~GILLocker();
};


class ScopedGILRelease {
public:
  ScopedGILRelease();
  ~ScopedGILRelease();
};


class ScopedGILLock {
public:
  ScopedGILLock();
  ~ScopedGILLock();
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
  NDArray();
  NDArray(int n);
  NDArray(int n1, int n2);
  NDArray(const Type* data, int n);
  NDArray(const Type* data, int n1, int n2);
  ~NDArray();
  void create(int n, bool fillZeros = SLOWER_AND_CLEARNER);
  void create(int n1, int n2, bool fillZeros = SLOWER_AND_CLEARNER);

  bool isset();
  void wrap(Type* data, int n);
  void wrap(Type* data, int n1, int n2);
  void copy(const Type* data, int n);
  void copy(const Type* data, size_t n1, size_t n2);

    %extend {
        Type __getitem__(int i) { return (*$self)(i); }
        Type __getitem__(int i,h) { return (*$self)(i,j); }

        void __setitem__(int i, const Type& value) { (*$self)(i) = value; }
        void __setitem__(int i, int j, const Type& value) { (*$self)(i,j) = value; }
    }
  int nd() const;
  int dim(size_t n);
  int dim1() const;
  int dim2() const;
  Type* getData();
};
}
%inline %{

inline Var execScript(const std::string& path) {
  std::ifstream t(path);
  
  std::string script((std::istreambuf_iterator<char>(t)),
                std::istreambuf_iterator<char>());
  return exec(script.c_str());
}

%}