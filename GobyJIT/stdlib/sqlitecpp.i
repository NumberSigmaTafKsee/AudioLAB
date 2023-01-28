%module sqlitecpp
%{
#include "SQLiteCpp/SQLiteCpp.h"
using namespace SQLite;
%}



//rm%import "SQLiteCpp/SQLiteCpp.h"
#undef  SQLITECPP_PURE_FUNC
#define SQLITECPP_PURE_FUNC
%include <SQLiteCpp/Assertion.h>
%include <SQLiteCpp/Exception.h>
%include <SQLiteCpp/Column.h>
%include <SQLiteCpp/Database.h>
%include <SQLiteCpp/Statement.h>
%include <SQLiteCpp/Transaction.h>


