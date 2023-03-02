%module jsoncpp
%{
#include <json/json.h>

using namespace std;
using namespace Json;
#include "StdJson.h"

%}

%include "stdint.i"
%include "std_vector.i"
%include "std_string.i"
%include "lua_fnptr.i"

%import "json/json.h"

enum ValueType {
  nullValue = 0, ///< 'null' value
  intValue,      ///< signed integer value
  uintValue,     ///< unsigned integer value
  realValue,     ///< double value
  stringValue,   ///< UTF-8 string value
  booleanValue,  ///< bool value
  arrayValue,    ///< array value (ordered list)
  objectValue    ///< object value (collection of name/value pairs).
};

%include "StdJson.h"

