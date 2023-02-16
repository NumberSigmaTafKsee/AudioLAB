
%module json
%{
#define JSON_HAS_INT64
#include "json.h"
using namespace Json;
%}

%include "std_string.i"

namespace Json
{

using String = std::basic_string<char, std::char_traits<char>, Allocator<char>>;
using Int = int;
using UInt = unsigned int;
using Int64 = int64_t;
using UInt64 = uint64_t;
using LargestInt = Int64;
using LargestUInt = UInt64;
using ArrayIndex = unsigned int;

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

enum CommentPlacement {
  commentBefore = 0,      ///< a comment placed on the line before a value
  commentAfterOnSameLine, ///< a comment just after a value on the same line
  commentAfter, ///< a comment on the line after a value (only make sense for
  /// root value)
  numberOfCommentPlacement
};

/** \brief Type of precision for formatting of real values.
 */
enum PrecisionType {
  significantDigits = 0, ///< we set max number of significant digits in string
  decimalPlaces          ///< we set max number of digits after "." in string
};

class Features {
public:

  static Features all();
  static Features strictMode();

  Features();
  bool allowComments_{true};
  bool strictRoot_{false};
  bool allowDroppedNullPlaceholders_{false};
  bool allowNumericKeys_{false};
};

class Value {
public:
  using Members = std::vector<String>;
  using iterator = ValueIterator;
  using const_iterator = ValueConstIterator;
  using UInt = Json::UInt;
  using Int = Json::Int;
  using UInt64 = Json::UInt64;
  using Int64 = Json::Int64;
  using LargestInt = Json::LargestInt;
  using LargestUInt = Json::LargestUInt;
  using ArrayIndex = Json::ArrayIndex;
  using value_type = std::string;

  static Value const& nullSingleton();
  static constexpr LargestInt minLargestInt = LargestInt(~(LargestUInt(-1) / 2));  
  static constexpr LargestInt maxLargestInt = LargestInt(LargestUInt(-1) / 2);  
  static constexpr LargestUInt maxLargestUInt = LargestUInt(-1);
  static constexpr Int minInt = Int(~(UInt(-1) / 2));
  static constexpr Int maxInt = Int(UInt(-1) / 2);
  static constexpr UInt maxUInt = UInt(-1);

  static constexpr Int64 minInt64 = Int64(~(UInt64(-1) / 2));
  static constexpr Int64 maxInt64 = Int64(UInt64(-1) / 2);
  static constexpr UInt64 maxUInt64 = UInt64(-1);
  static constexpr UInt defaultRealPrecision = 17;
  static constexpr double maxUInt64AsDouble = 18446744073709551615.0;

public:

  Value(ValueType type = nullValue);
  Value(Int value);
  Value(UInt value);
  Value(Int64 value);
  Value(UInt64 value);
  Value(double value);
  Value(const char* value); ///< Copy til first 0. (NULL causes to seg-fault.)
  Value(const char* begin, const char* end); ///< Copy all, incl zeroes.
  //Value(const StaticString& value);
  Value(const String& value);
  Value(bool value);  
  Value(const Value& other);
  
  ~Value();

  Value& operator=(const Value& other);
  
  void swap(Value& other);
  //void swapPayload(Value& other);
  void copy(const Value& other);
  //void copyPayload(const Value& other);

  ValueType type() const;

  bool operator<(const Value& other) const;
  bool operator<=(const Value& other) const;
  bool operator>=(const Value& other) const;
  bool operator>(const Value& other) const;
  bool operator==(const Value& other) const;
  bool operator!=(const Value& other) const;

  int compare(const Value& other) const;

  const char* asCString() const; ///< Embedded zeroes could cause you trouble!

  String asString() const; ///< Embedded zeroes are possible.
  bool getString(char const** begin, char const** end) const;
  Int asInt() const;
  UInt asUInt() const;
  Int64 asInt64() const;
  UInt64 asUInt64() const;
  LargestInt asLargestInt() const;
  LargestUInt asLargestUInt() const;
  float asFloat() const;
  double asDouble() const;
  bool asBool() const;


  bool isNull() const;
  bool isBool() const;
  bool isInt() const;
  bool isInt64() const;
  bool isUInt() const;
  bool isUInt64() const;
  bool isIntegral() const;
  bool isDouble() const;
  bool isNumeric() const;
  bool isString() const;
  bool isArray() const;
  bool isObject() const;

  bool isConvertibleTo(ValueType other) const;
  ArrayIndex size() const;
  bool empty() const;
  explicit operator bool() const;
  void clear();
  void resize(ArrayIndex newSize);
 
  %extend {
        Value __getitem__(ArrayIndex index) { return (*$self)[index]; }
        void  __setitem__(ArrayIndex indx, Value v) { (*$self)[indx] = v; }        
        Value __getitem__(int index)  { return (*$self)[index]; }
        void  __setitem__(int indx, Value v) { (*$self)[indx] = v; }        
        Value __getitem__(const char * index) { return (*$self)[index]; }
        void  __setitem__(const char * indx, Value v) { (*$self)[indx] = v; }        
        //Value __getitem__(const String& index) { return (*$self)[index]; }
        //void  __setitem__(const String& indx, Value v) { (*$self)[indx] = v; }        
        //Value __getitem__(const StaticString& key) { return (*$self)[index]; }
        //void  __setitem__(const StaticString& indx, Value v) { (*$self)[indx] = v; }        

    }
    
  
  Value get(ArrayIndex index, const Value& defaultValue) const;
  bool isValidIndex(ArrayIndex index) const;
  Value& append(const Value& value);
  Value& append(Value&& value);
  bool insert(ArrayIndex index, const Value& newValue);
  bool insert(ArrayIndex index, Value&& newValue);

  Value get(const char* key, const Value& defaultValue) const;
  Value get(const char* begin, const char* end,
            const Value& defaultValue) const;
  //Value get(const String& key, const Value& defaultValue) const;
  Value const* find(char const* begin, char const* end) const;
  
  
  //Value* demand(char const* begin, char const* end);
  void removeMember(const char* key);
  //void removeMember(const String& key);
  bool removeMember(const char* key, Value* removed);
  //bool removeMember(String const& key, Value* removed);
  bool removeMember(const char* begin, const char* end, Value* removed);
  bool removeIndex(ArrayIndex index, Value* removed);
  bool isMember(const char* key) const;
  //bool isMember(const String& key) const;
  bool isMember(const char* begin, const char* end) const;

  Members getMemberNames() const;
  
  //void setComment(const char* comment, CommentPlacement placement);  
  //void setComment(const char* comment, size_t len, CommentPlacement placement);  
  //void setComment(String comment, CommentPlacement placement);
  //bool hasComment(CommentPlacement placement) const;  
  //String getComment(CommentPlacement placement) const;
  //String toStyledString() const;

  //const_iterator begin() const;
  //const_iterator end() const;

  //iterator begin();
  //iterator end();

  void setOffsetStart(ptrdiff_t start);
  void setOffsetLimit(ptrdiff_t limit);
  ptrdiff_t getOffsetStart() const;
  ptrdiff_t getOffsetLimit() const;

};

class Reader {
public:
  using Char = char;
  using Location = const Char*;

  struct StructuredError;
  
  Reader();
  Reader(const Features& features);

  bool parse(const std::string& document, Value& root,
             bool collectComments = true);

  bool parse(const char* beginDoc, const char* endDoc, Value& root,
             bool collectComments = true);

  //bool parse(IStream& is, Value& root, bool collectComments = true);

  String getFormatedErrorMessages() const;
  String getFormattedErrorMessages() const;
  std::vector<StructuredError> getStructuredErrors() const;
  bool pushError(const Value& value, const String& message);
  bool pushError(const Value& value, const String& message, const Value& extra);
  bool good() const;
}; 

struct Reader::StructuredError {
    ptrdiff_t offset_start;
    ptrdiff_t offset_limit;
    String message;
};

class FastWriter {
public:
  
  FastWriter();
  ~FastWriter();

  void enableYAMLCompatibility();

  void dropNullPlaceholders();
  void omitEndingLineFeed();

  String write(const Value& root) override;
};

} // namespace Json


