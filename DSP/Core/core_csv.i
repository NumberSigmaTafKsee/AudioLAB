%module csv
%{
#include "rapidcsv.h"
using namespace rapidcsv;
%}
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"

namespace rapidcsv
{

    struct ConverterParams
    {

        explicit ConverterParams(const bool pHasDefaultConverter = false,
                                const long double pDefaultFloat = std::numeric_limits<long double>::signaling_NaN(),
                                const long long pDefaultInteger = 0,
                                const bool pNumericLocale = true);

        bool mHasDefaultConverter;
        long double mDefaultFloat;
        long long mDefaultInteger;
        bool mNumericLocale;
    };


    template<typename T>
    class Converter
    {
    public:

        Converter(const ConverterParams& pConverterParams);

        void ToStr(const T& pVal, std::string& pStr) const;                
        void ToVal(const std::string& pStr, T& pVal) const;            
    };

    struct LabelParams
    {

        explicit LabelParams(const int pColumnNameIdx = 0, const int pRowNameIdx = -1);

        int mColumnNameIdx;
        int mRowNameIdx;
    };

    struct SeparatorParams
    {

        explicit SeparatorParams(const char pSeparator = ',', const bool pTrim = false,
                                const bool pHasCR = sPlatformHasCR, const bool pQuotedLinebreaks = false,
                                const bool pAutoQuote = true);

        char mSeparator;
        bool mTrim;
        bool mHasCR;
        bool mQuotedLinebreaks;
        bool mAutoQuote;
    };

    struct LineReaderParams
    {
        explicit LineReaderParams(const bool pSkipCommentLines = false,
                                const char pCommentPrefix = '#',
                                const bool pSkipEmptyLines = false);

        bool mSkipCommentLines;
        char mCommentPrefix;
        bool mSkipEmptyLines;
    };

    class Document
    {
    public:
        
        explicit Document(const std::string& pPath = std::string(),
                        const LabelParams& pLabelParams = LabelParams(),
                        const SeparatorParams& pSeparatorParams = SeparatorParams(),
                        const ConverterParams& pConverterParams = ConverterParams(),
                        const LineReaderParams& pLineReaderParams = LineReaderParams());

        explicit Document(std::istream& pStream,
                        const LabelParams& pLabelParams = LabelParams(),
                        const SeparatorParams& pSeparatorParams = SeparatorParams(),
                        const ConverterParams& pConverterParams = ConverterParams(),
                        const LineReaderParams& pLineReaderParams = LineReaderParams());

        void Load(const std::string& pPath,
                const LabelParams& pLabelParams = LabelParams(),
                const SeparatorParams& pSeparatorParams = SeparatorParams(),
                const ConverterParams& pConverterParams = ConverterParams(),
                const LineReaderParams& pLineReaderParams = LineReaderParams());

        void Load(std::istream& pStream,
                const LabelParams& pLabelParams = LabelParams(),
                const SeparatorParams& pSeparatorParams = SeparatorParams(),
                const ConverterParams& pConverterParams = ConverterParams(),
                const LineReaderParams& pLineReaderParams = LineReaderParams());

        void Save(const std::string& pPath = std::string());
        void Save(std::ostream& pStream);
        void Clear();

        size_t GetColumnIdx(const std::string& pColumnName) const;
        
        %extend {
            std::vector<std::string> GetColumn(const size_t p) { return $self->GetColumn<std::string>(p); }
            std::vector<std::string> GetColumn(const std::string& p) { return $self->GetColumn<std::string>(p); }
        }

        void SetColumn(const size_t pColumnIdx, const std::vector<std::string>& pColumn);
        void SetColumn(const std::string& pColumnName, const std::vector<std::string>& pColumn);
        void RemoveColumn(const size_t pColumnIdx);
        void RemoveColumn(const std::string& pColumnName);

        %extend {
            void InsertColumn(const size_t pColumnIdx, const std::vector<std::string>& pColumn = std::vector<std::string>(),
                        const std::string& pColumnName = std::string()) { return $self->InsertColumn<std::string>(pColumnIdx, pColumn, pColumnName); }
        }

        size_t GetRowIdx(const std::string& pRowName) const;
        
        %extend {
            std::vector<std::string> GetRow(const size_t pRowIdx) const { return $self->GetRow<std::string>(pRowIdx); }
            std::vector<std::string> GetRow(const std::string& pRowName) const { return $self->GetRow<std::string>(pRowName); }
        }
        void SetRow(const size_t pRowIdx, const std::vector<std::string>& pRow);
        void SetRow(const std::string& pRowName, const std::vector<std::string>& pRow);

        void RemoveRow(const size_t pRowIdx);
        void RemoveRow(const std::string& pRowName);


        %extend  {
            void InsertRow(const size_t pRowIdx, const std::vector<std::string>& pRow = std::vector<std::string>(),const std::string& pRowName = std::string())
            {
                $self->InsertRow<std::string>(pRowIdx,pRow,pRowName);
            }
        }        
        size_t GetRowCount() const;
        
        %extend {
            std::string GetCell(const size_t pColumnIdx, const size_t pRowIdx) const { return $self->GetCell<std::string>(pColumnIdx, pRowIdx); }
            std::string GetCell(const std::string& pColumnName, const std::string& pRowName) const { return $self->GetCell<std::string>(pColumnName, pRowName); }
            std::string GetCell(const std::string& pColumnName, const size_t pRowIdx) const { return $self->GetCell<std::string>(pColumnName, pRowIdx); }
            std::string GetCell(const size_t pColumnIdx, const std::string& pRowName) const { return $self->GetCell<std::string>(pColumnIdx, pRowName); }
        }
        void SetCell(const size_t pColumnIdx, const size_t pRowIdx, const std::string& pCell);
        void SetCell(const std::string& pColumnName, const std::string& pRowName, const std::string& pCell);

        std::string GetColumnName(const size_t pColumnIdx);
        void SetColumnName(size_t pColumnIdx, const std::string& pColumnName);
        
        std::string GetRowName(const size_t pRowIdx);
        void SetRowName(size_t pRowIdx, const std::string& pRowName);

        %extend {
            std::vector<std::string> GetColumnNames() { return $self->GetColumnNames(); }
            std::vector<std::string> GetRowNames() { return $self->GetRowNames(); }
        }
    };
}

%template(string_vector) std::vector<std::string>;