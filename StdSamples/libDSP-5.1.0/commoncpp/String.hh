// $Id: String.hh,v 1.1 2005/12/17 15:50:11 sonarnerd Exp $

/*

    Class for string handling
    Copyright (C) 2005 Jussi Laako

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/


#ifndef STRING_HH
#define STRING_HH

#include <cstddef>
#include <string>

#include <commoncpp/Exception.hh>
#include <commoncpp/Compilers.hh>

#define STRING_MAX_ERROR_LENGTH         255


    /**
        Exception class for clString.
    */
    class clXString : public clException
    {
        public:
            clXString (const char *cpErrorMsg) :
                clException(cpErrorMsg)
                { }
            clXString (const char *cpErrorMsg, int iErrorCode) :
                clException(cpErrorMsg, iErrorCode)
                { }
    };


    /**
        Class for string handling.
    */
    class clString
    {
            void Init ()
                {
                    sizeAlloc = 0;
                    sizeStr = 0;
                    cpStr = NULL;
                    cpTokStr = NULL;
                    cpTokTmp = NULL;
                }
            size_t RoundSize (size_t);
            void CheckSize (size_t);
            void SetSize (size_t);
            int Compare (const char *) const;

        protected:

            size_t sizeAlloc;   ///< Allocation size
            size_t sizeStr;     ///< String size
            char *cpStr;        ///< String pointer
            char *cpTokStr;     ///< Temporary pointer for Tok
            char *cpTokTmp;     ///< Temporary pointer for Tok

        public:

            /**
                Constructor.
                
                \throw clXString
            */
            clString ();
            /**
                Copy constructor.
                
                \throw clXString
            */
            clString (const clString &CopySrc)
                {
                    Init();
                    *this = CopySrc;
                }
            /// \overload
            clString (const char *, size_t);
            /// \overload
            clString (const char *cpNewStr)
                {
                    Init();
                    this->operator=(cpNewStr);
                }
            /// \overload
            clString (const std::string strNew)
                {
                    Init();
                    this->operator=(strNew.c_str());
                }
            /// \overload
            clString (char cValue)
                {
                    Init();
                    Int((long) cValue);
                }
            /// \overload
            clString (unsigned char ucValue)
                {
                    Init();
                    UInt((size_t) ucValue);
                }
            /// \overload
            clString (short iValue)
                {
                    Init();
                    Int((long) iValue);
                }
            /// \overload
            clString (unsigned short uValue)
                {
                    Init();
                    UInt((size_t) uValue);
                }
            /// \overload
            clString (int iValue)
                {
                    Init();
                    Int((long) iValue);
                }
            /// \overload
            clString (unsigned int uValue)
                {
                    Init();
                    UInt((size_t) uValue);
                }
            /// \overload
            clString (long lValue)
                {
                    Init();
                    Int(lValue);
                }
            /// \overload
            clString (unsigned long ulValue)
                {
                    Init();
                    UInt(ulValue);
                }
            /// \overload
            clString (float fValue)
                {
                    Init();
                    Float((double) fValue);
                }
            /// \overload
            clString (double dValue)
                {
                    Init();
                    Float(dValue);
                }
            /**
                Destructor.
            */
            virtual ~clString ();

            /**
                Set string from integer.
                
                \param lValue Value
            */
            void Int (long);
            // \overload
            void UInt (unsigned long);
            /**
                Set string from float.
                
                \param dValue Value
            */
            void Float (double);
            /**
                Get string as C "const char *".
                
                \return Pointer to C string
            */
            const char *c_str () const
                { return cpStr; }
            /**
                Get length of string.
                
                \return Length
            */
            size_t Length () const
                { return sizeStr; }
            /**
                Get size of string (including terminatin NULL).
                
                \return Size
            */
            size_t Size () const
                { return (sizeStr + 1); }
            /**
                Check if string is empty.
                
                \return String is empty?
            */
            bool IsEmpty () const
                { return ((sizeStr == 0) ? true : false); }
            /**
                Clear (empty) string.
            */
            void Clear ();

            /**
                Assignment operator.
            */
            clString & operator = (const char *);
            /// \overload
            clString & operator = (const clString &);
            clString & operator += (const char *);
            clString & operator += (const clString &);
            clString operator + (const char *);
            clString operator + (const clString &);

            operator const char * () const
                { return cpStr; }
            /**
                Access specified character in string.
                
                \param sizeIdx Index to character
                \throw clXString
            */
            char & operator [] (size_t);

            /**
                Case sensitive compare.
            */
            bool operator == (const clString &CmpStr) const
                {
                    if (Compare(CmpStr.c_str()) == 0) return true;
                    else return false;
                }
            bool operator == (const char *cpCmpStr) const
                {
                    if (Compare(cpCmpStr) == 0) return true;
                    else return false;
                }
            bool operator != (const clString &CmpStr) const
                {
                    if (Compare(CmpStr.c_str()) != 0) return true;
                    else return false;
                }
            bool operator != (const char *cpCmpStr) const
                {
                    if (Compare(cpCmpStr) != 0) return true;
                    else return false;
                }
            bool operator < (const clString &CmpStr) const
                {
                    if (Compare(CmpStr.c_str()) < 0) return true;
                    else return false;
                }
            bool operator < (const char *cpCmpStr) const
                {
                    if (Compare(cpCmpStr) < 0) return true;
                    else return false;
                }
            bool operator > (const clString &CmpStr) const
                {
                    if (Compare(CmpStr.c_str()) > 0) return true;
                    else return false;
                }
            bool operator > (const char *cpCmpStr) const
                {
                    if (Compare(cpCmpStr) > 0) return true;
                    else return false;
                }
            bool operator <= (const clString &CmpStr) const
                {
                    if (Compare(CmpStr.c_str()) <= 0) return true;
                    else return false;
                }
            bool operator <= (const char *cpCmpStr) const
                {
                    if (Compare(cpCmpStr) <= 0) return true;
                    else return false;
                }
            bool operator >= (const clString &CmpStr) const
                {
                    if (Compare(CmpStr.c_str()) >= 0) return true;
                    else return false;
                }
            bool operator >= (const char *cpCmpStr) const
                {
                    if (Compare(cpCmpStr) >= 0) return true;
                    else return false;
                }
            bool operator () (const clString &CmpStr1,
                const clString &CmpStr2) const
                {
                    return (CmpStr1 < CmpStr2);
                }

            /**
                Case insensitive compare.
                
                This compare is otherwise equal to ::operator== .
            */
            int CaseCompare (const clString &CmpStr) const
                { return CaseCompare(CmpStr.c_str()); }
            int CaseCompare (const char *) const;

            /**
                Case sensitively compare specified range in strings.
                
                \param CmpStr String to compare against
                \param uiIdx Start index
                \param uiLen Length of compare
                \return Result [-1, 0, 1]
                \throw clXString
            */
            int CompareRange (const clString &CmpStr,
                unsigned int uiIdx, unsigned int uiLen) const
                { return CompareRange(CmpStr.c_str(), uiIdx, uiLen); }
            /// \overload
            int CompareRange (const char *, unsigned int, unsigned int) const;
            /**
                Case insensitively compare specified range in strings.
                
                \param CmpStr String to compare against
                \param uiIdx Start index
                \param uiLen Length of compare
                \return Result [-1, 0, 1]
                \throw clXString
            */
            int CaseCompareRange (const clString &CmpStr,
                unsigned int uiIdx, unsigned int uiLen) const
                { return CaseCompareRange(CmpStr.c_str(), uiIdx, uiLen); }
            /// \overload
            int CaseCompareRange (const char *, unsigned int,
                unsigned int) const;
            /**
                Collate compare.
                
                \note See strcoll(3)
                
                \param String to compare against
                \return Result [-1, 0, 1]
            */
            int CollCompare (const clString &CmpStr) const
                { return CollCompare(CmpStr.c_str()); }
            /// \overload
            int CollCompare (const char *) const;
            /**
                Find first occurrence of character in string.
                
                \param cChar Character
                \return Position
            */
            long Pos (char) const;
            /**
                Find last occurrence of character in string.
                
                \param cChar Character
                \return Position
            */
            long RPos (char) const;
            /**
                Find first occurrence of string within string.
                
                \param PosStr String to search for
                \return Resulting string
            */
            clString Pos (clString &PosStr) const
                { return Pos(PosStr.c_str()); }
            /// \overload
            clString Pos (const char *cpPosStr) const;
            /// \overload
            const char * PosPtr (const char *cpPosStr) const;
            /**
                Scan string for a charater.
                
                \param BrkChr Characters to search for
                \return Position
            */
            long Break (const clString &BrkChr) const
                { return Break(BrkChr.c_str()); }
            /// \overload
            long Break (const char *) const;
            /**
                Get length of substring.
                
                \note See strspn(3)
                
                \param SpnChr Characters to search for
                \return Length
            */
            long Spn (const clString &SpnChr) const
                { return Spn(SpnChr.c_str()); }
            /// \overload
            long Spn (const char *) const;
            /**
                Get length of complementary substring.
                
                \note See strcspn(3)
                
                \param SpnChr Characters to search for
                \return Length
            */
            long CSpn (const clString &SpnChr) const
                { return CSpn(SpnChr.c_str()); }
            /// \overload
            long CSpn (const char *) const;
            /**
                Split string into tokens.
                
                \param TokChr Delimiter characters
                \param bFirst Start from the beginning
                \return Token
                \throw clXString
            */
            clString Tok (const clString &TokChr, bool bFirst = false)
                { return Tok(TokChr.c_str(), bFirst); }
            /// \overload
            clString Tok (const char *, bool = false);
            /**
                End Tok().
                
                This frees up the temporary storage allocated for Tok().
            */
            void TokEnd ();

            /**
                Get substring.
                
                \param uIdx Start index
                \return Substring
                \throw clXString
            */
            clString SubString (unsigned int) const;
            /**
                Get substring.
                
                \param uIdx Start index
                \param uCount Length of substring
                \return Substring
                \throw clXString
            */
            clString SubString (unsigned int, unsigned int) const;
            /**
                Insert string into the string.
                
                \param InsStr String to insert
                \param uIdx Index of insert position
                \throw clXString
            */
            void Insert (const clString &, unsigned int);
            /**
                Delete characters.
                
                \param uIdx Index of start position
                \param uCount Number of characters to delete
                \throw clXString
            */
            void Delete (unsigned int, unsigned int);

            /**
                Convert string to lower case.
            */
            void ToLower ();
            /**
                Convert string to upper case.
            */
            void ToUpper ();
            /**
                Remove leading and trailing whitespace.
            */
            void Clean ();

            /**
                Convert string to integer value.
                
                \return Integer value
                \throw clXString
            */
            long ToInt ();
            /// \overload
            unsigned long ToUInt ();
            /**
                Convert string to floating point value.
                
                \return Floating point value
                \throw clXString
            */
            float ToFloat ();
            /// \overload
            double ToDouble ();
            /// \overload
            long double ToLongDouble ();

            /**
                Set string from errno(3) value.
            */
            void StrError (int);
    };


#endif
