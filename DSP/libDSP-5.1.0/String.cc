// $Id: String.cc,v 1.4 2005/12/17 15:50:11 sonarnerd Exp $

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


#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cerrno>
#include <cctype>

// strcasecmp() & co
#include <strings.h>

#include "commoncpp/String.hh"


#define STRING_CONVBUF_SIZE     32


inline size_t clString::RoundSize (size_t sizeExact)
{
    size_t sizeRound;

    sizeRound = (size_t) pow(2, ceil(log(sizeExact + 1) / log(2)));
    if (unlikely(sizeRound < sizeof(void *)))
        sizeRound = sizeof(void *);
    return sizeRound;
}


inline void clString::CheckSize (size_t sizeNew)
{
    if (sizeNew >= sizeAlloc)
    {
        SetSize(RoundSize(sizeNew));
    }
    else if (sizeNew < (sizeAlloc / 2))
    {
        SetSize(RoundSize(sizeNew));
    }
}


void clString::SetSize (size_t sizeNew)
{
    sizeAlloc = sizeNew;
    if (cpStr != NULL)
    {
        cpStr = (char *) realloc(cpStr, sizeAlloc);
        if (unlikely(cpStr == NULL))
            throw clXString("clString::CheckSize(): realloc()", errno);
    }
    else
    {
        cpStr = (char *) malloc(sizeAlloc);
        if (unlikely(cpStr == NULL))
            throw clXString("clString::CheckSize(): malloc()", errno);
    }
}


int clString::Compare (const char *cpCmpStr) const
{
    return strcmp(cpStr, cpCmpStr);
}


clString::clString ()
{
    Init();
    sizeAlloc = sizeof(void *);
    cpStr = (char *) malloc(sizeAlloc);
    if (unlikely(cpStr == NULL))
        throw clXString("clString::clString(): malloc()", errno);
    sizeStr = 0;
    cpStr[0] = '\0';
}


clString::clString (const char *cpNewStr, size_t sizeNew)
{
    Init();
    CheckSize(sizeNew);
    if (cpNewStr != NULL)
    {
        strncpy(cpStr, cpNewStr, sizeNew);
        cpStr[sizeNew - 1] = '\0';
        sizeStr = strlen(cpStr);
    }
}


clString::~clString ()
{
    if (likely(cpStr != NULL))
        free(cpStr);
}


void clString::Int (long lValue)
{
    size_t sizeValue;

    // value + possible sign
    sizeValue = (size_t) ceil(log10(abs(lValue + 1))) + 1;
    CheckSize(sizeValue);
    snprintf(cpStr, sizeValue + 1, "%ld", lValue);
    sizeStr = sizeValue;
}


void clString::UInt (unsigned long ulValue)
{
    size_t sizeValue;

    sizeValue = (size_t) ceil(log10(ulValue + 1));
    CheckSize(sizeValue);
    snprintf(cpStr, sizeValue + 1, "%lu", ulValue);
    sizeStr = sizeValue;
}


void clString::Float (double dValue)
{
    size_t sizeValue;

    // sign + mantissa + decimal dot + 'e' + sign + exponent == 24
    /*sizeValue = (size_t) 1 + 
        ceil(log10(pow(2, 52))) + 1 + 2 + ceil(log10(pow(2, 11)));*/
    sizeValue = 32;
    SetSize(sizeValue);
    snprintf(cpStr, sizeValue, "%g", dValue);
    sizeStr = strlen(cpStr);
}


void clString::Clear ()
{
    SetSize(sizeof(void *));
    sizeStr = 0;
    cpStr[0] = '\0';
}


clString & clString::operator = (const char *cpNewStr)
{
    sizeStr = strlen(cpNewStr);
    CheckSize(sizeStr);
    strcpy(cpStr, cpNewStr);
    return (*this);
}


clString & clString::operator = (const clString &NewStr)
{
    SetSize(NewStr.sizeAlloc);
    sizeStr = NewStr.sizeStr;
    strcpy(cpStr, NewStr.cpStr);
    return (*this);
}


clString & clString::operator += (const char *cpAddStr)
{
    size_t sizeRes;
    size_t sizeAdd;

    sizeAdd = strlen(cpAddStr);
    sizeRes = sizeStr + sizeAdd;
    CheckSize(sizeRes);
    /*strcat(cpStr, cpAddStr);*/
    // optimization, as we already know lengths of both strings
    memcpy(&cpStr[sizeStr], cpAddStr, sizeAdd + 1);
    sizeStr = sizeRes;
    return (*this);
}


clString & clString::operator += (const clString &AddStr)
{
    size_t sizeNew;

    sizeNew = sizeStr + AddStr.sizeStr;
    CheckSize(sizeNew);
    /*strcat(cpStr, AddStr.cpStr);*/
    // optimization, as we already know lengths of both strings
    memcpy(&cpStr[sizeStr], AddStr.cpStr, AddStr.sizeStr + 1);
    sizeStr = sizeNew;
    return (*this);
}


clString clString::operator + (const char *cpAddStr)
{
    size_t sizeRes;
    size_t sizeAdd;

    sizeAdd = strlen(cpAddStr);
    sizeRes = sizeStr + sizeAdd;

    clString Res(NULL, sizeRes);

    Res.sizeStr = sizeRes;
    /*strcpy(Res.cpStr, cpStr);
    strcat(Res.cpStr, cpAddStr);*/
    // optimization, as we already know lengths
    memcpy(Res.cpStr, cpStr, sizeStr);
    memcpy(&Res.cpStr[sizeStr], cpAddStr, sizeAdd + 1);
    return Res;
}


clString clString::operator + (const clString &AddStr)
{
    size_t sizeRes;

    sizeRes = sizeStr + AddStr.sizeStr;

    clString Res(NULL, sizeRes);
    
    Res.sizeStr = sizeRes;
    /*strcpy(Res.cpStr, cpStr);
    strcat(Res.cpStr, AddStr.cpStr);*/
    // optimization, as we already know lengths
    memcpy(Res.cpStr, cpStr, sizeStr);
    memcpy(&Res.cpStr[sizeStr], AddStr.cpStr, AddStr.sizeStr + 1);
    return Res;
}


char & clString::operator [] (size_t sizeIdx)
{
    if (sizeIdx >= sizeStr)
        throw clXString("clString::operator[](): sizeIdx >= sizeStr");
    return cpStr[sizeIdx];
}


int clString::CaseCompare (const char *cpCmpStr) const
{
    return strcasecmp(cpStr, cpCmpStr);
}


int clString::CompareRange (const char *cpCmpStr, unsigned int uiIdx,
    unsigned int uiLen) const
{
    if (uiIdx + uiLen > sizeStr)
        throw clXString("clString::CompareRange(): uiIdx + uiLen > sizeStr");
    return strncmp(cpStr + uiIdx, cpCmpStr, uiLen);
}


int clString::CaseCompareRange (const char *cpCmpStr, unsigned int uiIdx,
    unsigned int uiLen) const
{
    if (uiIdx + uiLen > sizeStr)
        throw clXString(
            "clString::CaseCompareRange(): uiIdx + uiLen > sizeStr");
    return strncasecmp(cpStr + uiIdx, cpCmpStr, uiLen);
}


int clString::CollCompare (const char *cpCmpStr) const
{
    return strcoll(cpStr, cpCmpStr);
}


long clString::Pos (char cChar) const
{
    char *cpPos;

    cpPos = strchr(cpStr, cChar);
    if (unlikely(cpPos == NULL))
        return -1;
    return (long) (cpPos - cpStr);
}


long clString::RPos (char cChar) const
{
    char *cpPos;

    cpPos = strrchr(cpStr, cChar);
    if (unlikely(cpPos == NULL))
        return -1;
    return (long) (cpPos - cpStr);
}


clString clString::Pos (const char *cpPosStr) const
{
    char *cpPos;
    clString Res;

    cpPos = strstr(cpStr, cpPosStr);
    if (likely(cpPos != NULL))
        Res = cpPos;
    return Res;
}


const char * clString::PosPtr (const char *cpPosStr) const
{
    return strstr(cpStr, cpPosStr);
}


long clString::Break (const char *cpBrkChr) const
{
    char *cpPos;

    cpPos = strpbrk(cpStr, cpBrkChr);
    if (unlikely(cpPos == NULL))
        return -1;
    return (long) (cpPos - cpStr);
}


long clString::Spn (const char *cpSpnChr) const
{
    long lRes = -1;

    lRes = (long) strspn(cpStr, cpSpnChr);
    return lRes;
}


long clString::CSpn (const char *cpSpnChr) const
{
    long lRes = -1;

    lRes = (long) strcspn(cpStr, cpSpnChr);
    return lRes;
}


clString clString::Tok (const char *cpTokChr, bool bFirst)
{
    char *cpRes;
    clString Res;

    if (bFirst)
    {
        if (cpTokStr != NULL)
            TokEnd();
        cpTokStr = strdup(cpStr);
        if (cpTokStr == NULL)
            throw clXString("clString::Tok(): strdup()", errno);
        cpRes = strtok_r(cpTokStr, cpTokChr, &cpTokTmp);
        if (cpRes != NULL)
            Res = cpRes;
    }
    else
    {
        if (cpTokStr == NULL)
        {
            cpTokStr = strdup(cpStr);
            if (cpTokStr == NULL)
                throw clXString("clString::Tok(): strdup()", errno);
            cpRes = strtok_r(cpTokStr, cpTokChr, &cpTokTmp);
        }
        else
        {
            cpRes = strtok_r(NULL, cpTokChr, &cpTokTmp);
        }
        if (cpRes != NULL)
            Res = cpRes;
    }
    return Res;
}


void clString::TokEnd ()
{
    if (cpTokStr != NULL)
    {
        free(cpTokStr);
        cpTokStr = NULL;
    }
    cpTokTmp = NULL;
}


clString clString::SubString (unsigned int uIdx) const
{
    unsigned int uiCount;

    if (unlikely(uIdx >= sizeStr))
        throw clXString("clString::SubString(): uIdx >= sizeStr");
    uiCount = sizeStr - uIdx;
    return SubString(uIdx, uiCount);
}


clString clString::SubString (unsigned int uIdx, unsigned int uCount) const
{
    clString Res;

    if (unlikely(uIdx >= sizeStr))
        throw clXString("clString::SubString(): uIdx >= sizeStr");
    if (unlikely((uIdx + uCount) > sizeStr))
        throw clXString("clString::SubString(): (uIdx + uCount) > sizeStr");
    Res.CheckSize(uCount);
    Res.sizeStr = uCount;
    memcpy(Res.cpStr, cpStr + uIdx, uCount);
    Res.cpStr[uCount] = '\0';
    return Res;
}


void clString::Insert (const clString &InsStr, unsigned int uIdx)
{
    if (unlikely(uIdx > sizeStr))
        throw clXString("clString::Insert(): uIdx > sizeStr");
    CheckSize(sizeStr + InsStr.Length());
    // take the NULL char into account, thus +1
    memmove(cpStr + uIdx + InsStr.Length(), cpStr + uIdx,
        sizeStr - uIdx + 1);
    memcpy(cpStr + uIdx, InsStr.c_str(), InsStr.Length());
    sizeStr += InsStr.Length();
}


void clString::Delete (unsigned int uIdx, unsigned int uCount)
{
    if (unlikely(uIdx >= sizeStr))
        throw clXString("clString::Delete(): uIdx >= sizeStr");
    if (unlikely((uIdx + uCount) > sizeStr))
        throw clXString("clString::Delete(): (uIdx + uCount) > sizeStr");
    // take the NULL char into account, thus +1
    memmove(cpStr + uIdx, cpStr + uIdx + uCount,
        (sizeStr + 1) - (uIdx + uCount));
    sizeStr -= uCount;
    CheckSize(sizeStr);
}


void clString::ToLower ()
{
    char *cpTemp = cpStr;

    while (*cpTemp != '\0')
    {
        *cpTemp = tolower(*cpTemp);
        cpTemp++;
    }
}


void clString::ToUpper ()
{
    char *cpTemp = cpStr;

    while (*cpTemp != '\0')
    {
        *cpTemp = toupper(*cpTemp);
        cpTemp++;
    }
}


void clString::Clean ()
{
    char *cpBegin = cpStr;
    char *cpEnd;
    size_t sizeRes;

    while (*cpBegin != '\0')
    {
        if (!isblank(*cpBegin))
            break;
        cpBegin++;
    }
    cpEnd = cpBegin;
    while (*cpEnd != '\0')
    {
        if (isblank(*cpEnd))
            break;
        cpEnd++;
    }
    sizeRes = cpEnd - cpBegin;
    if (likely(sizeRes > 0))
        memmove(cpStr, cpBegin, sizeRes);
    cpStr[sizeRes] = '\0';
    sizeStr = sizeRes;
    CheckSize(sizeStr);
}


long clString::ToInt ()
{
    long lRes;
    char *cpEnd;

    lRes = strtol(cpStr, &cpEnd, 0);
    if (likely((cpEnd != cpStr) && (*cpEnd == '\0')))
        return lRes;
    throw clXString("clString::ToInt(): not an integer");
}


unsigned long clString::ToUInt ()
{
    unsigned long ulRes;
    char *cpEnd;

    ulRes = strtoul(cpStr, &cpEnd, 0);
    if (likely((cpEnd != cpStr) && (*cpEnd == '\0')))
        return ulRes;
    throw clXString("clString::ToUInt(): not an integer");
}


float clString::ToFloat ()
{
    float fRes;
    char *cpEnd;

    fRes = strtof(cpStr, &cpEnd);
    if (likely((cpEnd != cpStr) && (*cpEnd == '\0')))
        return fRes;
    throw clXString("clString::ToFloat(): not a float");
}


double clString::ToDouble ()
{
    double dRes;
    char *cpEnd;

    dRes = strtod(cpStr, &cpEnd);
    if (likely((cpEnd != cpStr) && (*cpEnd == '\0')))
        return dRes;
    throw clXString("clString::ToDouble(): not a double");
}


long double clString::ToLongDouble ()
{
    long double ldRes;
    char *cpEnd;

    ldRes = strtold(cpStr, &cpEnd);
    if (likely((cpEnd != cpStr) && (*cpEnd == '\0')))
        return ldRes;
    throw clXString("clString::ToLongDouble(): not a long double");
}


void clString::StrError (int iErrorCode)
{
    /* strerror_r() seems to be buggy _still_ */
    /*char cpTemp[STRING_MAX_ERROR_LENGTH + 1];

    if (strerror_r(iErrorCode, cpTemp, STRING_MAX_ERROR_LENGTH) != 0)
        throw clXString("clString::StrError()");
    cpTemp[STRING_MAX_ERROR_LENGTH] = '\0';
    *this = cpTemp;*/

    char *cpError;

    cpError = strerror(iErrorCode);
    if (cpError == NULL)
        throw clXString("clString::StrError()");
    *this = cpError;
}
