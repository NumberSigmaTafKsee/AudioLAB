// $Id: DateTime.cc,v 1.2 2005/12/17 15:50:11 sonarnerd Exp $

/*

    Class for time handling
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


#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <cfloat>

#include <sys/time.h>

#include "commoncpp/DateTime.hh"


datetime_t DateTime_Now ()
{
    datetime_t time;
    struct timeval timevalNow;

    memset(&timevalNow, 0x00, sizeof(timevalNow));
    if (unlikely(gettimeofday(&timevalNow, NULL) != 0))
        throw clXDateTime("DateTime_Now(): gettimeofday()", errno);
    time = (datetime_t) timevalNow.tv_sec +
        (datetime_t) timevalNow.tv_usec * 1.0e-6;
    return time;
}


datetime_t DateTime_HighResNow ()
{
    datetime_t time;
    struct timespec timespecNow;

    if (unlikely(clock_gettime(CLOCK_REALTIME, &timespecNow) != 0))
        throw clXDateTime("DateTime_HighResNow(): clock_gettime()", errno);
    time = (datetime_t) timespecNow.tv_sec +
        (datetime_t) timespecNow.tv_nsec * 1.0e-9;
    return time;
}


datetime_t DateTime_RelativeNow ()
{
    datetime_t time;
    struct timespec timespecNow;

    if (unlikely(clock_gettime(CLOCK_REALTIME, &timespecNow) != 0))
        throw clXDateTime("DateTime_RelativeNow(): clock_gettime()", errno);
    // relative clock will wrap around once per 30 days
    time = (datetime_t) (timespecNow.tv_sec % 2592000) +
        (datetime_t) timespecNow.tv_nsec * 1.0e-9;
    return time;
}


datetime_t DateTime_Monotonic ()
{
    datetime_t time;
    struct timespec timespecNow;

    if (unlikely(clock_gettime(CLOCK_MONOTONIC, &timespecNow) != 0))
        throw clXDateTime("DateTime_Monotonic(): clock_gettime()", errno);
    time = (datetime_t) timespecNow.tv_sec +
        (datetime_t) timespecNow.tv_nsec * 1.0e-9;
    return time;
}


datetime_t DateTime_ProcessNow ()
{
    datetime_t time;
    clockid_t clockidClock;
    struct timespec timespecNow;

    if (unlikely(clock_getcpuclockid(0, &clockidClock) != 0))
        throw clXDateTime("DateTime_ProcessNow(): clock_getcpuclockid()",
            errno);
    if (unlikely(clock_gettime(clockidClock, &timespecNow) != 0))
        throw clXDateTime("DateTime_ProcessNow(): clock_gettime()", errno);
    time = (datetime_t) timespecNow.tv_sec +
        (datetime_t) timespecNow.tv_nsec * 1.0e-9;
    return time;
}


datetime_t DateTime_ThreadNow ()
{
    datetime_t time;
    clockid_t clockidClock;
    struct timespec timespecNow;

    if (unlikely(pthread_getcpuclockid(pthread_self(), &clockidClock) != 0))
        throw clXDateTime("DateTime_ThreadNow(): pthread_getcpuclockid()",
            errno);
    if (unlikely(clock_gettime(clockidClock, &timespecNow) != 0))
        throw clXDateTime("DateTime_ThreadNow(): clock_gettime()", errno);
    time = (datetime_t) timespecNow.tv_sec +
        (datetime_t) timespecNow.tv_nsec * 1.0e-9;
    return time;
}


void clDateTime::Init ()
{
    memset(&sTM, 0x00, sizeof(sTM));
}


clDateTime & clDateTime::operator = (const clDateTime &NewTime)
{
    dTime = NewTime.dTime;
    memcpy(&sTM, &NewTime.sTM, sizeof(struct tm));
    return (*this);
}


clDateTime & clDateTime::operator = (struct tm *spNewTime)
{
    time_t timeNew;

    Init();
    if (unlikely(spNewTime == NULL))
        throw clXDateTime(
            "clDateTime::operator=(struct tm*): spNewTime == NULL");
    memcpy(&sTM, spNewTime, sizeof(struct tm));
    timeNew = mktime(spNewTime);
    if (timeNew == (time_t) -1)
        throw clXDateTime("clDateTime::operator=(struct tm *): mktime()",
            errno);
    this->operator=(timeNew);
    return (*this);
}


datetime_t clDateTime::Diff (const clDateTime &Op2) const
{
    datetime_t diff;

    diff = fabs(dTime - Op2.dTime);
    return diff;
}


struct tm * clDateTime::LocalTime ()
{
    time_t time;

    time = Time();
    if (unlikely(localtime_r(&time, &sTM) != &sTM))
        throw clXDateTime("clDateTime::LocalTime(): localtime_r()", errno);
    return (&sTM);
}


struct tm * clDateTime::GmTime ()
{
    time_t time;

    time = Time();
    if (unlikely(gmtime_r(&time, &sTM) != &sTM))
        throw clXDateTime("clDateTime::GmTime(): gmtime_r()", errno);
    return (&sTM);
}


int clDateTime::Millisecond () const
{
    int iMSecs;

    iMSecs = (int) ((dTime - floor(dTime)) * 1000.0);
    return iMSecs;
}


int clDateTime::Microsecond () const
{
    int iUSecs;

    iUSecs = (int) ((dTime - floor(dTime)) * 1.0e6);
    return iUSecs;
}


clString clDateTime::AsString (const char *cpFormat)
{
    char cpConvBuf[DATETIME_MAX_CONV_LENGTH + 1];

    if (unlikely(cpFormat == NULL))
        throw clXDateTime("clDateTime::AsString(): cpFormat == NULL");
    if (unlikely(strftime(cpConvBuf, DATETIME_MAX_CONV_LENGTH,
        cpFormat, LocalTime()) == 0))
        throw clXDateTime("clDateTime::AsString(): strftime()", errno);
    return clString(cpConvBuf);
}


clString clDateTime::GmAsString (const char *cpFormat)
{
    char cpConvBuf[DATETIME_MAX_CONV_LENGTH + 1];

    if (unlikely(cpFormat == NULL))
        throw clXDateTime("clDateTime::GmAsString(): cpFormat == NULL");
    if (unlikely(strftime(cpConvBuf, DATETIME_MAX_CONV_LENGTH,
        cpFormat, GmTime()) == 0))
        throw clXDateTime("clDateTime::GmAsString(): strftime()", errno);
    return clString(cpConvBuf);
}


void clDateTime::Parse (const char *cpParse, const char *cpFormat)
{
    struct tm sTMParse;

    if (unlikely(cpParse == NULL))
        throw clXDateTime("clDateTime::Parse(): cpParse == NULL");
    if (unlikely(cpFormat == NULL))
        throw clXDateTime("clDateTime::Parse(): cpFormat == NULL");
    if (unlikely(strptime(cpParse, cpFormat, &sTMParse) == NULL))
        throw clXDateTime("clDateTime::Parse(): strptime()", errno);
    this->operator=(&sTMParse);
}
