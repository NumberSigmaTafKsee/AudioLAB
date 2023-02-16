// $Id: DateTime.hh,v 1.1 2005/12/17 15:50:11 sonarnerd Exp $

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


#ifndef DATETIME_HH
#define DATETIME_HH


#include <cstddef>
#include <ctime>

#include <commoncpp/Exception.hh>
#include <commoncpp/Compilers.hh>
#include <commoncpp/String.hh>

#define DATETIME_MAX_CONV_LENGTH        255


typedef double datetime_t;  ///< Internal time representation


datetime_t DateTime_Now ();
datetime_t DateTime_HighResNow ();
datetime_t DateTime_RelativeNow ();
datetime_t DateTime_Monotonic ();
datetime_t DateTime_ProcessNow ();
datetime_t DateTime_ThreadNow ();


    /**
        Exception class for clDateTime.
    */
    class clXDateTime : public clException
    {
        public:
            clXDateTime (const char *cpErrorMsg) :
                clException(cpErrorMsg)
                { }
            clXDateTime (const char *cpErrorMsg, int iErrorCode) :
                clException(cpErrorMsg, iErrorCode)
                { }
    };


    /**
        Class for date and time handling.
    */
    class clDateTime
    {
        private:

            void Init ();

        protected:

            datetime_t dTime;   ///< Internal time storage
            struct tm sTM;      ///< Temporary storage for struct tm types

        public:

            /**
                Constructor.
            */
            clDateTime ()
                {
                    Init();
                    dTime = 0.0;
                }
            /**
                Copy constructor.
            */
            clDateTime (const clDateTime &CopyTime)
                { this->operator=(CopyTime); }
            /// \overload
            clDateTime (const datetime_t &dNewTime)
                { this->operator=(dNewTime); }
            /// \overload
            clDateTime (const time_t &timeNew)
                { this->operator=(timeNew); }
            /// \overload
            clDateTime (struct tm *spNewTime)
                { this->operator=(spNewTime); }
            /**
                Destructor.
            */
            virtual ~clDateTime ()
                { }

            /**
                Get current time.
                
                \throw clXDateTime
            */
            static clDateTime Now ()
                { return clDateTime(DateTime_Now()); }
            /**
                Get current time in high resolution format.
                
                \throw clXDateTime
            */
            static clDateTime HighResNow ()
                { return clDateTime(DateTime_HighResNow()); }
            /**
                Get current time in high resolution format, only for delta use.
                
                \throw clXDateTime
            */
            static clDateTime RelativeNow ()
                { return clDateTime(DateTime_RelativeNow()); }
            /**
                Get current time from monotonic clock.
                
                \throw clXDateTime
            */
            static clDateTime Monotonic ()
                { return clDateTime(DateTime_Monotonic()); }
            /**
                Get time spent by current process.
                
                \throw clXDateTime
            */
            static clDateTime ProcessNow ()
                { return clDateTime(DateTime_ProcessNow()); }
            /**
                Get time spent by current thread.
                
                \throw clXDateTime
            */
            static clDateTime ThreadNow ()
                { return clDateTime(DateTime_ThreadNow()); }
            /**
                Update internal time storage to current time.
                
                \throw clXDateTime
            */
            void Update ()
                { dTime = DateTime_Now(); }

            /**
                Assignment operator.
            */
            clDateTime & operator = (const clDateTime &NewTime);
            /// \overload
            clDateTime & operator = (const datetime_t &dNewTime)
                {
                    Init();
                    dTime = dNewTime;
                    return (*this);
                }
            /// \overload
            clDateTime & operator = (const time_t &timeNew)
                {
                    Init();
                    dTime = (datetime_t) timeNew;
                    return (*this);
                }
            /**
                \overload
                
                \throw clXDateTime
            */
            clDateTime & operator = (struct tm *);

            clDateTime & operator += (const clDateTime &Op2)
                {
                    dTime += Op2.dTime;
                    return (*this);
                }
            clDateTime & operator -= (const clDateTime &Op2)
                {
                    dTime -= Op2.dTime;
                    return (*this);
                }
            clDateTime & operator *= (const clDateTime &Op2)
                {
                    dTime *= Op2.dTime;
                    return (*this);
                }
            clDateTime & operator /= (const clDateTime &Op2)
                {
                    dTime /= Op2.dTime;
                    return (*this);
                }
            clDateTime operator + (const clDateTime &Op2) const
                {
                    clDateTime Res(dTime + Op2.dTime);
                    return Res;
                }
            clDateTime operator - (const clDateTime &Op2) const
                {
                    clDateTime Res(dTime - Op2.dTime);
                    return Res;
                }
            clDateTime operator * (const clDateTime &Op2) const
                {
                    clDateTime Res(dTime * Op2.dTime);
                    return Res;
                }
            clDateTime operator / (const clDateTime &Op2) const
                {
                    clDateTime Res(dTime / Op2.dTime);
                    return Res;
                }

            bool operator == (const clDateTime &Op2) const
                { return (dTime == Op2.dTime); }
            bool operator != (const clDateTime &Op2) const
                { return (dTime != Op2.dTime); }
            bool operator <= (const clDateTime &Op2) const
                { return (dTime <= Op2.dTime); }
            bool operator >= (const clDateTime &Op2) const
                { return (dTime >= Op2.dTime); }
            bool operator < (const clDateTime &Op2) const
                { return (dTime < Op2.dTime); }
            bool operator > (const clDateTime &Op2) const
                { return (dTime > Op2.dTime); }
            bool operator () (const clDateTime &Op1,
                const clDateTime &Op2) const
                { return (Op1 < Op2); }

            /**
                Get absolute time difference.
                
                \param Op2 Second time value
                \return Absolute time difference
            */
            datetime_t Diff (const clDateTime &) const;

            /**
                Get internal "struct tm *".
                
                \return Internal "struct tm *"
            */
            struct tm * tm_ptr ()
                { return (&sTM); }
            /**
                Get time as "time_t".
                
                \return Time as "time_t"
            */
            time_t Time ()
                { return ((time_t) (dTime + 0.5)); }
            /**
                Get time as local "struct tm *" time.
                
                \note Also sets internal "struct tm".
                
                \return Time as local "struct tm *"
                \throw clXDateTime
            */
            struct tm * LocalTime ();
            /**
                Get time as UTC "struct tm *" time.
                
                \note Also sets internal "struct tm".
                
                \return Time as UTC "struct tm *"
                \throw clXDateTime
            */
            struct tm * GmTime ();

            // you have to call LocalTime() or GmTime() first

            /**
                Get day part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Days
            */
            int Day () const
                { return sTM.tm_mday; }
            /**
                Set day part of internal "struct tm".
                
                \param iDay Days
            */
            void Day (int iDay)
                { sTM.tm_mday = iDay; }
            /**
                Get month part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Months
            */
            int Month () const
                { return sTM.tm_mon + 1; }
            /**
                Set month part of internal "struct tm".
                
                \param iMonth Monts
            */
            void Month (int iMonth)
                { sTM.tm_mon = iMonth - 1; }
            /**
                Get year part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Years
            */
            int Year () const
                { return sTM.tm_year + 1900; }
            /**
                Set year part of internal "struct tm".
                
                \param iYear Years
            */
            void Year (int iYear)
                { sTM.tm_year = iYear - 1900; }
            /**
                Get hour part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Hours
            */
            int Hour () const
                { return sTM.tm_hour; }
            /**
                Set hour part of internal "struct tm".
                
                \param iHour Hours
            */
            void Hour (int iHour)
                { sTM.tm_hour = iHour; }
            /**
                Get minute part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Minutes
            */
            int Minute () const
                { return sTM.tm_min; }
            /**
                Set minute part of internal "struct tm".
                
                \param iMinute Minutes
            */
            void Minute (int iMinute)
                { sTM.tm_min = iMinute; }
            /**
                Get seconds part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Seconds
            */
            int Second () const
                { return sTM.tm_sec; }
            /**
                Set seconds part of internal "struct tm".
                
                \param iSecond Seconds
            */
            void Second (int iSecond)
                { sTM.tm_sec = iSecond; }
            /**
                Get milliseconds part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Milliseconds
            */
            int Millisecond () const;
            /**
                Get microseconds part of internal "struct tm".
                
                \note You have to call LocalTime() or GmTime() first!
                
                \return Microseconds
            */
            int Microsecond () const;

            /**
                Get time as "double".
                
                This is similar in representation to time_t, except having
                subsecond resolution in fractional part.
                
                \return Time as "double"
            */
            double AsDouble () const
                { return dTime; }
            /**
                Get local time as string.
                
                \param cpFormat Format as for strftime()
                \return Time string
                \throw clXDateTime
            */
            clString AsString (const char *);
            /**
                Get UTC time as string.
                
                \param cpFormat Format as for strftime()
                \return Time string
                \throw clXDateTime
            */
            clString GmAsString (const char *);
            /**
                Parse time string.
                
                \param cpParse Time string
                \param cpFormat Format as for strptime()
                \throw clXDateTime
            */
            void Parse (const char *, const char *);
    };


#endif
