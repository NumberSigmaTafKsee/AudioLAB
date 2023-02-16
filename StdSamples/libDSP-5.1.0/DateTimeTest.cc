// $Id: DateTimeTest.cc,v 1.2 2005/12/17 15:50:11 sonarnerd Exp $

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "commoncpp/DateTime.hh"


int main (int argc, char *argv[])
{
    try
    {
        datetime_t now;
        datetime_t highres;
        datetime_t relative;
        datetime_t monotonic;
        clDateTime TimeNow;

        now = DateTime_Now();
        highres = DateTime_HighResNow();
        relative = DateTime_RelativeNow();
        monotonic = DateTime_Monotonic();
        TimeNow = clDateTime::Now();
        printf("DateTime_Now() = %g\n", now);
        printf("DateTime_HighResNow() = %g\n", highres);
        printf("DateTime_RelativeNow() = %g\n", relative);
        printf("DateTime_Monotonic() = %g\n", monotonic);
        printf("highresnow / now delta: %g us\n", fabs(highres - now) * 1.0e6);

        printf("DateTime_ProcessNow() = %g us\n",
            DateTime_ProcessNow() * 1.0e6);
        printf("DateTime_ThreadNow() = %g us\n",
            DateTime_ThreadNow() * 1.0e6);

        clDateTime StringTime;
        clDateTime StringTime2;
        clString TimeStr;
        StringTime.Update();
        TimeStr = StringTime.AsString("%d.%m.%Y %T");
        puts(TimeStr);
        StringTime2.Parse(TimeStr, "%d.%m.%Y %T");
        if (StringTime.Time() == StringTime2.Time())
        {
            puts(".AsString() && .Parse() OK");
        }
        else
        {
            puts(".AsString() || .Parse() not OK!");
            return 1;
        }
    }
    catch (clException &X)
    {
        puts(X.what());
        return 1;
    }
    catch (...)
    {
        puts("unknown exception");
        return 1;
    }

    return 0;
}
