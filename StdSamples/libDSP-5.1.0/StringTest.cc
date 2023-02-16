// $Id: StringTest.cc,v 1.2 2005/12/17 15:50:11 sonarnerd Exp $

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>

#include "commoncpp/String.hh"


bool fp_comp (double d1, double d2)
{
    if (fabs(d1 - d2) < DBL_EPSILON)
        return true;
    return false;
}


int main (int argc, char *argv[])
{
    try
    {
        clString Hello("Hello world!");
        if (strcmp(Hello.c_str(), "Hello world!") == 0)
        {
            printf("%s: OK\n", Hello.c_str());
        }
        else
        {
            puts("Not OK!");
            return 1;
        }

        clString Value(-1048576);
        if (atoi(Value.c_str()) == -1048576)
        {
            printf("%s: OK\n", Value.c_str());
        }
        else
        {
            puts(".Int() not OK!");
            return 1;
        }
        
        clString Value2(1048576u);
        if (atoi(Value2.c_str()) == 1048576)
        {
            printf("%s: OK\n", Value2.c_str());
        }
        else
        {
            puts(".UInt() not OK!");
            return 1;
        }

        clString Hello1("Hello ");
        clString Hello2("world!");
        clString Addr1("127.0.0.1");
        clString Addr2(1234u);
        if (strcmp((Hello1 + Hello2).c_str(), "Hello world!") == 0 &&
            (strcmp((Addr1 + ":" + Addr2).c_str(), "127.0.0.1:1234") == 0))
        {
            puts("operator+ OK");
        }
        else
        {
            puts("operator+ not OK!");
            return 1;
        }
        Hello1 += Hello2;
        if (strcmp(Hello1.c_str(), "Hello world!") == 0)
        {
            puts("operator+= OK");
        }
        else
        {
            puts("operator+= not OK!");
            return 1;
        }

        clString SubStr("---|XXX|---");
        clString SubRes;
        SubRes = SubStr.SubString(4, 3);
        if (strcmp(SubRes.c_str(), "XXX") == 0)
        {
            puts(".SubString() OK");
        }
        else
        {
            puts(".SubString() not OK!");
            return 1;
        }

        clString Insert("===||===");
        Insert.Insert(clString("###"), 4);
        if (strcmp(Insert.c_str(), "===|###|===") == 0)
        {
            puts(".Insert() OK");
        }
        else
        {
            puts(".Insert() not OK!");
            return 1;
        }

        clString Delete("abcdef");
        Delete.Delete(4, 2);
        Delete.Delete(0, 2);
        if (strcmp(Delete.c_str(), "cd") == 0)
        {
            puts(".Delete() OK");
        }
        else
        {
            puts(".Delete() not OK!");
            return 1;
        }

        clString ToLower("AbCdEf");
        clString ToUpper("AbCdEf");
        ToLower.ToLower();
        ToUpper.ToUpper();
        if (strcmp(ToLower.c_str(), "abcdef") == 0)
        {
            puts(".ToLower() OK");
        }
        else
        {
            puts(".ToLower() not OK!");
            return 1;
        }
        if (strcmp(ToUpper.c_str(), "ABCDEF") == 0)
        {
            puts(".ToUpper() OK");
        }
        
        clString Clean1("\t   ...---XXX---...\t   \n\r");
        clString Clean2("   \t...---XXX---...");
        clString Clean3("...---XXX---...   \t\n\r");
        clString Clean4("XXX");
        Clean1.Clean();
        Clean2.Clean();
        Clean3.Clean();
        Clean4.Clean();
        if ((strcmp(Clean1.c_str(), "...---XXX---...") == 0) &&
            (strcmp(Clean2.c_str(), "...---XXX---...") == 0) &&
            (strcmp(Clean3.c_str(), "...---XXX---...") == 0) &&
            (strcmp(Clean4.c_str(), "XXX") == 0))
        {
            puts(".Clean() OK");
        }
        else
        {
            puts(".Clean() not OK!");
            return 1;
        }
        
        clString Pos("...XXX...");
        if (Pos.Pos('X') == 3)
        {
            puts(".Pos() OK");
        }
        else
        {
            puts(".Pos() not OK!");
            return 1;
        }
        if (Pos.RPos('X') == 5)
        {
            puts(".RPos() OK");
        }
        else
        {
            puts(".RPos() not OK!");
            return 1;
        }
        clString PosRes = Pos.Pos("XXX");
        if (strcmp(PosRes.c_str(), "XXX...") == 0)
        {
            puts(".Pos() OK");
        }
        else
        {
            puts(".Pos() not OK!");
            return 1;
        }
        
        clString Int("12345");
        clString Float("3.14e-3");
        if (Int.ToInt() == 12345)
        {
            puts(".ToInt() OK");
        }
        else
        {
            puts(".ToInt() not OK!");
            return 1;
        }
        if (fp_comp(Float.ToDouble(), 3.14e-3))
        {
            puts(".ToFloat() OK");
        }
        else
        {
            puts(".ToFloat() not OK!");
            return 1;
        }

        clString Tok("\t    Hello \t, world");
        clString Tok1, Tok2;
        Tok1 = Tok.Tok(" \t\n\r");
        Tok2 = Tok.Tok(", \t\n\r");
        if (strcmp(Tok1.c_str(), "Hello") == 0 &&
            strcmp(Tok2.c_str(), "world") == 0)
        {
            puts(".Tok() OK");
        }
        else
        {
            puts(".Tok() not OK!");
            return 1;
        }

        clString StrError;
        StrError.StrError(0);
        puts(StrError.c_str());
    }
    catch (clException &X)
    {
        fprintf(stderr, "exception: %s\n", X.what());
        return 10;
    }
    catch (...)
    {
        fprintf(stderr, "unknown exception\n");
        return 20;
    }
    
    return 0;
}
