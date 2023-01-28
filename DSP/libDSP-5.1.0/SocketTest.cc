// $Id: SocketTest.cc,v 1.2 2005/12/17 15:50:11 sonarnerd Exp $

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "commoncpp/SocketAddr.hh"
#include "commoncpp/SocketOp.hh"


int main (int argc, char *argv[])
{
    try
    {
        clSocketAddrLocal AddrL;
        clSocketAddrIPv4 Addr4;
        clSocketAddrIPv6 Addr6;

        puts("--- UNIX");

        AddrL.Set("/tmp/test.socket");
        puts(AddrL.AsString());

        puts("--- IPv4");

        Addr4.Set("localhost", 1234);
        puts(Addr4.AsStringP());
        puts(Addr4.AsHostNameP());

        puts("--- IPv6");

        Addr6.Set("ipv6-localhost", 1234);
        puts(Addr6.AsStringP());
        puts(Addr6.AsHostNameP());

        puts("-----");

        puts(clSocketAddrDetect::Get(AddrL.Get(), AddrL.Size())->AsString());
        puts(clSocketAddrDetect::Get(Addr4.Get(), Addr4.Size())->AsString());
        puts(clSocketAddrDetect::Get(Addr6.Get(), Addr6.Size())->AsString());
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
