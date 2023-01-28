// $Id: SocketAddr.cc,v 1.3 2005/12/17 15:50:11 sonarnerd Exp $

/*

    Classes for address type abstraction
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

#include <arpa/inet.h>
#include <netdb.h>

#include "commoncpp/SocketAddr.hh"


// --- clSocketAddrLocal


void clSocketAddrLocal::Init ()
{
    memset(&sAddr, 0x00, sizeof(sAddr));
}


void clSocketAddrLocal::Set (const char *cpFileName)
{
    /*if (strlen(cpFileName) >= UNIX_PATH_MAX)
        throw clXSocketAddr("clSocketAddrLocal::Set(): socket path too long");*/
    sAddr.sun_family = AF_UNIX;
    strncpy(sAddr.sun_path, cpFileName, UNIX_PATH_MAX - 1);
    sAddr.sun_path[UNIX_PATH_MAX - 1] = '\0';
}


void clSocketAddrLocal::Set (sockaddrptr_t spNewAddr, socklen_t uiAddrSize)
{
    if (spNewAddr == NULL)
        throw clXSocketAddr("clSocketAddrLocal::Set(): spNewAddr == NULL");
    if (uiAddrSize != sizeof(sAddr))
        throw clXSocketAddr("clSocketAddrLocal::Set(): size mismatch");
    memcpy(&sAddr, spNewAddr, sizeof(sAddr));
}


clString clSocketAddrLocal::AsString () const
{
    return clString(sAddr.sun_path);
}


clSocketAddrLocal & clSocketAddrLocal::operator = (
    const clSocketAddrLocal &Src)
{
    memcpy(&sAddr, &Src.sAddr, sizeof(sAddr));
    return (*this);
}


// --- clSocketAddrIPv4


void clSocketAddrIPv4::Init ()
{
    memset(&sAddr, 0x00, sizeof(sAddr));
}


void clSocketAddrIPv4::Set (const char *cpAddress, unsigned short uiPort)
{
    int iEC;
    struct addrinfo sAddrHint;
    struct addrinfo *spAddrInfo;

    sAddr.sin_family = AF_INET;
    sAddr.sin_port = htons(uiPort);
    if (cpAddress != NULL)
    {
        iEC = inet_pton(AF_INET, cpAddress, &sAddr.sin_addr.s_addr);
        if (iEC < 0)
            throw clXSocketAddr("clSocketAddrIPv4::Set(): inet_pton()", errno);
        else if (iEC == 0)
        {
            memset(&sAddrHint, 0x00, sizeof(sAddrHint));
            sAddrHint.ai_family = AF_INET;
            if (getaddrinfo(cpAddress, NULL, &sAddrHint, &spAddrInfo) != 0)
                throw clXSocketAddr("clSocketAddrIPv4::Set(): getaddrinfo()",
                    errno);
            sAddr.sin_addr.s_addr = 
                ((struct sockaddr_in *) spAddrInfo->ai_addr)->sin_addr.s_addr;
            freeaddrinfo(spAddrInfo);
        }
    }
    else
    {
        sAddr.sin_addr.s_addr = htonl(INADDR_ANY);
    }
}


void clSocketAddrIPv4::Set (sockaddrptr_t spNewAddr, socklen_t uiAddrSize)
{
    if (spNewAddr == NULL)
        throw clXSocketAddr("clSocketAddrIPv4::Set(): spNewAddr == NULL");
    if (uiAddrSize != sizeof(sAddr))
        throw clXSocketAddr("clSocketAddrIPv4::Set(): size mismatch");
    memcpy(&sAddr, spNewAddr, sizeof(sAddr));
}


clString clSocketAddrIPv4::AsString () const
{
    char cpAddrStr[INET_ADDRSTRLEN + 1];

    if (inet_ntop(AF_INET, &sAddr.sin_addr.s_addr,
        cpAddrStr, INET_ADDRSTRLEN) == NULL)
        throw clXSocketAddr("clSocketAddrIPv4::AsString(): inet_ntop()");
    return clString(cpAddrStr);
}


clString clSocketAddrIPv4::AsStringP () const
{
    clString PortStr(ntohs(sAddr.sin_port));
    return (AsString() + ":" + PortStr);
}


clString clSocketAddrIPv4::AsHostName () const
{
    char cpHostName[SOCKETADDR_MAX_HOSTNAME_LENGTH + 1];

    if (getnameinfo((struct sockaddr *) &sAddr, (socklen_t) sizeof(sAddr),
        cpHostName, SOCKETADDR_MAX_HOSTNAME_LENGTH, NULL, 0, 0) != 0)
        throw clXSocketAddr("clSocketAddrIPv4::AsHostName(): getnameinfo()",
            errno);
    return clString(cpHostName);
}


clString clSocketAddrIPv4::AsHostNameP () const
{
    clString PortStr(ntohs(sAddr.sin_port));
    return (AsHostName() + ":" + PortStr);
}


struct ip_mreq * clSocketAddrIPv4::GetMcastReq (const clSocketAddrIPv4 &Addr,
    const clSocketAddrIPv4 *pIf)
{
    const struct sockaddr_in *spAddr;
    const struct sockaddr_in *spIf = NULL;

    memset(&sMreq, 0x00, sizeof(sMreq));
    spAddr = (const struct sockaddr_in *) Addr.Get();
    memcpy(&sMreq.imr_multiaddr, &spAddr->sin_addr, sizeof(struct in_addr));
    if (pIf != NULL)
    {
        spIf = (const struct sockaddr_in *) pIf->Get();
        memcpy(&sMreq.imr_interface, &spIf->sin_addr, sizeof(struct in_addr));
    }
    else
    {
        sMreq.imr_interface.s_addr = htonl(INADDR_ANY);
    }
    return (&sMreq);
}


clSocketAddrIPv4 & clSocketAddrIPv4::operator = (
    const clSocketAddrIPv4 &Src)
{
    memcpy(&sAddr, &Src.sAddr, sizeof(sAddr));
    return (*this);
}


// --- clSocketAddrIPv6


void clSocketAddrIPv6::Init ()
{
    memset(&sAddr, 0x00, sizeof(sAddr));
}


void clSocketAddrIPv6::Set (const char *cpAddress, unsigned short uiPort)
{
    int iEC;
    struct addrinfo sAddrHint;
    struct addrinfo *spAddrInfo;

    sAddr.sin6_family = AF_INET6;
    sAddr.sin6_port = htons(uiPort);
    if (cpAddress != NULL)
    {
        iEC = inet_pton(AF_INET6, cpAddress, &sAddr.sin6_addr.in6_u);
        if (iEC < 0)
            throw clXSocketAddr("clSocketAddrIPv6::Set(): inet_pton()", errno);
        else if (iEC == 0)
        {
            memset(&sAddrHint, 0x00, sizeof(sAddrHint));
            sAddrHint.ai_family = AF_INET6;
            if (getaddrinfo(cpAddress, NULL, &sAddrHint, &spAddrInfo) != 0)
                throw clXSocketAddr("clSocketAddrIPv6::Set(): getaddrinfo()",
                    errno);
            memcpy(&sAddr.sin6_addr,
                &(((struct sockaddr_in6 *) spAddrInfo->ai_addr)->sin6_addr),
                sizeof(struct in6_addr));
            freeaddrinfo(spAddrInfo);
        }
    }
    else
    {
        memcpy(&sAddr.sin6_addr, &in6addr_any, sizeof(struct in6_addr));
    }
}


void clSocketAddrIPv6::Set (sockaddrptr_t spNewAddr, socklen_t uiAddrSize)
{
    if (spNewAddr == NULL)
        throw clXSocketAddr("clSocketAddrIPv6::Set(): spNewAddr == NULL");
    if (uiAddrSize != sizeof(sAddr))
        throw clXSocketAddr("clSocketAddrIPv6::Set(): size mismatch");
    memcpy(&sAddr, spNewAddr, sizeof(sAddr));
}


clString clSocketAddrIPv6::AsString () const
{
    char cpAddrStr[INET6_ADDRSTRLEN + 1];

    if (inet_ntop(AF_INET6, &sAddr.sin6_addr.in6_u,
        cpAddrStr, INET6_ADDRSTRLEN) == NULL)
        throw clXSocketAddr("clSocketAddrIPv6::AsString(): inet_ntop()");
    return clString(cpAddrStr);
}


clString clSocketAddrIPv6::AsStringP () const
{
    clString PortStr(ntohs(sAddr.sin6_port));
    return (AsString() + ":" + PortStr);
}


clString clSocketAddrIPv6::AsHostName () const
{
    char cpHostName[SOCKETADDR_MAX_HOSTNAME_LENGTH + 1];

    if (getnameinfo((struct sockaddr *) &sAddr, (socklen_t) sizeof(sAddr),
        cpHostName, SOCKETADDR_MAX_HOSTNAME_LENGTH, NULL, 0, 0) != 0)
        throw clXSocketAddr("clSocketAddrIPv6::AsHostName(): getnameinfo()",
            errno);
    return clString(cpHostName);
}


clString clSocketAddrIPv6::AsHostNameP () const
{
    clString PortStr(ntohs(sAddr.sin6_port));
    return (AsHostName() + ":" + PortStr);
}


struct ipv6_mreq * clSocketAddrIPv6::GetMcastReq (const clSocketAddrIPv6 &Addr,
    unsigned int uiIf)
{
    const struct sockaddr_in6 *spAddr;

    memset(&sMreq, 0x00, sizeof(sMreq));
    spAddr = (const struct sockaddr_in6 *) Addr.Get();
    memcpy(&sMreq.ipv6mr_multiaddr, &spAddr->sin6_addr,
        sizeof(struct in6_addr));
    sMreq.ipv6mr_interface = uiIf;
    return (&sMreq);
}


clSocketAddrIPv6 & clSocketAddrIPv6::operator = (
    const clSocketAddrIPv6 &Src)
{
    memcpy(&sAddr, &Src.sAddr, sizeof(sAddr));
    return (*this);
}


// --- clSocketAddrDetect


clSocketAddr * clSocketAddrDetect::Get (sockaddrptr_t spAddr,
    socklen_t uiAddrSize)
{
    if (uiAddrSize == sizeof(struct sockaddr_un))
        return ((clSocketAddr *) new clSocketAddrLocal(spAddr, uiAddrSize));
    else if (uiAddrSize == sizeof(struct sockaddr_in))
        return ((clSocketAddr *) new clSocketAddrIPv4(spAddr, uiAddrSize));
    else if (uiAddrSize == sizeof(struct sockaddr_in6))
        return ((clSocketAddr *) new clSocketAddrIPv6(spAddr, uiAddrSize));
    throw clXSocketAddr("clSocketAddrDetect::GetAddr(): no match");
}
