// $Id: SocketOp.cc,v 1.6 2005/12/17 15:50:11 sonarnerd Exp $

/*

    Class for socket handling
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

#include <fcntl.h>
#include <sys/select.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netinet/udp.h>
#include <arpa/inet.h>

#include "commoncpp/String.hh"
#include "commoncpp/SocketOp.hh"


bool clSocket::GetLine (clString &RecvString)
{
    long lNLIdx = -1;
    char *cpRecvBuf;
    size_t uIdx;
    
    cpRecvBuf = RecvBuf;
    for (uIdx = 0; uIdx < ulRecvCount; uIdx++)
    {
        if (cpRecvBuf[uIdx] == '\n' || cpRecvBuf[uIdx] == '\0')
        {
            lNLIdx = (long) uIdx;
            break;
        }
    }
    if (lNLIdx < 0)
        return false;
    
    cpRecvBuf[lNLIdx] = '\0';
    RecvString = cpRecvBuf;
    
    return true;
}


clSocket::~clSocket ()
{
    try
    {
        if (bCloseOnDestruct)
            Close();
    }
    catch (...)
    { }
}


void clSocket::Create (int iAddrFamily, int iType, int iProtocol)
{
    iHandle = socket(iAddrFamily, iType, iProtocol);
    if (SOCKET_INVALID_HANDLE(iHandle))
        throw clXSocket("clSocket::Create(): socket()", errno);
}


socketh_t clSocket::CreatePair (int iType)
{
    socketh_t ipSockH[2];

    if (socketpair(AF_UNIX, iType, 0, ipSockH) != 0)
        throw clXSocket("clSocket::CreatePair(): socketpair()", errno);
    iHandle = ipSockH[0];
    return ipSockH[1];
}


void clSocket::Bind (const clSocketAddr &Addr)
{
    if (bind(iHandle, Addr.Get(), Addr.Size()) != 0)
        throw clXSocket("clSocket::Bind(): bind()", errno);
}


void clSocket::Connect (const clSocketAddr &Addr)
{
    if (connect(iHandle, Addr.Get(), Addr.Size()) != 0)
        throw clXSocket("clSocket::Connect(): connect()", errno);
}


bool clSocket::Connect (const clSocketAddr &Addr, unsigned int uiTimeout)
{
    bool bWaitRes;

    SetBlocking(false);
    if (connect(iHandle, Addr.Get(), Addr.Size()) != 0)
    {
        if (errno != EINPROGRESS)
            throw clXSocket("clSocket::Connect(): connect()", errno);
    }
    bWaitRes = WaitRead(uiTimeout);
    SetBlocking(true);
    if (!bWaitRes)
        return false;
    return true;
}


void clSocket::Listen (int iBacklog)
{
    if (listen(iHandle, iBacklog) != 0)
        throw clXSocket("clSocket::Listen(): listen()", errno);
}


socketh_t clSocket::Accept (sockaddrptr_t spAddr, socklen_t *uipAddrSize)
{
    socketh_t iSockH;

    iSockH = accept(iHandle, spAddr, uipAddrSize);
    if (iSockH < 0)
        throw clXSocket("clSocket::Accept(): accept()", errno);
    return iSockH;
}


socketh_t clSocket::Accept (clSocketAddr &Addr)
{
    socketh_t iSockH;
    socklen_t iAddrSize;

    iAddrSize = Addr.Size();
    iSockH = accept(iHandle, Addr.Get(), &iAddrSize);
    if (iSockH < 0)
        throw clXSocket("clSocket::Accept(): accept()", errno);
    return iSockH;
}


socketh_t clSocket::Accept (clSocketAddr **ppAddr)
{
    socketh_t iSockH;
    // we try to reserve enough space for any possible address
    unsigned char cpAddrBuf[SOCKETOP_MAX_ADDR_SIZE];
    socklen_t uiAddrSize;

    uiAddrSize = SOCKETOP_MAX_ADDR_SIZE;
    memset(cpAddrBuf, 0x00, sizeof(cpAddrBuf));
    iSockH = accept(iHandle, (struct sockaddr *) cpAddrBuf, &uiAddrSize);
    if (iSockH < 0)
        clXSocket("clSocket::Accept(): accept()", errno);
    if (ppAddr != NULL)
    {
        *ppAddr =
            clSocketAddrDetect::Get((sockaddrptr_t) cpAddrBuf, uiAddrSize);
    }
    return iSockH;
}


bool clSocket::WaitRead (unsigned int uiTimeout)
{
    int iSelectRes;
    fd_set fdsetRead;
    fd_set fdsetError;
    struct timeval sTimeout;

    FD_ZERO(&fdsetRead);
    FD_SET(iHandle, &fdsetRead);
    FD_ZERO(&fdsetError);
    FD_SET(iHandle, &fdsetError);
    sTimeout.tv_sec = uiTimeout / 1000u;
    sTimeout.tv_usec = uiTimeout % 1000u * 1000u;
    iSelectRes = select(iHandle + 1, &fdsetRead, NULL, &fdsetError, &sTimeout);
    if (iSelectRes == 0)
        return false;
    else if (iSelectRes < 0)
        throw clXSocket("clSocket::WaitRead(): select()", errno);
    if (FD_ISSET(iHandle, &fdsetError))
        throw clXSocketError("clSocket::WaitRead(): socket error");
    return true;
}


bool clSocket::WaitWrite (unsigned int uiTimeout)
{
    int iSelectRes;
    fd_set fdsetWrite;
    fd_set fdsetError;
    struct timeval sTimeout;

    FD_ZERO(&fdsetWrite);
    FD_SET(iHandle, &fdsetWrite);
    FD_ZERO(&fdsetError);
    FD_SET(iHandle, &fdsetError);
    sTimeout.tv_sec = uiTimeout / 1000u;
    sTimeout.tv_usec = uiTimeout % 1000u * 1000u;
    iSelectRes = select(iHandle + 1, NULL, &fdsetWrite, &fdsetError, &sTimeout);
    if (iSelectRes == 0)
        return false;
    else if (iSelectRes < 0)
        throw clXSocket("clSocket::WaitWrite(): select()", errno);
    if (FD_ISSET(iHandle, &fdsetError))
        throw clXSocketError("clSocket::WaitWrite(): socket error");
    return true;
}


bool clSocket::WaitReadWrite (unsigned int uiTimeout)
{
    int iSelectRes;
    fd_set fdsetRead;
    fd_set fdsetWrite;
    fd_set fdsetError;
    struct timeval sTimeout;

    FD_ZERO(&fdsetRead);
    FD_SET(iHandle, &fdsetRead);
    FD_ZERO(&fdsetWrite);
    FD_SET(iHandle, &fdsetWrite);
    FD_ZERO(&fdsetError);
    FD_SET(iHandle, &fdsetError);
    sTimeout.tv_sec = uiTimeout / 1000u;
    sTimeout.tv_usec = uiTimeout % 1000u * 1000u;
    iSelectRes = select(iHandle + 1, &fdsetRead, &fdsetWrite, &fdsetError,
        &sTimeout);
    if (iSelectRes == 0)
        return false;
    else if (iSelectRes < 0)
        throw clXSocket("clSocket::WaitReadWrite(): select()", errno);
    if (FD_ISSET(iHandle, &fdsetError))
        throw clXSocketError("clSocket::WaitReadWrite(): socket error");
    return true;
}


unsigned long clSocket::WaitRead (const socketh_t *ipSockH,
    unsigned int uiHandleCount, unsigned int uiTimeout)
{
    int iSelectRes;
    unsigned int uiIdx;
    unsigned long ulRet = 0;
    socketh_t iSockH;
    socketh_t iMaxH = 0;
    fd_set fdsetRead;
    fd_set fdsetError;
    struct timeval sTimeout;

    if (ipSockH == NULL)
        throw clXSocket("clSocketOp::WaitRead(): ipSockH == NULL");
    FD_ZERO(&fdsetRead);
    FD_ZERO(&fdsetError);
    for (uiIdx = 0; uiIdx < uiHandleCount; uiIdx++)
    {
        iSockH = ipSockH[uiIdx];
        FD_SET(iSockH, &fdsetRead);
        FD_SET(iSockH, &fdsetError);
        if (iSockH > iMaxH)
            iMaxH = iSockH;
    }
    sTimeout.tv_sec = uiTimeout / 1000u;
    sTimeout.tv_usec = uiTimeout % 1000u * 1000u;
    iSelectRes = select(iMaxH + 1, &fdsetRead, NULL, &fdsetError, &sTimeout);
    if (iSelectRes == 0)
        return 0ul;
    else if (iSelectRes < 0)
        throw clXSocket("clSocket::WaitRead(): select()", errno);
    for (uiIdx = 0; uiIdx < uiHandleCount; uiIdx++)
    {
        iSockH = ipSockH[uiIdx];
        if (FD_ISSET(iSockH, &fdsetError))
            throw clXSocketError("clSocket::WaitRead(): socket error");
        if (FD_ISSET(iSockH, &fdsetRead))
            ulRet |= (1 << uiIdx);
    }
    return ulRet;
}


unsigned long clSocket::WaitWrite (const socketh_t *ipSockH,
    unsigned int uiHandleCount, unsigned int uiTimeout)
{
    int iSelectRes;
    unsigned int uiIdx;
    unsigned long ulRet = 0;
    socketh_t iSockH;
    socketh_t iMaxH = 0;
    fd_set fdsetWrite;
    fd_set fdsetError;
    struct timeval sTimeout;

    if (ipSockH == NULL)
        throw clXSocket("clSocketOp::WaitWrite(): ipSockH == NULL");
    FD_ZERO(&fdsetWrite);
    FD_ZERO(&fdsetError);
    for (uiIdx = 0; uiIdx < uiHandleCount; uiIdx++)
    {
        iSockH = ipSockH[uiIdx];
        FD_SET(iSockH, &fdsetWrite);
        FD_SET(iSockH, &fdsetError);
        if (iSockH > iMaxH)
            iMaxH = iSockH;
    }
    sTimeout.tv_sec = uiTimeout / 1000u;
    sTimeout.tv_usec = uiTimeout % 1000u * 1000u;
    iSelectRes = select(iMaxH + 1, NULL, &fdsetWrite, &fdsetError, &sTimeout);
    if (iSelectRes == 0)
        return 0ul;
    else if (iSelectRes < 0)
        throw clXSocket("clSocket::WaitWrite(): select()", errno);
    for (uiIdx = 0; uiIdx < uiHandleCount; uiIdx++)
    {
        iSockH = ipSockH[uiIdx];
        if (FD_ISSET(iSockH, &fdsetError))
            throw clXSocketError("clSocket::WaitWrite(): socket error");
        if (FD_ISSET(iSockH, &fdsetWrite))
            ulRet |= (1 << uiIdx);
    }
    return ulRet;
}


unsigned long clSocket::WaitReadWrite (const socketh_t *ipSockH,
    unsigned int uiHandleCount, unsigned int uiTimeout)
{
    int iSelectRes;
    unsigned int uiIdx;
    unsigned long ulRet = 0;
    socketh_t iSockH;
    socketh_t iMaxH = 0;
    fd_set fdsetRead;
    fd_set fdsetWrite;
    fd_set fdsetError;
    struct timeval sTimeout;

    if (ipSockH == NULL)
        throw clXSocket("clSocketOp::WaitReadWrite(): ipSockH == NULL");
    FD_ZERO(&fdsetRead);
    FD_ZERO(&fdsetWrite);
    FD_ZERO(&fdsetError);
    for (uiIdx = 0; uiIdx < uiHandleCount; uiIdx++)
    {
        iSockH = ipSockH[uiIdx];
        FD_SET(iSockH, &fdsetRead);
        FD_SET(iSockH, &fdsetWrite);
        FD_SET(iSockH, &fdsetError);
        if (iSockH > iMaxH)
            iMaxH = iSockH;
    }
    sTimeout.tv_sec = uiTimeout / 1000u;
    sTimeout.tv_usec = uiTimeout % 1000u * 1000u;
    iSelectRes = select(iMaxH + 1, &fdsetRead, &fdsetWrite, &fdsetError,
        &sTimeout);
    if (iSelectRes == 0)
        return 0ul;
    else if (iSelectRes < 0)
        throw clXSocket("clSocket::WaitReadWrite(): select()", errno);
    for (uiIdx = 0; uiIdx < uiHandleCount; uiIdx++)
    {
        iSockH = ipSockH[uiIdx];
        if (FD_ISSET(iSockH, &fdsetError))
            throw clXSocketError("clSocket::WaitReadWrite(): socket error");
        if (FD_ISSET(iSockH, &fdsetRead) || FD_ISSET(iSockH, &fdsetWrite))
            ulRet |= (1 << uiIdx);
    }
    return ulRet;
}


ssize_t clSocket::Read (void *vpBuffer, size_t sizeBuf)
{
    ssize_t iReadRes;

    iReadRes = read(iHandle, vpBuffer, sizeBuf);
    if (iReadRes < 0)
        throw clXSocket("clSocket::Read(): read()", errno);
    return iReadRes;
}


ssize_t clSocket::Write (const void *vpBuffer, size_t sizeBuf)
{
    ssize_t iWriteRes;

    iWriteRes = write(iHandle, vpBuffer, sizeBuf);
    if (iWriteRes < 0)
        throw clXSocket("clSocket::Write(): write()", errno);
    return iWriteRes;
}


ssize_t clSocket::Recv (void *vpBuffer, size_t sizeBuf, int iFlags)
{
    ssize_t iRecvRes;

    iRecvRes = recv(iHandle, vpBuffer, sizeBuf, iFlags);
    if (iRecvRes < 0)
        throw clXSocket("clSocket::Recv(): recv()", errno);
    return iRecvRes;
}


ssize_t clSocket::Send (const void *vpBuffer, size_t sizeBuf, int iFlags)
{
    ssize_t iSendRes;

    iSendRes = send(iHandle, vpBuffer, sizeBuf, iFlags);
    if (iSendRes < 0)
        throw clXSocket("clSocket::Send(): send()", errno);
    return iSendRes;
}


ssize_t clSocket::RecvLine (char *cpBuffer, size_t sizeBuf)
{
    ssize_t iRecvRes;
    clString RecvString;
    
    iRecvRes = RecvLine(RecvString);
    strncpy(cpBuffer, RecvString.c_str(), sizeBuf);
    cpBuffer[sizeBuf - 1] = '\0';
    return ((iRecvRes < (ssize_t) sizeBuf) ? iRecvRes : sizeBuf);
}


ssize_t clSocket::RecvLine (clString &RecvString)
{
    ssize_t iRecvRes;
    char *cpRecvBuf;
    
    if (RecvBuf.Size() != SOCKETOP_READBUF_SIZE)
        RecvBuf.SetSize(SOCKETOP_READBUF_SIZE);
    cpRecvBuf = RecvBuf;
    while (true)
    {
        iRecvRes = Recv(cpRecvBuf, SOCKETOP_READBUF_SIZE, MSG_PEEK);
        ulRecvCount = iRecvRes;
        if (GetLine(RecvString))
        {
            iRecvRes = Recv(cpRecvBuf, RecvString.Length(), 0);
            break;
        }
        else if (ulRecvCount >= SOCKETOP_READBUF_SIZE)
            throw clXSocket(
                "clSocket::RecvLine(): line receive buffer overflow");
    }
    return RecvString.Length();
}


ssize_t clSocket::SendLine (const char *cpBuffer)
{
    ssize_t iSendRes;
    clString Line(cpBuffer);

    if (Line[Line.Length() - 1] != '\n')
        Line += "\n";
    iSendRes = Send(Line.c_str(), Line.Length(), 0);
    if (iSendRes != (ssize_t) Line.Length())
        throw clXSocket(
            "clSocket::SendLine(): short Send()");
    return iSendRes;
}


ssize_t clSocket::SendLine (const clString &SendString)
{
    return SendLine(SendString.c_str());
}


ssize_t clSocket::RecvFrom (void *vpBuffer, size_t sizeBuf, int iFlags,
    sockaddrptr_t spAddr, socklen_t *sizepAddr)
{
    ssize_t iRecvRes;

    iRecvRes = recvfrom(iHandle, vpBuffer, sizeBuf, iFlags,
        spAddr, sizepAddr);
    if (iRecvRes < 0)
        throw clXSocket("clSocket::RecvFrom(): recvfrom()", errno);
    return iRecvRes;
}


ssize_t clSocket::SendTo (const void *vpBuffer, size_t sizeBuf, int iFlags,
    const sockaddrptr_t spAddr, socklen_t sizeAddr)
{
    ssize_t iSendRes;

    iSendRes = sendto(iHandle, vpBuffer, sizeBuf, iFlags,
        spAddr, sizeAddr);
    if (iSendRes < 0)
        throw clXSocket("clSocket::SendTo(): sendto()", errno);
    return iSendRes;
}


ssize_t clSocket::RecvMsg (struct msghdr *spMsg, int iFlags)
{
    ssize_t iRecvRes;

    iRecvRes = recvmsg(iHandle, spMsg, iFlags);
    if (iRecvRes < 0)
        throw clXSocket("clSocket::RecvMsg(): recvmsg()", errno);
    return iRecvRes;
}


ssize_t clSocket::SendMsg (const struct msghdr *spMsg, int iFlags)
{
    ssize_t iSendRes;

    iSendRes = sendmsg(iHandle, spMsg, iFlags);
    if (iSendRes < 0)
        throw clXSocket("clSocket::SendMsg(): sendmsg()", errno);
    return iSendRes;
}


void clSocket::Close ()
{
    if (SOCKET_HANDLE_OK(iHandle))
    {
        if (close(iHandle) != 0)
            throw clXSocket("clSocket::Close(): close()", errno);
    }
}


void clSocket::SetBlocking (bool bBlocking)
{
    int iFlags;

    iFlags = fcntl(iHandle, F_GETFL);
    if (iFlags < 0)
        throw clXSocket("clSocket::SetBlocking(): fcntl(F_GETFL)");
    if (bBlocking)
        iFlags |= O_NONBLOCK;
    else
        iFlags &= ~O_NONBLOCK;
    if (fcntl(iHandle, F_SETFL, iFlags) < 0)
        throw clXSocket("clSocket::SetBlocking(): fcntl(F_SETFL)");
}


void clSocket::SetOption (int iLevel, int iOption, const void *vpValue,
    socklen_t iValueSize)
{
    if (setsockopt(iHandle, iLevel, iOption, vpValue, iValueSize) < 0)
        throw clXSocket("clSocket::SetOption(): setsockopt()", errno);
}


void clSocket::GetOption (int iLevel, int iOption, void *vpValue,
    socklen_t *ipValueSize)
{
    if (getsockopt(iHandle, iLevel, iOption, vpValue, ipValueSize) < 0)
        throw clXSocket("clSocket::GetOption(): getsockopt()", errno);
}


void clSocket::SetNagle (bool bNagle)
{
    int iNagle;

    iNagle = (bNagle) ? 0 : 1;
    try
    {
        SetOption(IPPROTO_TCP, TCP_NODELAY, &iNagle, sizeof(iNagle));
    }
    catch (clXSocket &XS)
    {
        throw clXSocket(
            clString("clSocket::SetNagle(): ") +
            XS.what());
    }
}


void clSocket::GetSockName (clSocketAddr &Addr)
{
    socklen_t iAddrSize;

    iAddrSize = Addr.Size();
    if (getsockname(iHandle, Addr.Get(), &iAddrSize) < 0)
        throw clXSocket("clSocket::GetSockName(): getsockname()", errno);
}


clSocketAddr * clSocket::GetSockName ()
{
    // we try to reserve enough space for any possible address
    unsigned char cpAddrBuf[SOCKETOP_MAX_ADDR_SIZE];
    socklen_t uiAddrSize;

    uiAddrSize = SOCKETOP_MAX_ADDR_SIZE;
    memset(cpAddrBuf, 0x00, sizeof(cpAddrBuf));
    if (getsockname(iHandle, (struct sockaddr *) cpAddrBuf, &uiAddrSize) < 0)
        throw clXSocket("clSocket::GetSockName(): getsockname()", errno);
    return clSocketAddrDetect::Get((sockaddrptr_t) cpAddrBuf, uiAddrSize);
}


void clSocket::GetPeerName (clSocketAddr &Addr)
{
    socklen_t iAddrSize;

    iAddrSize = Addr.Size();
    if (getpeername(iHandle, Addr.Get(), &iAddrSize) < 0)
        throw clXSocket("clSocket::GetPeerName(): getpeername()", errno);
}


clSocketAddr * clSocket::GetPeerName ()
{
    // we try to reserve enough space for any possible address
    unsigned char cpAddrBuf[SOCKETOP_MAX_ADDR_SIZE];
    socklen_t uiAddrSize;

    uiAddrSize = SOCKETOP_MAX_ADDR_SIZE;
    memset(cpAddrBuf, 0x00, sizeof(cpAddrBuf));
    if (getpeername(iHandle, (struct sockaddr *) cpAddrBuf, &uiAddrSize) < 0)
        throw clXSocket("clSocket::GetPeerName(): getpeername()", errno);
    return clSocketAddrDetect::Get((sockaddrptr_t) cpAddrBuf, uiAddrSize);
}


bool clSocket::Connected ()
{
    unsigned char ucByte;
    int iError;
    socklen_t iOptSize;

    if (SOCKET_INVALID_HANDLE(iHandle))
        return false;
    try
    {
        iOptSize = sizeof(iError);
        GetOption(SOL_SOCKET, SO_ERROR, &iError, &iOptSize);
        if (iError != 0) return false;
        if (WaitRead(0))
        {
            if (Recv(&ucByte, sizeof(ucByte), MSG_PEEK) <= 0)
                return false;
        }
    }
    catch (clXSocket &XS)
    {
        return false;
    }
    return true;
}


void clSocket::McastJoin (const clSocketAddr &Addr, const clSocketAddr *pIf)
{
    if (clString(Addr.GetType()) == "ipv4")
    {
        clSocketAddrIPv4 *pAddr4 = (clSocketAddrIPv4 *) &Addr;
        clSocketAddrIPv4 *pIf4 = (clSocketAddrIPv4 *) pIf;
        clSocketAddrIPv4 McastReqAddr;
    
        SetOption(IPPROTO_IP, IP_ADD_MEMBERSHIP,
            McastReqAddr.GetMcastReq(*pAddr4, pIf4),
            McastReqAddr.GetMcastReqSize());
    }
    else if (clString(Addr.GetType()) == "ipv6")
    {
        clSocketAddrIPv6 *pAddr6 = (clSocketAddrIPv6 *) &Addr;
        clSocketAddrIPv6 McastReqAddr;

        SetOption(IPPROTO_IPV6, IP_ADD_MEMBERSHIP,
            McastReqAddr.GetMcastReq(*pAddr6, (unsigned int) pIf),
            McastReqAddr.GetMcastReqSize());
    }
}


void clSocket::McastLeave (const clSocketAddr &Addr, const clSocketAddr *pIf)
{
    if (clString(Addr.GetType()) == "ipv4")
    {
        clSocketAddrIPv4 *pAddr4 = (clSocketAddrIPv4 *) &Addr;
        clSocketAddrIPv4 *pIf4 = (clSocketAddrIPv4 *) pIf;
        clSocketAddrIPv4 McastReqAddr;
    
        SetOption(IPPROTO_IP, IP_DROP_MEMBERSHIP,
            McastReqAddr.GetMcastReq(*pAddr4, pIf4),
            McastReqAddr.GetMcastReqSize());
    }
    else if (clString(Addr.GetType()) == "ipv6")
    {
        clSocketAddrIPv6 *pAddr6 = (clSocketAddrIPv6 *) &Addr;
        clSocketAddrIPv6 McastReqAddr;

        SetOption(IPPROTO_IPV6, IP_DROP_MEMBERSHIP,
            McastReqAddr.GetMcastReq(*pAddr6, (unsigned int) pIf),
            McastReqAddr.GetMcastReqSize());
    }
}


void clSocket::McastSetInterface (const clSocketAddr &If)
{
    if (clString(If.GetType()) == "ipv4")
    {
        struct sockaddr_in *spAddr;

        spAddr = (struct sockaddr_in *) If.Get();
        SetOption(IPPROTO_IP, IP_MULTICAST_IF, &spAddr->sin_addr,
            sizeof(struct in_addr));
    }
    else if (clString(If.GetType()) == "ipv6")
    {
        struct sockaddr_in6 *spAddr;

        spAddr = (struct sockaddr_in6 *) If.Get();
        SetOption(IPPROTO_IPV6, IP_MULTICAST_IF, &spAddr->sin6_addr,
            sizeof(struct in6_addr));
    }
}


void clSocket::McastSetTTL (const clSocketAddr &Addr, unsigned char ucTTL)
{
    if (clString(Addr.GetType()) == "ipv4")
    {
        SetOption(IPPROTO_IP, IP_MULTICAST_TTL,
            &ucTTL, sizeof(ucTTL));
    }
    else if (clString(Addr.GetType()) == "ipv6")
    {
        SetOption(IPPROTO_IPV6, IP_MULTICAST_TTL,
            &ucTTL, sizeof(ucTTL));
    }
}


void clSocket::McastSetLoop (const clSocketAddr &Addr, bool bLoop)
{
    unsigned char ucLoop;

    ucLoop = (bLoop) ? 1 : 0;
    if (clString(Addr.GetType()) == "ipv4")
    {
        SetOption(IPPROTO_IP, IP_MULTICAST_LOOP,
            &ucLoop, sizeof(ucLoop));
    }
    else if (clString(Addr.GetType()) == "ipv6")
    {
        SetOption(IPPROTO_IPV6, IP_MULTICAST_LOOP,
            &ucLoop, sizeof(ucLoop));
    }
}


clSocket & clSocket::operator = (const clSocket &Src)
{
    bCloseOnDestruct = Src.bCloseOnDestruct;
    iHandle = dup(Src.GetHandle());
    return (*this);
}
