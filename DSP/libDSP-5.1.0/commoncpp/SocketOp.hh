// $Id: SocketOp.hh,v 1.1 2005/12/17 15:50:11 sonarnerd Exp $

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


#ifndef SOCKETOP_HH
#define SOCKETOP_HH

#include <cstddef>
#include <sys/socket.h>

#include <commoncpp/Exception.hh>
#include <commoncpp/Compilers.hh>
#include <commoncpp/Alloc.hh>
#include <commoncpp/String.hh>
#include <commoncpp/SocketTypes.hh>
#include <commoncpp/SocketAddr.hh>

#define SOCKETOP_MAX_ADDR_SIZE          4096
#define SOCKETOP_READBUF_SIZE           4096


    /**
        Exception class for clSocket.
    */
    class clXSocket : public clException
    {
        public:
            clXSocket (const char *cpErrorMsg) :
                clException(cpErrorMsg)
                { }
            clXSocket (const char *cpErrorMsg, int iErrorCode) :
                clException(cpErrorMsg)
                { }
    };


    /**
        Exception class for handling errors from Wait*() methods.
    */
    class clXSocketError : public clXSocket
    {
        public:
            clXSocketError (const char *cpErrorMsg) :
                clXSocket(cpErrorMsg)
                { }
            clXSocketError (const char *cpErrorMsg, int iErrorCode) :
                clXSocket(cpErrorMsg, iErrorCode)
                { }
    };


    /**
        Class of socket operations.
        
        Also a container object for single socket.
    */
    class clSocket
    {
        private:

            unsigned long ulRecvCount;  ///< Receive count for line based I/O
            clAlloc RecvBuf;  ///< Receive buffer for line based I/O

            bool GetLine (clString &);

        protected:

            bool bCloseOnDestruct;  ///< Close on destruct
            socketh_t iHandle;  ///< Socket handle

        public:

            /**
                Constructor.
            */
            clSocket ()
                {
                    ulRecvCount = 0;
                    bCloseOnDestruct = true;
                    iHandle = SOCKET_INVALID_SOCKET;
                }
            /**
                Copy constructor.
            */
            clSocket (const clSocket &CopySrc)
                { this->operator=(CopySrc); }
            /// \overload
            clSocket (socketh_t iNewH)
                {
                    bCloseOnDestruct = true;
                    this->operator=(iNewH);
                }
            /**
                Destructor.
            */
            virtual ~clSocket ();
            /**
                Create socket based on address family (eg. AF_UNIX, AF_INET or
                AF_INET6), type (eg. SOCK_STREAM, SOCK_DGRAM or SOCK_SEQPACKET)
                and protocol (eg. IPPROTO_TCP or IPPROTO_UDP).
                
                \note Protocol 0 means default protocol for adress family/type.
                
                \param iAddrFamily Address family
                \param iType Socket type
                \param iProtocol Protocol
                \throw clXSocket
            */
            void Create (int, int, int);
            /**
                Create socket pair of connected local (UNIX) domain sockets.
                
                Returns one instance, another instance is "this".
                
                \param iType Socket type
                \return Handle to the other socket
                \throw clXSocket
            */
            socketh_t CreatePair (int);
            /**
                Bind socket to a specified address.
                
                \param Addr Address
                \throw clXSocket
            */
            void Bind (const clSocketAddr &);
            /**
                Connect to a specified address.
                
                \param Addr Address
                \throw clXSocket
            */
            void Connect (const clSocketAddr &);
            /**
                Connect to a specified address.
                
                \param Addr Address
                \param uiTimeout Connect timeout (ms)
                \return Connection succeeded?
                \throw clXSocket
            */
            bool Connect (const clSocketAddr &, unsigned int);
            /**
                Start listening for connections.
                
                \param iBacklog Number of connections to keep in backlog
                \throw clXSocket
            */
            void Listen (int = SOMAXCONN);
            /**
                Accept incoming connection.
                
                \param spAddr Generic address pointer
                \param uipAddrSize Size of generic address struct
                \return Handle to a new connection
                \throw clXSocket
            */
            socketh_t Accept (sockaddrptr_t = NULL, socklen_t * = NULL);
            /**
                \overload
                
                \param Addr Address class instance
                \return Handle to a new connection
                \throw clXSocket
            */
            socketh_t Accept (clSocketAddr &);
            /**
                \overload
                
                \param ppAddr Source address of the incoming connection
                \return Handle to a new connection
                \throw clXSocket
            */
            socketh_t Accept (clSocketAddr **);
            /**
                Wait for socket to become readable.
                
                \param uiTimeout Timeout (ms)
                \return Readable?
                \throw clXSocket
                \throw clXSocketError
            */
            bool WaitRead (unsigned int);
            /**
                Wait for socket to become writable.
                
                \param uiTimeout Timeout (ms)
                \return Writable?
                \throw clXSocket
                \throw clXSocketError
            */
            bool WaitWrite (unsigned int);
            /**
                Wait for socket to become readable or writable.
                
                \param uiTimeout Timeout (ms)
                \return Readable or writable?
                \throw clXSocket
                \throw clXSocketError
            */
            bool WaitReadWrite (unsigned int);
            /**
                Wait for some of the specified sockets to become readable.
                
                \note Each bit in return value represent single handle.
                
                \param ipSockH Socket handles
                \param uiHandleCount Number of handles
                \param uiTimeout Timeout (ms)
                \return Bit mask of readable handles
                \throw clXSocket
                \throw clXSocketError
            */
            static unsigned long WaitRead (const socketh_t *, unsigned int,
                unsigned int);
            /**
                Wait for some of the specified sockets to become writable.
                
                \note Each bit in return value represent single handle.
                
                \param ipSockH Socket handles
                \param uiHandleCount Number of handles
                \param uiTimeout Timeout (ms)
                \return Bit mask of writable handles
                \throw clXSocket
                \throw clXSocketError
            */
            static unsigned long WaitWrite (const socketh_t *, unsigned int,
                unsigned int);
            /**
                Wait for some of the specified sockets to become readable or writable.
                
                \note Each bit in return value represent single handle.
                
                \param ipSockH Socket handles
                \param uiHandleCount Number of handles
                \param uiTimeout Timeout (ms)
                \return Bit mask of readable or writable handles
                \throw clXSocket
                \throw clXSocketError
            */
            static unsigned long WaitReadWrite (const socketh_t *, unsigned int,
                unsigned int);
            /**
                Read data from socket.
                
                \param vpBuffer Pointer to buffer
                \param sizeBuf Size of buffer
                \return Number of bytes read
                \throw clXSocket
            */
            ssize_t Read (void *, size_t);
            /**
                Write data to socket.
                
                \param vpBuffer Pointer to buffer
                \param sizeBuf Size of data in buffer
                \return Number of bytes written
                \throw clXSocket
            */
            ssize_t Write (const void *, size_t);
            /**
                Receive data from socket.
                
                \note Flags can be OR'ed set of MSG_PEEK, MSG_OOB, MSG_WAITALL.
                
                \param vpBuffer Pointer to buffer
                \param sizeBuf Size of buffer
                \param iFlags Flags
                \return Number of bytes received
                \throw clXSocket
            */
            ssize_t Recv (void *, size_t, int);
            /**
                Send data to socket.
                
                \note Flags can be OR'ed set of MSG_EOR, MSG_OOB.
                
                \param vpBuffer Pointer to buffer
                \param sizeBuf Size of data in buffer
                \param iFlags Flags
                \return Number of bytes sent
                \throw clXSocket
            */
            ssize_t Send (const void *, size_t, int);
            /**
                Receive line of text.
                
                \note Terminating new line is removed.
                
                \param cpBuffer Pointer to buffer
                \param sizeBuf Size of buffer
                \return Number of characters received
                \throw clXSocket
            */
            ssize_t RecvLine (char *, size_t);
            /**
                \overload
                
                \param RecvString Reference to the destination string
                \return Number of characters received
                \throw clXSocket
            */
            ssize_t RecvLine (clString &);
            /**
                Send line of text.
                
                \note Terminating new line is appended if missing.
                
                \param cpBuffer Pointer to buffer
                \return Number of characters sent
                \throw clXSocket
            */
            ssize_t SendLine (const char *);
            /**
                \overload
                
                \param SendString Reference to the source string
                \return Number of characters sent
                \throw clXSocket
            */
            ssize_t SendLine (const clString &);
            /**
                Receive data with address.
                
                \note Flags can be OR'ed set of MSG_PEEK, MSG_OOB, MSG_WAITALL.
                
                \param vpBuffer Pointer to buffer
                \param sizeBuf Size of buffer
                \param iFlags Flags
                \param spAddr Pointer to address struct
                \param sizepAddr Pointer to size of address struct, size of address on return
                \throw clXSocket
            */
            ssize_t RecvFrom (void *, size_t, int,
                sockaddrptr_t, socklen_t *);
            /**
                Send data to address.
                
                \note Flags can be OR'ed set of MSG_EOR, MSG_OOB.
                
                \param vpBuffer Pointer to buffer
                \param iBufSize Size of data in buffer
                \param iFlags Flags
                \param spAddr Pointer to address struct
                \param sizeAddr Size of address struct
                \throw clXSocket
            */
            ssize_t SendTo (const void *, size_t, int,
                const sockaddrptr_t, socklen_t);
            /**
                Receive message.
                
                \note Flags can be OR'ed set of MSG_OOB, MSG_PEEK, MSG_WAITALL.
                
                \param spMsg Message
                \param iFlags Flags
                \throw clXSocket
            */
            ssize_t RecvMsg (struct msghdr *, int);
            /**
                Send message.
                
                \note Flags can be OR'ed set of MSG_EOR, MSG_OOB.
                
                \param spMsg Message
                \param iFlags Flags
                \throw clXSocket
            */
            ssize_t SendMsg (const struct msghdr *, int);
            /**
                Close socket.
                
                \throw clXSocket
            */
            void Close ();
            /**
                Get reference to the socket handle.
                
                \return Reference to socket handle.
            */
            socketh_t & Handle ()
                { return iHandle; }
            /**
                Get socket handle.
                
                \return Socket handle
            */
            socketh_t GetHandle () const
                { return iHandle; }
            /**
                Set close on destruct behaviour.
                
                \param bClose Close on destruct?
            */
            void SetCloseOnDestruct (bool bClose)
                { bCloseOnDestruct = bClose; }
            /**
                Set socket to blocking/non-blocking mode.
                
                \param bBlocking Blocking mode?
                \throw clXSocket
            */
            void SetBlocking (bool);
            /**
                Set socket option.
                
                \param iLevel Level
                \param iOption Option
                \param vpValue Option value
                \param iValueSize Size of option value
                \throw clXSocket
            */
            void SetOption (int, int, const void *, socklen_t);
            /**
                Get socket option.
                
                \param iLevel Level
                \param iOption Option
                \param vpValue Option value
                \param ipValueSize Size of option value
                \throw clXSocket
            */
            void GetOption (int, int, void *, socklen_t *);
            /**
                Enable/disable TCP Nagle-algorithm.
                
                \param bNagle Nagle state
                \throw clXSocket
            */
            void SetNagle (bool);
            /**
                Get socket name.
                
                \param Addr Address class instance
                \throw clXSocket
            */
            void GetSockName (clSocketAddr &);
            /**
                \overload
                
                \return Address class instance
                \throw clXSocket
            */
            clSocketAddr * GetSockName ();
            /**
                Get peer name.
                
                \param Addr Address class instance
                \throw clXSocket
            */
            void GetPeerName (clSocketAddr &);
            /**
                \overload
                
                \return Address class instance
                \throw clXSocket
            */
            clSocketAddr * GetPeerName ();
            /**
                Is socket connected without errors?
                
                \return Connection OK?
            */
            bool Connected ();

            /**
                Join multicast group.
                
                \note NULL as interface address means any address.
                
                \param Addr Group address
                \param pIf Interface address
                \throw clXSocket
            */
            void McastJoin (const clSocketAddr &, const clSocketAddr * = NULL);
            /**
                Leave multicast group.
                
                \note NULL as interface address means any address.
                
                \param Addr Group address
                \param pIf Interface address
                \throw clXSocket
            */
            void McastLeave (const clSocketAddr &, const clSocketAddr * = NULL);
            /**
                Set multicast interface.
                
                \param If Interface address
                \throw clXSocket
            */
            void McastSetInterface (const clSocketAddr &);
            /**
                Set multicast TTL.
                
                \param Addr Group address
                \param ucTLL Max hop count
                \throw clXSocket
            */
            void McastSetTTL (const clSocketAddr &, unsigned char);
            /**
                Set multicast loopback.
                
                \param Addr Group address
                \param bLoop Multicast loopback?
                \throw clXSocket
            */
            void McastSetLoop (const clSocketAddr &, bool);
                

            /**
                Assign socket.
                
                \note Handle is copied!
            */
            clSocket & operator = (const clSocket &);
            /**
                Assign socket handle.
                
                \note Handle is not copied!
            */
            clSocket & operator = (socketh_t iNewH)
                {
                    iHandle = iNewH;
                    return (*this);
                }
    };


#endif
