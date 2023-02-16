// $Id: SocketAddr.hh,v 1.1 2005/12/17 15:50:11 sonarnerd Exp $

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


#ifndef SOCKETADDR_HH
#define SOCKETADDR_HH


#include <cstddef>
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>

#include <commoncpp/Exception.hh>
#include <commoncpp/Compilers.hh>
#include <commoncpp/String.hh>
#include <commoncpp/SocketTypes.hh>

#define SOCKETADDR_MAX_HOSTNAME_LENGTH  255

#ifndef UNIX_PATH_MAX
#define UNIX_PATH_MAX                   108
#endif
#if (UNIX_PATH_MAX != 108)
#error Check the size of sockaddr_un.sun_path!
#endif


typedef struct sockaddr * sockaddrptr_t;


    /**
        Exception class for clSocketAddr.
    */
    class clXSocketAddr : public clException
    {
        public:
            clXSocketAddr (const char *cpErrorMsg) :
                clException(cpErrorMsg)
                { }
            clXSocketAddr (const char *cpErrorMsg, int iErrorCode) :
                clException(cpErrorMsg, iErrorCode)
                { }
    };


    /**
        Abstract base class for different address families.
    */
    class clSocketAddr
    {
        public:

            /**
                Constructor.
            */
            clSocketAddr ()
                { }
            /**
                Destructor.
            */
            virtual ~clSocketAddr ()
                { }
            /**
                Get address type string.
                
                \return Address type ("null")
            */
            virtual const char * GetType () const
                { return "null"; }
            /**
                Method for getting generic address pointer.

                \return Generic address pointer
            */
            virtual sockaddrptr_t Get () const
                { return NULL; }
            /**
                Method for getting size of generic address.

                \return Size of address
            */
            virtual socklen_t Size () const
                { return 0; }
            /**
                Set address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
                \throw clXSocketAddr
            */
            virtual void Set (sockaddrptr_t, socklen_t)
                { }
            /**
                Get address as string, without port number.
                
                \return Address as string w/o port
                \throw clXSocketAddr
            */
            virtual clString AsString () const
                { return clString(); }
            /**
                Get address as string, with port number.
                
                \return Address as string w/ port
                \throw clXSocketAddr
            */
            virtual clString AsStringP () const
                { return AsString(); }
            /**
                Get address as hostname, without port number.
                
                \return Address as hostname w/o port
                \throw clXSocketAddr
            */
            virtual clString AsHostName () const
                { return clString(); }
            /**
                Get address as hostname, with port number.
                
                \return Address as hostname w/ port
                \throw clXSocketAddr
            */
            clString AsHostNameP () const
                { return AsHostName(); }

            virtual clSocketAddr & operator = (const clSocketAddr &Src)
                { return (*this); }
    };


    /**
        Class for local (UNIX) address family.
    */
    class clSocketAddrLocal : public clSocketAddr
    {
        private:

            void Init ();

        protected:

            struct sockaddr_un sAddr;

        public:

            /**
                Constructor.
            */
            clSocketAddrLocal ()
                { Init(); }
            /**
                Copy constructor.
            */
            clSocketAddrLocal (const clSocketAddrLocal &CopySrc)
                { this->operator=(CopySrc); }
            /**
                Constructor from address string.
                
                \param cpFileName Socket filename
            */
            clSocketAddrLocal (const char *cpFileName)
                { Set(cpFileName); }
            /**
                Constructor which sets the address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
            */
            clSocketAddrLocal (sockaddrptr_t spNewAddr, socklen_t uiAddrSize)
                { Set(spNewAddr, uiAddrSize); }
            /**
                Destructor.
            */
            virtual ~clSocketAddrLocal ()
                { }
            /**
                Get address type string.
                
                \return Address type ("local")
            */
            virtual const char * GetType () const
                { return "local"; }
            /**
                Method to get generic address pointer for local address.
                
                \return Generic address pointer
            */
            virtual sockaddrptr_t Get () const
                { return ((sockaddrptr_t) &sAddr); }
            /**
                Method for getting size of generic address.

                \return Size of address
            */
            virtual socklen_t Size () const
                { return (sizeof(sAddr)); }
            /**
                Method for setting address from string.
                
                \param cpFileName Filename for the socket
                \throw clXSocketAddr
            */
            void Set (const char *);
            /**
                Set address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
                \throw clXSocketAddr
            */
            virtual void Set (sockaddrptr_t, socklen_t);
            /**
                Get address as string.
                
                \return Address as string
                \throw clXSocketAddr
            */
            virtual clString AsString () const;

            virtual clSocketAddrLocal & operator = (const clSocketAddrLocal &);
    };


    /**
        Class for IPv4 address family.
    */
    class clSocketAddrIPv4 : public clSocketAddr
    {
        private:

            void Init ();

        protected:

            struct sockaddr_in sAddr;  ///< IPv4 address
            struct ip_mreq sMreq;  ///< IPv4 multicast request

        public:

            /**
                Constructor.
            */
            clSocketAddrIPv4 ()
                { Init(); }
            /**
                Copy constructor.
            */
            clSocketAddrIPv4 (const clSocketAddrIPv4 &CopySrc)
                { this->operator=(CopySrc); }
            /**
                Constructor which sets the address from address string.
                
                \param cpAddress Address string
                \param uiPort Port number
                \throw clXSocketAddr
            */
            clSocketAddrIPv4 (const char *cpAddress, unsigned short uiPort)
                { Set(cpAddress, uiPort); }
            /**
                Constructor which sets the address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
            */
            clSocketAddrIPv4 (sockaddrptr_t spNewAddr, socklen_t uiAddrSize)
                { Set(spNewAddr, uiAddrSize); }
            /**
                Destructor.
            */
            virtual ~clSocketAddrIPv4 ()
                { }
            /**
                Get address type string.
                
                \return Address type ("ipv4")
            */
            virtual const char * GetType () const
                { return "ipv4"; }
            /**
                Method to get generic address pointer for IPv4 address.
                
                \return Generic address pointer
            */
            virtual sockaddrptr_t Get () const
                { return ((sockaddrptr_t) &sAddr); }
            /**
                Method for getting size of generic address.

                \return Size of address
            */
            virtual socklen_t Size () const
                { return (sizeof(sAddr)); }
            /**
                Method for setting address from string.
                
                \note Set Address to NULL and port to 0 for any address.
                
                \param cpAddress Address string
                \param uiPort Port number
                \throw clXSocketAddr
            */
            void Set (const char * = NULL, unsigned short = 0);
            /**
                Set address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
                \throw clXSocketAddr
            */
            virtual void Set (sockaddrptr_t, socklen_t);
            /**
                Get address as string, without port number.
                
                \return Address as string w/o port
                \throw clXSocketAddr
            */
            virtual clString AsString () const;
            /**
                Get address as string, with port number.
                
                \return Address as string w/ port
                \throw clXSocketAddr
            */
            virtual clString AsStringP () const;
            /**
                Get address as hostname, without port number.
                
                \return Address as hostname w/o port
                \throw clXSocketAddr
            */
            virtual clString AsHostName () const;
            /**
                Get address as hostname, with port number.
                
                \return Address as hostname w/ port
                \throw clXSocketAddr
            */
            virtual clString AsHostNameP () const;
            /**
                Get multicast request.
                
                \param Addr Multicast address
                \param pIf Multicast interface (NULL = any address)
                \throw clXSocketAddr
            */
            struct ip_mreq * GetMcastReq (const clSocketAddrIPv4 &,
                const clSocketAddrIPv4 * = NULL);
            /**
                Get size of multicast request;
                
                \return Size of multicast request
            */
            socklen_t GetMcastReqSize ()
                { return (sizeof(sMreq)); }

            virtual clSocketAddrIPv4 & operator = (const clSocketAddrIPv4 &);
    };


    /**
        Class for IPv6 address family.
    */
    class clSocketAddrIPv6 : public clSocketAddr
    {
        private:

            void Init ();

        protected:

            struct sockaddr_in6 sAddr;  ///< IPv6 address
            struct ipv6_mreq sMreq;  ///< IPv6 multicast request

        public:

            /**
                Constructor.
            */
            clSocketAddrIPv6 ()
                { Init(); }
            /**
                Copy constructor.
            */
            clSocketAddrIPv6 (const clSocketAddrIPv6 &CopySrc)
                { this->operator=(CopySrc); }
            /**
                Constructor which sets the address from address string.
                
                \param cpAddress Address string
                \param uiPort Port number
                \throw clXSocketAddr
            */
            clSocketAddrIPv6 (const char *cpAddress, unsigned short uiPort)
                { Set(cpAddress, uiPort); }
            /**
                Constructor which sets the address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
            */
            clSocketAddrIPv6 (sockaddrptr_t spNewAddr, socklen_t uiAddrSize)
                { Set(spNewAddr, uiAddrSize); }
            /**
                Destructor.
            */
            virtual ~clSocketAddrIPv6 ()
                { }
            /**
                Get address type string.
                
                \return Address type ("ipv6")
            */
            virtual const char * GetType () const
                { return "ipv6"; }
            /**
                Method to get generic address pointer for IPv4 address.
                
                \return Generic address pointer
            */
            virtual sockaddrptr_t Get () const
                { return ((sockaddrptr_t) &sAddr); }
            /**
                Method for getting size of generic address.

                \return Size of address
            */
            virtual socklen_t Size () const
                { return (sizeof(sAddr)); }
            /**
                Method for setting address from string.
                
                \note Set Address to NULL and port to 0 for any address.
                
                \param cpAddress Address string
                \param uiPort Port number
                \throw clXSocketAddr
            */
            void Set (const char * = NULL, unsigned short = 0);
            /**
                Set address from generic address pointer.
                
                \param spNewAddr Generic address pointer
                \param uiAddrSize Size of address
                \throw clXSocketAddr
            */
            virtual void Set (sockaddrptr_t, socklen_t);
            /**
                Get address as string, without port number.
                
                \return Address as string w/o port
                \throw clXSocketAddr
            */
            virtual clString AsString () const;
            /**
                Get address as string, with port number.
                
                \return Address as string w/ port
                \throw clXSocketAddr
            */
            virtual clString AsStringP () const;
            /**
                Get address as hostname, without port number.
                
                \return Address as hostname w/o port
                \throw clXSocketAddr
            */
            virtual clString AsHostName () const;
            /**
                Get address as hostname, with port number.
                
                \return Address as hostname w/ port
                \throw clXSocketAddr
            */
            virtual clString AsHostNameP () const;
            /**
                Get multicast request.
                
                \param Addr Multicast address
                \param uiIf Multicast interface
                \throw clXSocketAddr
            */
            struct ipv6_mreq * GetMcastReq (const clSocketAddrIPv6 &,
                unsigned int = 0);
            /**
                Get size of multicast request;
                
                \return Size of multicast request
            */
            socklen_t GetMcastReqSize ()
                { return (sizeof(sMreq)); }

            virtual clSocketAddrIPv6 & operator = (const clSocketAddrIPv6 &);
    };


    class clSocketAddrDetect
    {
        public:
            /**
                Get address class from generic address pointer.
                
                \param spAddr Generic address pointer
                \param uiAddrSize Size of address
                \throw clXSocketAddr
            */
            static clSocketAddr * Get (sockaddrptr_t, socklen_t);
    };


#endif
