// $Id: SocketTypes.hh,v 1.1 2005/12/17 15:50:11 sonarnerd Exp $

/*

    Socket base definitions
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


#ifndef SOCKETTYPES_HH
#define SOCKETTYPES_HH

#define SOCKET_INVALID_SOCKET           -1
#define SOCKET_HANDLE_OK(x)             (x >= 0)
#define SOCKET_INVALID_HANDLE(x)        (x < 0)


typedef int socketh_t;


#endif
