/*

    DSP types
    Copyright (C) 2002 Jussi Laako

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#ifndef DSPTYPES_H
    #define DSPTYPES_H

    #pragma pack(1)

    /**
        Single precision complex datatype. 
        Type prefix is left out to prevent need of source code changes when 
        precision is changed.
    */
    typedef struct _sSCplx
    {
        float R;  ///< Real
        float I;  ///< Imaginary
    } stSCplx, *stpSCplx;

    /**
        Double precision complex datatype.
    */
    typedef struct _sDCplx
    {
        double R;  ///< Real
        double I;  ///< Imaginary
    } stDCplx, *stpDCplx;

    /**
        Single precision polar datatype.
    */
    typedef struct _sSPolar
    {
        float M;  ///< Magnitude
        float P;  ///< Phase
    } stSPolar, *stpSPolar;

    /**
        Double precision polar datatype.
    */
    typedef struct _sDPolar
    {
        double M;  ///< Magnitude
        double P;  ///< Phase
    } stDPolar, *stpDPolar;

    /**
        Single precision cartesian/polar datatype.
    */
    typedef union _uSCoord
    {
        stSCplx C;   ///< Cartesian
        stSPolar P;  ///< Polar
    } utSCoord, *utpSCoord;

    /**
        Double precision cartesian/polar datatype.
    */
    typedef union _uDCoord
    {
        stDCplx C;   ///< Cartesian
        stDPolar P;  ///< Polar
    } utDCoord, *utpDCoord;

    #pragma pack()

#endif
