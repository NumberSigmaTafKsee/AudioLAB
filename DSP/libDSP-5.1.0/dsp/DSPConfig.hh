/** \mainpage

    DSP operations
    
    Copyright (C) 1998-2004 Jussi Laako

    \author Jussi Laako
    $Date: 2005/12/17 21:29:49 $

    \par Overview
    This set of classes implement routines commonly used in digital signal
    processing.

    \par Licensing
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    \par
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    \par
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#ifndef DSPCONFIG_HH
#define DSPCONFIG_HH

    // If neither is defined, then split-radix transforms are used
    #undef DSP_USE_RADIX4
    #undef DSP_USE_RADIX8

    // Define to parallelize large FFTs
    #undef USE_CDFT_THREADS

    // Define to use FFTW
    #undef DSP_USE_FFTW
    
    #define DSP_FILT_DEF_OVERLAPF   0.75f
    #define DSP_FILT_DEF_OVERLAP    0.75
    #define DSP_FILT_DEF_BETAF      14.96454265f
    #define DSP_FILT_DEF_BETA       14.96454265

#endif
