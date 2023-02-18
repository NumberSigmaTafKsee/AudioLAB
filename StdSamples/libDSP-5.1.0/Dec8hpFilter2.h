/*

    Decimation filter for rate 8
    Copyright (C) 2001 Jussi Laako

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


#ifdef _MSC_VER
#pragma warning(disable:4305)
#endif


static const long lDec8hpFilterLen = 103l;

static const float fDec8hpFilterGain = 1.0f;
static const float fpDec8hpFilterCoeffs[] = {
    -1.65190872271921E-5,
    1.15333328006912E-5,
    5.05216408135781E-6,
    -3.66969878882348E-5,
    8.56114367969341E-5,
    -0.00015190812611348,
    0.000232810442646285,
    -0.000322050222020069,
    0.000409592149747397,
    -0.000481809852296002,
    0.000522206795238215,
    -0.000512724532045543,
    0.000435614639584877,
    -0.000275774400256122,
    2.33676459091394E-5,
    0.000323519650787846,
    -0.000756494713748482,
    0.00125507937875774,
    -0.00178586889793034,
    0.00230291757924441,
    -0.00274955980510631,
    0.00306175830622723,
    -0.00317292682499798,
    0.00302001324743053,
    -0.00255046689760477,
    0.00172956688156007,
    -0.000547474333516481,
    -0.000974694188754834,
    0.002780479802087,
    -0.00477602540768809,
    0.0068308631131305,
    -0.00878186176247722,
    0.0104404052586862,
    -0.0116025944645654,
    0.0120619793181363,
    -0.0116240596241431,
    0.0101215692794949,
    -0.00742940360263099,
    0.00347798222999672,
    0.00173612708535739,
    -0.00814336195237868,
    0.0155989001163807,
    -0.023885463033512,
    0.0327213005634033,
    -0.0417730653922508,
    0.050672887600792,
    -0.059038601938188,
    0.0664957952544986,
    -0.0727001560425869,
    0.0773585408161661,
    -0.0802472316605976,
    0.0812260431843405,
    -0.0802472316605976,
    0.0773585408161661,
    -0.0727001560425869,
    0.0664957952544986,
    -0.059038601938188,
    0.050672887600792,
    -0.0417730653922508,
    0.0327213005634033,
    -0.023885463033512,
    0.0155989001163807,
    -0.00814336195237868,
    0.00173612708535739,
    0.00347798222999672,
    -0.00742940360263099,
    0.0101215692794949,
    -0.0116240596241431,
    0.0120619793181363,
    -0.0116025944645654,
    0.0104404052586862,
    -0.00878186176247722,
    0.0068308631131305,
    -0.00477602540768809,
    0.002780479802087,
    -0.000974694188754834,
    -0.000547474333516481,
    0.00172956688156007,
    -0.00255046689760477,
    0.00302001324743053,
    -0.00317292682499798,
    0.00306175830622723,
    -0.00274955980510631,
    0.00230291757924441,
    -0.00178586889793034,
    0.00125507937875774,
    -0.000756494713748482,
    0.000323519650787846,
    2.33676459091394E-5,
    -0.000275774400256122,
    0.000435614639584877,
    -0.000512724532045543,
    0.000522206795238215,
    -0.000481809852296002,
    0.000409592149747397,
    -0.000322050222020069,
    0.000232810442646285,
    -0.00015190812611348,
    8.56114367969341E-5,
    -3.66969878882348E-5,
    5.05216408135781E-6,
    1.15333328006912E-5,
    -1.65190872271921E-5
};


static const double dDec8hpFilterGain = 1.0;
static const double dpDec8hpFilterCoeffs[] = {
    -1.65190872271921E-5,
    1.15333328006912E-5,
    5.05216408135781E-6,
    -3.66969878882348E-5,
    8.56114367969341E-5,
    -0.00015190812611348,
    0.000232810442646285,
    -0.000322050222020069,
    0.000409592149747397,
    -0.000481809852296002,
    0.000522206795238215,
    -0.000512724532045543,
    0.000435614639584877,
    -0.000275774400256122,
    2.33676459091394E-5,
    0.000323519650787846,
    -0.000756494713748482,
    0.00125507937875774,
    -0.00178586889793034,
    0.00230291757924441,
    -0.00274955980510631,
    0.00306175830622723,
    -0.00317292682499798,
    0.00302001324743053,
    -0.00255046689760477,
    0.00172956688156007,
    -0.000547474333516481,
    -0.000974694188754834,
    0.002780479802087,
    -0.00477602540768809,
    0.0068308631131305,
    -0.00878186176247722,
    0.0104404052586862,
    -0.0116025944645654,
    0.0120619793181363,
    -0.0116240596241431,
    0.0101215692794949,
    -0.00742940360263099,
    0.00347798222999672,
    0.00173612708535739,
    -0.00814336195237868,
    0.0155989001163807,
    -0.023885463033512,
    0.0327213005634033,
    -0.0417730653922508,
    0.050672887600792,
    -0.059038601938188,
    0.0664957952544986,
    -0.0727001560425869,
    0.0773585408161661,
    -0.0802472316605976,
    0.0812260431843405,
    -0.0802472316605976,
    0.0773585408161661,
    -0.0727001560425869,
    0.0664957952544986,
    -0.059038601938188,
    0.050672887600792,
    -0.0417730653922508,
    0.0327213005634033,
    -0.023885463033512,
    0.0155989001163807,
    -0.00814336195237868,
    0.00173612708535739,
    0.00347798222999672,
    -0.00742940360263099,
    0.0101215692794949,
    -0.0116240596241431,
    0.0120619793181363,
    -0.0116025944645654,
    0.0104404052586862,
    -0.00878186176247722,
    0.0068308631131305,
    -0.00477602540768809,
    0.002780479802087,
    -0.000974694188754834,
    -0.000547474333516481,
    0.00172956688156007,
    -0.00255046689760477,
    0.00302001324743053,
    -0.00317292682499798,
    0.00306175830622723,
    -0.00274955980510631,
    0.00230291757924441,
    -0.00178586889793034,
    0.00125507937875774,
    -0.000756494713748482,
    0.000323519650787846,
    2.33676459091394E-5,
    -0.000275774400256122,
    0.000435614639584877,
    -0.000512724532045543,
    0.000522206795238215,
    -0.000481809852296002,
    0.000409592149747397,
    -0.000322050222020069,
    0.000232810442646285,
    -0.00015190812611348,
    8.56114367969341E-5,
    -3.66969878882348E-5,
    5.05216408135781E-6,
    1.15333328006912E-5,
    -1.65190872271921E-5
};