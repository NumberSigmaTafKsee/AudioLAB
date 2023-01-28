/*

    IIR Decimation filter for rate 2 (highpass)
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


#ifdef _MSC_VER
#pragma warning(disable:4305)
#endif


/* f1: 0.25,  f2: 0.275, r: 0.1, a: 96, t: chebyshev */

static const long lDec2hpIIRSize = 12l;

static const float fpDec2hpIIRCoeffs[][5] = {
{  // 0
    0.42064423485298,
    -0.84128846970596,
    0.42064423485298,
    -0.303652893838038,
    -0.986229833249959
},
{  // 1
    0.40672194829703,
    -0.813443896594059,
    0.40672194829703,
    -0.332018621461108,
    -0.958906414649227
},
{  // 2
    0.384827952051726,
    -0.769655904103452,
    0.384827952051726,
    -0.392022091101754,
    -0.931333899308658
},
{  // 3
    0.355023287751185,
    -0.71004657550237,
    0.355023287751185,
    -0.482883089499192,
    -0.902976240503933
},
{  // 4
    0.317511649981645,
    -0.63502329996329,
    0.317511649981645,
    -0.603428632062601,
    -0.873475231989181
},
{  // 5
    0.272878000762055,
    -0.54575600152411,
    0.272878000762055,
    -0.751232186542352,
    -0.842744189590572
},
{  // 6
    0.222412684936652,
    -0.444825369873304,
    0.222412684936652,
    -0.921453544293719,
    -0.811104284040327
},
{  // 7
    0.168485436189186,
    -0.336970872378372,
    0.168485436189186,
    -1.10550548746505,
    -0.779447232221797
},
{  // 8
    0.114848408005692,
    -0.229696816011384,
    0.114848408005692,
    -1.28997369440667,
    -0.74936732642944
},
{  // 9
    0.0666347442085557,
    -0.133269488417111,
    0.0666347442085557,
    -1.45661356441776,
    -0.723152541251987
},
{  // 10
    0.0297697124561515,
    -0.0595394249123029,
    0.0297697124561515,
    -1.58442226880081,
    -0.703501118625415
},
{  // 11
    0.00968490905923375,
    -0.0193698181184675,
    0.00968490905923375,
    -1.65417213149347,
    -0.692911767730408
}
};


static const double dpDec2hpIIRCoeffs[][5] = {
{  // 0
    0.42064423485298,
    -0.84128846970596,
    0.42064423485298,
    -0.303652893838038,
    -0.986229833249959
},
{  // 1
    0.40672194829703,
    -0.813443896594059,
    0.40672194829703,
    -0.332018621461108,
    -0.958906414649227
},
{  // 2
    0.384827952051726,
    -0.769655904103452,
    0.384827952051726,
    -0.392022091101754,
    -0.931333899308658
},
{  // 3
    0.355023287751185,
    -0.71004657550237,
    0.355023287751185,
    -0.482883089499192,
    -0.902976240503933
},
{  // 4
    0.317511649981645,
    -0.63502329996329,
    0.317511649981645,
    -0.603428632062601,
    -0.873475231989181
},
{  // 5
    0.272878000762055,
    -0.54575600152411,
    0.272878000762055,
    -0.751232186542352,
    -0.842744189590572
},
{  // 6
    0.222412684936652,
    -0.444825369873304,
    0.222412684936652,
    -0.921453544293719,
    -0.811104284040327
},
{  // 7
    0.168485436189186,
    -0.336970872378372,
    0.168485436189186,
    -1.10550548746505,
    -0.779447232221797
},
{  // 8
    0.114848408005692,
    -0.229696816011384,
    0.114848408005692,
    -1.28997369440667,
    -0.74936732642944
},
{  // 9
    0.0666347442085557,
    -0.133269488417111,
    0.0666347442085557,
    -1.45661356441776,
    -0.723152541251987
},
{  // 10
    0.0297697124561515,
    -0.0595394249123029,
    0.0297697124561515,
    -1.58442226880081,
    -0.703501118625415
},
{  // 11
    0.00968490905923375,
    -0.0193698181184675,
    0.00968490905923375,
    -1.65417213149347,
    -0.692911767730408
}
};


/* f1: 0.25,  f2: 0.2625, r: 0.1, a: 96, t: chebyshev */

/*static const long lDec2hpIIRSize = 17l;

static const float fpDec2hpIIRCoeffs[][5] = {
{ // 0
    0.460065603944421,
    -0.920131207888842,
    0.460065603944421,
    -0.152782879484544,
    -0.993045295262228
},
{ // 1
    0.45267708191174,
    -0.90535416382348,
    0.45267708191174,
    -0.168469260361288,
    -0.979177588008248
},
{ // 2
    0.441112442803256,
    -0.882224885606513,
    0.441112442803256,
    -0.20074959065465,
    -0.965199361867676
},
{ // 3
    0.425325013440645,
    -0.85065002688129,
    0.425325013440645,
    -0.249623822437774,
    -0.950923876200354
},
{ // 4
    0.405250618433441,
    -0.810501236866882,
    0.405250618433441,
    -0.315180111963534,
    -0.936182585697298
},
{ // 5
    0.38083671020161,
    -0.761673420403221,
    0.38083671020161,
    -0.397485743203235,
    -0.920832584009678
},
{ // 6
    0.352084328513632,
    -0.704168657027264,
    0.352084328513632,
    -0.496431542620658,
    -0.904768856675186
},
{ // 7
    0.319106376813331,
    -0.638212753626662,
    0.319106376813331,
    -0.611517181968269,
    -0.887942689221593
},
{ // 8
    0.282204668082882,
    -0.564409336165765,
    0.282204668082882,
    -0.741568467790158,
    -0.870387140121687
},
{ // 9
    0.241964751293878,
    -0.483929502587757,
    0.241964751293878,
    -0.884390212388789,
    -0.852249217564303
},
{ // 10
    0.199360268487163,
    -0.398720536974326,
    0.199360268487163,
    -1.03638470035386,
    -0.833825774302511
},
{ // 11
    0.155846576060795,
    -0.311693152121591,
    0.155846576060795,
    -1.19220950813359,
    -0.815595812376769
},
{ // 12
    0.113407883707899,
    -0.226815767415798,
    0.113407883707899,
    -1.34460479387256,
    -0.798236328704152
},
{ // 13
    0.0745092969272763,
    -0.149018593854553,
    0.0745092969272763,
    -1.48456702536589,
    -0.782604213074993
},
{ // 14
    0.0419071362403782,
    -0.0838142724807565,
    0.0419071362403782,
    -1.60203888241858,
    -0.769667427380098
},
{ // 15
    0.0183020342306161,
    -0.0366040684612323,
    0.0183020342306161,
    -1.68717174215049,
    -0.76037987907295
},
{ // 16
    0.00588292592699223,
    -0.0117658518539845,
    0.00588292592699223,
    -1.73198555043042,
    -0.755517254138386
}
};


static const double dpDec2hpIIRCoeffs[][5] = {
{ // 0
    0.460065603944421,
    -0.920131207888842,
    0.460065603944421,
    -0.152782879484544,
    -0.993045295262228
},
{ // 1
    0.45267708191174,
    -0.90535416382348,
    0.45267708191174,
    -0.168469260361288,
    -0.979177588008248
},
{ // 2
    0.441112442803256,
    -0.882224885606513,
    0.441112442803256,
    -0.20074959065465,
    -0.965199361867676
},
{ // 3
    0.425325013440645,
    -0.85065002688129,
    0.425325013440645,
    -0.249623822437774,
    -0.950923876200354
},
{ // 4
    0.405250618433441,
    -0.810501236866882,
    0.405250618433441,
    -0.315180111963534,
    -0.936182585697298
},
{ // 5
    0.38083671020161,
    -0.761673420403221,
    0.38083671020161,
    -0.397485743203235,
    -0.920832584009678
},
{ // 6
    0.352084328513632,
    -0.704168657027264,
    0.352084328513632,
    -0.496431542620658,
    -0.904768856675186
},
{ // 7
    0.319106376813331,
    -0.638212753626662,
    0.319106376813331,
    -0.611517181968269,
    -0.887942689221593
},
{ // 8
    0.282204668082882,
    -0.564409336165765,
    0.282204668082882,
    -0.741568467790158,
    -0.870387140121687
},
{ // 9
    0.241964751293878,
    -0.483929502587757,
    0.241964751293878,
    -0.884390212388789,
    -0.852249217564303
},
{ // 10
    0.199360268487163,
    -0.398720536974326,
    0.199360268487163,
    -1.03638470035386,
    -0.833825774302511
},
{ // 11
    0.155846576060795,
    -0.311693152121591,
    0.155846576060795,
    -1.19220950813359,
    -0.815595812376769
},
{ // 12
    0.113407883707899,
    -0.226815767415798,
    0.113407883707899,
    -1.34460479387256,
    -0.798236328704152
},
{ // 13
    0.0745092969272763,
    -0.149018593854553,
    0.0745092969272763,
    -1.48456702536589,
    -0.782604213074993
},
{ // 14
    0.0419071362403782,
    -0.0838142724807565,
    0.0419071362403782,
    -1.60203888241858,
    -0.769667427380098
},
{ // 15
    0.0183020342306161,
    -0.0366040684612323,
    0.0183020342306161,
    -1.68717174215049,
    -0.76037987907295
},
{ // 16
    0.00588292592699223,
    -0.0117658518539845,
    0.00588292592699223,
    -1.73198555043042,
    -0.755517254138386
}
};*/