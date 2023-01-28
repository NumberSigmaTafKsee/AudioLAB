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


/* f1: 0.25,  f2: 0.25625, r: 0.1, a: 96, t: elliptic */

/*static const long lDec2hpIIRSize = 8l;

static const float fpDec2hpIIRCoeffs[][5] = {
{ // 0
    0.0636055269017405,
    -0.120237273165707,
    0.0636055269017405,
    -1.10059257576193,
    -0.34804090273112
},
{ // 1
    0.197395100757952,
    -0.251955931876418,
    0.197395100757952,
    -0.853316818943548,
    -0.500062952335869
},
{ // 2
    0.416867536074146,
    -0.291504346387221,
    0.416867536074146,
    -0.557073607566324,
    -0.682313026101835
},
{ // 3
    0.630732106186422,
    -0.220804780514429,
    0.630732106186422,
    -0.336225336671185,
    -0.818494329558458
},
{ // 4
    0.786029157621284,
    -0.12850692420049,
    0.786029157621284,
    -0.201609680048344,
    -0.902174919491402
},
{ // 5
    0.881212080327745,
    -0.059820066128943,
    0.881212080327745,
    -0.127415820574785,
    -0.949660047359217
},
{ // 6
    0.933417750217718,
    -0.0197780000165422,
    0.933417750217718,
    -0.0898376822944958,
    -0.976451182746475
},
{ // 7
    0.958297484517192,
    -0.00209738478453132,
    0.958297484517192,
    -0.0744385769165409,
    -0.993130930735456
}
};


static const double dpDec2hpIIRCoeffs[][5] = {
{ // 0
    0.0636055269017405,
    -0.120237273165707,
    0.0636055269017405,
    -1.10059257576193,
    -0.34804090273112
},
{ // 1
    0.197395100757952,
    -0.251955931876418,
    0.197395100757952,
    -0.853316818943548,
    -0.500062952335869
},
{ // 2
    0.416867536074146,
    -0.291504346387221,
    0.416867536074146,
    -0.557073607566324,
    -0.682313026101835
},
{ // 3
    0.630732106186422,
    -0.220804780514429,
    0.630732106186422,
    -0.336225336671185,
    -0.818494329558458
},
{ // 4
    0.786029157621284,
    -0.12850692420049,
    0.786029157621284,
    -0.201609680048344,
    -0.902174919491402
},
{ // 5
    0.881212080327745,
    -0.059820066128943,
    0.881212080327745,
    -0.127415820574785,
    -0.949660047359217
},
{ // 6
    0.933417750217718,
    -0.0197780000165422,
    0.933417750217718,
    -0.0898376822944958,
    -0.976451182746475
},
{ // 7
    0.958297484517192,
    -0.00209738478453132,
    0.958297484517192,
    -0.0744385769165409,
    -0.993130930735456
}
};*/


/* f1: 0.25,  f2: 0.25625, r: 0.1, a: 120, t: elliptic */

static const long lDec2hpIIRSize = 10l;

static const float fpDec2hpIIRCoeffs[][5] = {
{ // 0
    0.0860883552986557,
    -0.147454170614253,
    0.0860883552986557,
    -1.13368419454638,
    -0.453315075757941
},
{ // 1
    0.230138726053263,
    -0.265710323710514,
    0.230138726053263,
    -0.862775228186637,
    -0.588763004003677
},
{ // 2
    0.423332055996104,
    -0.290943741886257,
    0.423332055996104,
    -0.58841934809913,
    -0.726027201977596
},
{ // 3
    0.607542248579554,
    -0.233820176523996,
    0.607542248579554,
    -0.381064569119293,
    -0.829969242802397
},
{ // 4
    0.749928014882174,
    -0.154162195076208,
    0.749928014882174,
    -0.244691516068237,
    -0.898709740908793
},
{ // 5
    0.846567103431297,
    -0.0869479806425301,
    0.846567103431297,
    -0.161336299588397,
    -0.94141848709352
},
{ // 6
    0.906783436809245,
    -0.0409329790075081,
    0.906783436809245,
    -0.112953857197872,
    -0.96745370982387
},
{ // 7
    0.941531243802369,
    -0.0138152282315159,
    0.941531243802369,
    -0.0868712028407657,
    -0.98374891867702
},
{ // 8
    0.958970880382364,
    -0.001484413533216,
    0.958970880382364,
    -0.0756521632465086,
    -0.995078337544452
},
{ // 9
    0.186490176427825,
    -0.186490176427825,
    0,
    -0.627019647144349,
    0
}
};


static const double dpDec2hpIIRCoeffs[][5] = {
{ // 0
    0.0860883552986557,
    -0.147454170614253,
    0.0860883552986557,
    -1.13368419454638,
    -0.453315075757941
},
{ // 1
    0.230138726053263,
    -0.265710323710514,
    0.230138726053263,
    -0.862775228186637,
    -0.588763004003677
},
{ // 2
    0.423332055996104,
    -0.290943741886257,
    0.423332055996104,
    -0.58841934809913,
    -0.726027201977596
},
{ // 3
    0.607542248579554,
    -0.233820176523996,
    0.607542248579554,
    -0.381064569119293,
    -0.829969242802397
},
{ // 4
    0.749928014882174,
    -0.154162195076208,
    0.749928014882174,
    -0.244691516068237,
    -0.898709740908793
},
{ // 5
    0.846567103431297,
    -0.0869479806425301,
    0.846567103431297,
    -0.161336299588397,
    -0.94141848709352
},
{ // 6
    0.906783436809245,
    -0.0409329790075081,
    0.906783436809245,
    -0.112953857197872,
    -0.96745370982387
},
{ // 7
    0.941531243802369,
    -0.0138152282315159,
    0.941531243802369,
    -0.0868712028407657,
    -0.98374891867702
},
{ // 8
    0.958970880382364,
    -0.001484413533216,
    0.958970880382364,
    -0.0756521632465086,
    -0.995078337544452
},
{ // 9
    0.186490176427825,
    -0.186490176427825,
    0,
    -0.627019647144349,
    0
}
};
