/*
* Copyright (c) 2015-2020 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser    General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.    If not, see <http://www.gnu.org/licenses/>.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dspl.h"


#ifdef DOXYGEN_ENGLISH
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int array_scale_lin(double* x, int n,
                        double xmin, double xmax, double dx,
                        double h, double* y)
\brief Vector `x` linear transformation

Function transforms values \f$x(i)\f$, \f$i = 0,1,\ldots n\f$
to the \f$y(i)\f$, accordint to equation:

\f[
y(i) = k_x x(i) + d_x, \qquad k_x =
\frac{h}{x_{\textrm{max}} - x_{\textrm{min}}}.
\f]

All values of the vector `x` between
\f$x_{\textrm{min}}\f$ and \f$x_{\textrm{max}}\f$, transforms to
the vector `y` between \f$d_x\f$ and \f$h + d_x\f$.
Parameter \f$d_x\f$ sets mean shift of the vector `y`.

This function is convenient for translating values โโ
of different dimensions. For example it can be used
to transfer the values โโof the vector `x`
to the graph of the height of` h`, where the height can
be set in the number of pixels, in centimeters, etc.

\param[in] x
Pointer to the input vector `x`. \n
Vector size is `[n x 1]`. \n
\n

\param[in] n
Size of vector `x`. \n
\n

\param[in] xmin
Parameter \f$x_{\textrm{min}}\f$. \n
\n

\param[in] xmax
Parameter \f$x_{\textrm{min}}\f$. \n
Value `xmax` must be more than `xmin`. \n
\n

\param[in] dx
Displacement after transformation. \n
This parameter must have output vector `y`
dimensions (pixels, centimeters). \n
\n

\param[in] h
Height of vector `y` after transforming between `dx` and `h+dx`. \n
\n

\param[out] y
Pointer to the output vector `y`. \n
Vector size is `[n x 1]`. \n
Memory must be allocated. \n
\note
Pointer    `y` can be equal to `x`.
Velues of vector `x` will be rewritten in this case. \n
\n

\return
`RES_OK` if function returns successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif

#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int array_scale_lin(double* x, int n,
                        double xmin, double xmax, double dx,
                        double h, double* y)
\brief ะะธะฝะตะนะฝะพะต ัะฐัััะถะตะฝะธะต ะฒะตะบัะพัะฐ ะดะฐะฝะฝัั `x`
ะคัะฝะบัะธั ะฟัะพะธะทะฒะพะดะธั ะฟัะตะพะฑัะฐะทะพะฒะฐะฝะธะต ะทะฝะฐัะตะฝะธะน \f$x(i)\f$, \f$i = 0,1,\ldots n\f$
ะฒ ะทะฝะฐัะตะฝะธั \f$y(i)\f$, ะฒ ัะพะพัะฒะตัััะฒะธะธ ั ัะพัะผัะปะพะน:

\f[
y(i) = k_x x(i) + d_x, \qquad k_x =
\frac{h}{x_{\textrm{max}} - x_{\textrm{min}}}.
\f]

ะขะฐะบะธะผ ะพะฑัะฐะทะพะผ, ะฒัะต ะทะฝะฐัะตะฝะธั ะฒัะพะดะฝะพะณะพ ะฒะตะบัะพัะฐ `x` ะฒ ะดะธะฐะฟะฐะทะพะฝะต ะพั
\f$x_{\textrm{min}}\f$ ะดะพ \f$x_{\textrm{max}}\f$, ะปะธะฝะตะนะฝะพ ัะฐัััะณะธะฒะฐัััั ะฒ
ะทะฝะฐัะตะฝะธั ะฒะตะบัะพัะฐ `y` ะฒ ะดะธะฐะฟะฐะทะพะฝะต ะพั \f$d_x\f$ ะดะพ \f$h + d_x\f$.
ะะฐะผะตัะธะผ, ััะพ \f$d_x\f$ ะทะฐะดะฐะตั ะปะธะฝะตะนะฝะพะต ัะผะตัะตะฝะธะต ะทะฝะฐัะตะฝะธะน ะฒะตะบัะพัะฐ `y`.

ะะฐะฝะฝะฐั ััะฝะบัะธั ัะดะพะฑะฝะฐ ะดะปั ะฟะตัะตะฒะพะดะฐ ะฒะตะปะธัะธะฝ ัะฐะทะฝัั ัะฐะทะผะตัะฝะพััะตะน, ะฒ ัะฐััะฝะพััะธ,
ะดะปั ะฟะตัะตะฝะพัะฐ ะทะฝะฐัะตะฝะธะน ะฒะตะบัะพัะฐ `x` ะฝะฐ ะณัะฐัะธะบ ะฒััะพัั `h`, ะณะดะต ะฒััะพัะฐ ะผะพะถะตั
ะฑััั ะทะฐะดะฐะฝะฐ ะฒ ะบะพะปะธัะตััะฒะต ะฟะธะบัะตะปะตะน, ะฒ ัะฐะฝัะธะผะตััะฐั ะธ ั.ะด.

\param[in] x
ะฃะบะฐะทะฐัะตะปั ะฝะฐ ะฒะตะบัะพั ะฒัะพะดะฝัั ะทะฝะฐัะตะฝะธะน `x`. \n
ะ?ะฐะทะผะตั ะฒะตะบัะพัะฐ `[n x 1]`. \n
\n

\param[in] n
ะ?ะฐะทะผะตั ะฒะตะบัะพัะฐ `x`. \n
\n

\param[in] xmin
ะะธะถะฝัั ะณัะฐะฝะธัะฐ ะดะธะฐะฟะฐะทะพะฝะฐ ััะฐะฝััะพัะผะฐัะธะธ. \n
\n

\param[in] xmax
ะะตััะฝัั ะณัะฐะฝะธัะฐ ะดะธะฐะฟะฐะทะพะฝะฐ ััะฐะฝััะพัะผะฐัะธะธ. \n
ะะฝะฐัะตะฝะธะต `xmax` ะดะพะปะถะฝะพ ะฑััั ัััะพะณะพ ะฑะพะปััะต ะทะฝะฐัะตะฝะธั `xmin`. \n
\n

\param[in] dx
ะกะผะตัะตะฝะธะต ะฟะพัะปะต ััะฐะฝััะพัะผะฐัะธะธ. \n
ะะฐะฝะฝัะน ะฟะฐัะฐะผะตัั ะดะพะปะถะตะฝ ะธะผะตัั ัะฐะทะผะตัะฝะพััั ะฒััะพะดะฝะพะณะพ ะฒะตะบัะพัะฐ `y`. \n
\n

\param[in] h
ะะธะฐะฟะฐะทะพะฝ ะทะฝะฐัะตะฝะธะน ะฒะตะบัะพัะฐ `y` ะฟะพัะปะต ััะฐะฝััะพัะผะฐัะธะธ ะพั `dx` ะดะพ `h+dx`. \n
\n

\param[out] y
ะฃะบะฐะทะฐัะตะปั ะฝะฐ ะฒะตะบัะพัะฐ ะดะฐะฝะฝัั ะฟะพัะปะต ััะฐะฝััะพัะผะฐัะธะธ. \n
ะ?ะฐะทะผะตั ะฒะตะบัะพัะฐ `[n x 1]`. \n
ะะฐะผััั ะดะพะปะถะฝะฐ ะฑััั ะฒัะดะตะปะตะฝะฐ. \n
\note
ะฃะบะฐะทะฐัะตะปั `y` ะผะพะถะตั ัะพะฒะฟะฐะดะฐัั ั `x`, ะฒ ััะพะผ ัะปััะฐะต,
ะดะฐะฝะฝัะต ะฒะตะบัะพัะฐ `x` ะฑัะดัั ะฟะตัะตะทะฐะฟะธัะฐะฝั ะปะธะฝะตะนะฝะพ ะธะทะผะตะฝะตะฝะฝัะผะธ ะฒ ัะพะพัะฒะตัััะฒะธะธ
ั ัะพัะผัะปะพะน ะฒััะต. \n
\n

\return
`RES_OK` ะตัะปะธ ััะฝะบัะธั ะฒัะฟะพะปะฝะตะฝะฐ ััะฟะตัะฝะพ. \n
 ะ ะฟัะพัะธะฒะฝะพะผ ัะปััะฐะต \ref ERROR_CODE_GROUP "ะบะพะด ะพัะธะฑะบะธ".

\author
ะะฐัััะธะฝ ะกะตัะณะตะน
www.dsplib.org
***************************************************************************** */
#endif

int DSPL_API array_scale_lin(double* x,   int n,
                             double xmin, double xmax, double dx,
                             double h,    double* y)
{
    double kx;
    int k;
    if(!x)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(h<0.0)
        return ERROR_NEGATIVE;

    if(xmin >= xmax)
        return ERROR_MIN_MAX;

    kx = h / (xmax - xmin);

    for(k = 0; k < n; k++)
        y[k] = (x[k] - xmin) * kx + dx;

    return RES_OK;
}
