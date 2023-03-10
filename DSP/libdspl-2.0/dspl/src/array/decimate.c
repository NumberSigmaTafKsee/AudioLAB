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
\fn int decimate(double* x, int n, int d, double* y, int* cnt)
\brief
Real vector decimation

Function `d` times decimates real vector `x`. \n
Output vector `y` keeps values corresponds to:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
Pointer to the input real vector `x`. \n
Vector `x` size is `[n x 1]`. \n \n

\param[in] n
Size of input vector `x`. \n \n

\param[in] d
Decimation coefficient. \n
Each d-th vector will be copy from vector `x` to the
output vector `y`. \n \n

\param[out] y
Pointer to the output decimated vector `y`. \n
Output vector size is `[n/d x 1]` will be copy
to the address `cnt`. \n

\param[out] cnt
Address which will keep decimated vector `y` size. \n
Pointer can be `NULL`, vector `y` will not return
in this case. \n \n

\return
`RES_OK` if function calculated successfully. \n
Else \ref ERROR_CODE_GROUP "code error".

Two-times decimation example:
\code{.cpp}
double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double y[5];
int d = 2;
int cnt;

decimate(x, 10, d, y, &cnt);
\endcode
As result variable `cnt` will be written value 5 and
vector `y` will keep    array:
\verbatim
c = [0.0, 2.0, 4.0, 6.0, 8.0]
\endverbatim

\author Sergey Bakhurin www.dsplib.org
***************************************************************************** */
#endif
#ifdef DOXYGEN_RUSSIAN
/*! ****************************************************************************
\ingroup ARRAY_GROUP
\fn int decimate(double* x, int n, int d, double* y, int* cnt)
\brief ?????????????????? ?????????????????????????? ?????????????? ????????????

?????????????? ???????????????????? ?????????????????? ?????????????????????????? ?????????????? `x` ?? `d` ??????. \n
?? ???????????????????? ???????????????? ???????????? `y` ???????????????? ????????????????:
`y(k) = x(k*d), k = 0...n/d-1` \n

\param[in] x
?????????????????? ???? ???????????? ?????????????? ???????????? `x`. \n
???????????? ?????????????? `[n x 1]`. \n \n

\param[in] n
???????????? ???????????????? ?????????????? `x`. \n \n

\param[in] d
?????????????????????? ??????????????????. \n
?? ???????????????????? ?????????????????? ???? ?????????????? `x` ?????????? ?????????? ????????????
d-?? ??????????????. \n \n

\param[out] y
?????????????????? ???? ???????????????????????????? ???????????? `y`. \n
???????????? ?????????????????? ?????????????? ?????????? `[n/d x 1]`
?????????? ???????????????? ???? ???????????? `cnt`. \n
???????????? ???????????? ???????? ????????????????. \n \n

\param[out] cnt
?????????????????? ????????????????????, ?? ?????????????? ?????????? ????????????????
???????????? ?????????????????? ?????????????? ?????????? ??????????????????. \n
?????????????????? ?????????? ???????? `NULL`, ?? ???????? ????????????
???????????? ?????????????? `y` ???? ????????????????????????. \n \n

\return
`RES_OK` ???????? ?????????????? ?????????????????? ??????????????. \n
 ?? ?????????????????? ???????????? \ref ERROR_CODE_GROUP "?????? ????????????".

???????????? ?????????????????? ?????????????????????????? ?????????????? ???????????? ?? 2 ????????:
\code{.cpp}
double x[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
double y[5];
int d = 2;
int cnt;
decimate(x, 10, d, y, &cnt);
\endcode
?? ???????????????????? ?? ???????????????????? `cnt` ?????????? ?????????????? ???????????? 5,
?? ???????????? `y` ?????????? ?????????????? ???????????? ????????????:
\verbatim
c = [0.0, 2.0, 4.0, 6.0, 8.0]
\endverbatim

\author
?????????????? ????????????
www.dsplib.org
**************************************************************************** */
#endif
int DSPL_API decimate(double* x, int n, int d, double* y, int* cnt)
{
    int k = 0, i = 0;
    if(!x || !y)
        return ERROR_PTR;
    if(n < 1)
        return ERROR_SIZE;
    if(d < 1)
        return ERROR_NEGATIVE;

    k = i = 0;
    while(k + d <= n)
    {
        y[i] = x[k];
        k+=d;
        i++;
    }
    if(cnt)
        *cnt = i;

    return RES_OK;
}

