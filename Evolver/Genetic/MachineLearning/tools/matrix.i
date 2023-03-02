/** 
 @cond
 ############################################################################
 # LGPL License                                                             #
 #                                                                          #
 # This file is part of the Machine Learning Framework.                     #
 # Copyright (c) 2010-2012, Philipp Kraus, <philipp.kraus@flashpixx.de>     #
 # This program is free software: you can redistribute it and/or modify     #
 # it under the terms of the GNU Lesser General Public License as           #
 # published by the Free Software Foundation, either version 3 of the       #
 # License, or (at your option) any later version.                          #
 #                                                                          #
 # This program is distributed in the hope that it will be useful,          #
 # but WITHOUT ANY WARRANTY; without even the implied warranty of           #
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
 # GNU Lesser General Public License for more details.                      #
 #                                                                          #
 # You should have received a copy of the GNU Lesser General Public License #
 # along with this program. If not, see <http://www.gnu.org/licenses/>.     #
 ############################################################################
 @endcond
 **/

/** interface file for the matrix calls **/


#ifdef SWIGJAVA
%module "matrixmodule"
%include "../swig/java/java.i"
%rename(Matrix) matrix;

%typemap(javabody)           machinelearning::tools::matrix ""
%typemap(javafinalize)       machinelearning::tools::matrix ""
%typemap(javadestruct)       machinelearning::tools::matrix ""
#endif


%nodefaultctor               machinelearning::tools::matrix;
%nodefaultdtor               machinelearning::tools::matrix;


%include "matrix.hpp"
%template(max) machinelearning::tools::matrix::max<double>;
%template(min) machinelearning::tools::matrix::min<double>;
%template(mean) machinelearning::tools::matrix::mean<double>;
%template(variance) machinelearning::tools::matrix::variance<double>;
%template(sum) machinelearning::tools::matrix::sum<double>;
%template(trace) machinelearning::tools::matrix::trace<double>;
%template(cov) machinelearning::tools::matrix::cov<double>;