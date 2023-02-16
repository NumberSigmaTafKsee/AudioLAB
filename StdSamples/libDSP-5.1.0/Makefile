#
# Makefile for libDSP
# Copyright (C) 2000-2005 Jussi Laako
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

CC = gcc
CXX = g++

# gcc 3.0.x
#CFLAGS = -march=athlon -mcpu=athlon -O6 -fomit-frame-pointer -ffast-math -funroll-loops -malign-loops=5 -malign-jumps=5 -malign-functions=5 -mpreferred-stack-boundary=5 -mfancy-math-387 -Wall -Werror #-g
#OPTCFLAGS = -march=athlon -O3 -ffast-math -funroll-loops -malign-loops=5 -malign-jumps=5 -malign-functions=5 -mpreferred-stack-boundary=5 -mfancy-math-387 -Wall -Werror #-g
# gcc >= 3.1.x
#  -tbird, -4, -xp, -mp
CFLAGS = -march=pentium2 -mcpu=pentium4 -mmmx -msse -msse2
#CFLAGS = -march=athlon-tbird -mcpu=athlon-xp -m3dnow
#CFLAGS = -march=pentium3 -mcpu=pentium3 -msse -mfpmath=sse
#CFLAGS = -march=pentium4 -mcpu=pentium4 -msse2 -mfpmath=sse -m128bit-long-double
# gcc >= 3.4.x
#CFLAGS = -march=pentium2 -mtune=pentium4 -mmmx
#CFLAGS = -march=athlon-tbird -mtune=athlon-xp -m3dnow
#CFLAGS = -march=pentium3 -mtune=pentium3 -msse -mfpmath=sse
#CFLAGS = -march=pentium4 -mtune=pentium4 -msse2 -mfpmath=sse -m128bit-long-double
#CFLAGS = -march=pentium4 -O2 -mfpmath=387
#CFLAGS = -march=prescott -mtune=prescott -msse3 -mfpmath=sse -m128bit-long-double
CFLAGS += -O3 -ffast-math -funroll-loops -fprefetch-loop-arrays
CFLAGS += -Wall -Werror -g
#CFLAGS += -ftracer
# gcc >= 4.x
#CFLAGS += -fmodulo-sched
#CFLAGS += -ftree-vectorize
#
#CFLAGS += -ggdb -pg
#
#
# gcc >= 3.1.x
OPTCFLAGS = -march=pentium2 -mcpu=pentium4 -mmmx -msse -msse2
#OPTCFLAGS = -march=athlon-tbird -mcpu=athlon-xp -m3dnow 
#OPTCFLAGS = -march=pentium3 -mcpu=pentium3 -msse
#OPTCFLAGS = -march=pentium4 -mcpu=pentium4 -msse2 -m128bit-long-double
# gcc >= 3.4.x
#OPTCFLAGS = -march=pentium2 -mtune=pentium4 -mmmx
#OPTCFLAGS = -march=athlon-tbird -mtune=athlon-xp -m3dnow 
#OPTCFLAGS = -march=pentium3 -mtune=pentium3 -msse
#OPTCFLAGS = -march=pentium4 -mtune=pentium4 -msse2 -m128bit-long-double
#OPTCFLAGS = -march=prescott -mtune=prescott -msse3 -m128bit-long-double
OPTCFLAGS += -O1 -ffast-math -fprefetch-loop-arrays #-funroll-loops
OPTCFLAGS += -Wall -Werror -g
#OPTCFLAGS += -ftracer
#
#OPTCFLAGS += -ggdb -pg
CXXFLAGS = $(CFLAGS)

# thread safety
DEFS = -D_REENTRANT -D_THREAD_SAFE
# define if your compiler/c-library is ISO C 9x standard compliant
DEFS += -D_ISOC9X_SOURCE
# define if you are using GNU C++-compiler and GNU C-library
DEFS += -D_GNU_SOURCE
# use memmove() for copying data instead of built-in routine, for gcc inlining
#DEFS += -DUSE_MEMMOVE
# define for x86 platform specific optimizations
DEFS += -DDSP_X86
# define for Intel Integrated Performance Primitives support
#DEFS += -DDSP_IPP

INCS = -I. #-I/usr/local/include
#INCS += `pkg-config --cflags fftw3f fftw3`
#INCS += -I/opt/intel/ipp/include

PREFIX = /usr/local
TGPREFIX = $(PREFIX)

CSRC = c
CCSRC = cc
OBJ = o

SRCS = DSPOp.$(CCSRC) \
	FFTDecimator.$(CCSRC) \
	FFTInterpolator.$(CCSRC) \
	FFTMultiRate.$(CCSRC) \
	Filter.$(CCSRC) \
	FIRDecimator.$(CCSRC) \
	FIRInterpolator.$(CCSRC) \
	FIRMultiRate.$(CCSRC) \
	FlipBand.$(CCSRC) \
	Hankel.$(CCSRC) \
	IIRCascade.$(CCSRC) \
	IIRDecimator.$(CCSRC) \
	IIRInterpolator.$(CCSRC) \
	IIRMultiRate.$(CCSRC) \
	ReBuffer.$(CCSRC) \
	RecDecimator.$(CCSRC) \
	RecInterpolator.$(CCSRC) \
	Transform4.$(CCSRC) \
	Transform8.$(CCSRC) \
	TransformS.$(CCSRC) \
	cdspop.$(CCSRC) \
	X86.$(CSRC)

OBJS = DSPOp.$(OBJ) \
	FFTDecimator.$(OBJ) \
	FFTInterpolator.$(OBJ) \
	FFTMultiRate.$(OBJ) \
	Filter.$(OBJ) \
	FIRDecimator.$(OBJ) \
	FIRInterpolator.$(OBJ) \
	FIRMultiRate.$(OBJ) \
	FlipBand.$(OBJ) \
	Hankel.$(OBJ) \
	IIRCascade.$(OBJ) \
	IIRDecimator.$(OBJ) \
	IIRInterpolator.$(OBJ) \
	IIRMultiRate.$(OBJ) \
	ReBuffer.$(OBJ) \
	RecDecimator.$(OBJ) \
	RecInterpolator.$(OBJ) \
	Transform4.$(OBJ) \
	Transform8.$(OBJ) \
	TransformS.$(OBJ) \
	cdspop.$(OBJ) \
	X86.$(OBJ)

default: all

all: libdsp.a

.c.o: $(SRCS)
	$(CC) $(CFLAGS) $(DEFS) $(INCS) -c $<

.cc.o: $(SRCS)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCS) -c $<

X86.$(OBJ): X86.$(CSRC)
	$(CC) $(OPTCFLAGS) $(DEFS) $(INCS) -c $<

libdsp.a: $(OBJS)
	ar rc libdsp.a $(OBJS)
	ranlib libdsp.a

install: libdsp.a
	install -m 755 -d $(PREFIX)/include/dsp
	install -m 644 dsp/DSPConfig.hh $(PREFIX)/include/dsp
	install -m 644 dsp/DSPOp.hh $(PREFIX)/include/dsp
	install -m 644 dsp/DSPVector.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FFTDecimator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FFTInterpolator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FFTMultiRate.hh $(PREFIX)/include/dsp
	install -m 644 dsp/Filter.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FIRDecimator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FIRInterpolator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FIRMultiRate.hh $(PREFIX)/include/dsp
	install -m 644 dsp/FlipBand.hh $(PREFIX)/include/dsp
	install -m 644 dsp/Hankel.hh $(PREFIX)/include/dsp
	install -m 644 dsp/IIRCascade.hh $(PREFIX)/include/dsp
	install -m 644 dsp/IIRDecimator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/IIRInterpolator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/IIRMultiRate.hh $(PREFIX)/include/dsp
	install -m 644 dsp/ReBuffer.hh $(PREFIX)/include/dsp
	install -m 644 dsp/ReBufferT.hh $(PREFIX)/include/dsp
	install -m 644 dsp/RecDecimator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/RecInterpolator.hh $(PREFIX)/include/dsp
	install -m 644 dsp/RingBuffer.hh $(PREFIX)/include/dsp
	install -m 644 dsp/Transform4.hh $(PREFIX)/include/dsp
	install -m 644 dsp/Transform8.hh $(PREFIX)/include/dsp
	install -m 644 dsp/TransformS.hh $(PREFIX)/include/dsp
	install -m 644 dsp/dsptypes.h $(PREFIX)/include/dsp
	install -m 644 dsp/dspop.h $(PREFIX)/include/dsp
	install -m 644 dsp/X86.h $(PREFIX)/include/dsp
	install -m 755 -d $(PREFIX)/lib
	install -m 644 libdsp.a $(PREFIX)/lib

libdsp.dep: $(SRCS)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCS) -MM $(SRCS) >libdsp.dep

clean:
	rm -f *.$(OBJ) *~ libdsp.a core*

cleanall:
	rm -rf *.$(OBJ) *~ libdsp.dep libdsp.a core*

include libdsp.dep

