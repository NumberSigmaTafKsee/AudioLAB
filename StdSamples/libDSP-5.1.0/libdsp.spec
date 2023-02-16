#
# RPM package spec for SuSE
#
# Copyright (C) 2004-2005 Jussi Laako <jussi@sonarnerd.net>
#

%define version 5.0.2
%define release 1.suse

%define pbuildroot %{_tmppath}/libdsp-%{version}-root

Name: libdsp
Summary: C/C++ library of digital signal processing functions
Version: %{version}
Release: %{release}
License: GPL
Vendor: Jussi Laako <jussi@sonarnerd.net>
Distribution: research software
# Copyright: Copyright (C) 1998-2005 Jussi Laako
Packager: Jussi Laako <jussi@sonarnerd.net>
Source: libdsp-%{version}.tar.bz2
Group: Development/Libraries/C and C++
BuildRoot: %{pbuildroot}

BuildRequires: gcc-c++
BuildRequires: libtool

%description
libDSP is a library of commonly used digital signal processing functions for
use in C++ language. Wrapper is provided for use in C language.


%package devel
Summary: Development tools required for developing software using libDSP
Group: Development/Libraries/C and C++

Requires: gcc-c++

%description devel
The libdsp-devel package contains files required to develop programs that
use libDSP digital signal processing library.


%package doc
Summary: Documentation for the libDSP library
Group: Documentation

%description doc
The libdsp-doc package contains documentation for developing applications
using libDSP digital signal processing library.

%prep
%setup -q -n libdsp-%{version}

%build
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/%{_prefix}/include
mkdir -p $RPM_BUILD_ROOT/%{_prefix}/lib
make -C Inlines PREFIX=$RPM_BUILD_ROOT%{_prefix} install
make -C DynThreads PREFIX=$RPM_BUILD_ROOT%{_prefix} \
	PKGINC=-I$RPM_BUILD_ROOT%{_prefix}/include \
	CXX=g++ CXXFLAGS="$RPM_OPT_FLAGS -ffast-math"
make -C DynThreads PREFIX=$RPM_BUILD_ROOT%{_prefix} install
make -C libDSP PREFIX=$RPM_BUILD_ROOT%{_prefix} \
	PKGINC=-I$RPM_BUILD_ROOT%{_prefix}/include \
	TGPREFIX=%{_prefix} \
	CC=gcc CXX=g++ \
	CFLAGS="$RPM_OPT_FLAGS -ffast-math" \
	OPTCFLAGS="$RPM_OPT_FLAGS -ffast-math" \
	CXXFLAGS="$RPM_OPT_FLAGS -ffast-math"

%install
make -C libDSP PREFIX=$RPM_BUILD_ROOT%{_prefix} TGPREFIX=%{_prefix} install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc CONTACT COPYING CREDITS
%{_prefix}/lib/libdsp.so.*

%files devel
%defattr(-,root,root,-)
%{_prefix}/lib/libDynThreads.a
%{_prefix}/lib/libdsp.so
%{_prefix}/lib/libdsp.a
%{_prefix}/lib/libdsp.la
%{_prefix}/include/Alloc.hh
%{_prefix}/include/Compilers.hh
%{_prefix}/include/Condition.hh
%{_prefix}/include/Exception.hh
%{_prefix}/include/Generic.hh
%{_prefix}/include/Mutex.hh
%{_prefix}/include/PthCond.hh
%{_prefix}/include/PthMutex.hh
%{_prefix}/include/RWLock.hh
%{_prefix}/include/Semaphore.hh
%{_prefix}/include/DynThreads.hh
%{_prefix}/include/dsp/*

%files doc
%defattr(-,root,root,-)
%doc html libdsp-doc.ps libdsp-doc.pdf

%changelog
* Fri May 13 2005 Jussi Laako <jussi@sonarnerd.net>
- Updated to a new version

* Wed Feb 23 2005 Jussi Laako <jussi@sonarnerd.net>
- Updated to a new version

* Sun Nov 07 2004 Jussi Laako <jussi@sonarnerd.net>
- Some changes to .spec for SuSE

* Tue Mar 02 2004 Jussi Laako <jussi@sonarnerd.net>
- Separated the documentation to own package

* Sun Sep 21 2003 Jussi Laako <jussi@sonarnerd.net>
- Create .spec file for rpm package
- First rpm build

