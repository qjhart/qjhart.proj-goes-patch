Name:           proj
Version:        4.8.0
Release:        3%{?dist}
Summary:        Cartographic projection software (PROJ.4)

Group:          Applications/Engineering
License:        MIT
URL:            http://proj.osgeo.org
Source0:        http://download.osgeo.org/proj/proj-%{version}.tar.gz
Source1:        http://download.osgeo.org/proj/proj-datumgrid-1.5.zip
Patch0:		proj-4.8.0-removeinclude.patch
Patch1:		proj-4.8.0-add_goes.patch
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires:  libtool


%description
Proj and invproj perform respective forward and inverse transformation of
cartographic data to or from cartesian data with a wide range of selectable
projection functions.


%package devel
Summary:        Development files for PROJ.4
Group:          Development/Libraries
Requires:       %{name} = %{version}-%{release}

%description devel
This package contains libproj and the appropriate header files and man pages.


%package static
Summary:        Development files for PROJ.4
Group:          Development/Libraries

%description static
This package contains libproj static library.


%package nad
Summary:        US and Canadian datum shift grids for PROJ.4
Group:          Applications/Engineering
Requires:       %{name} = %{version}-%{release}

%description nad
This package contains additional US and Canadian datum shift grids.


%package epsg
Summary:        EPSG dataset for PROJ.4
Group:          Applications/Engineering
Requires:       %{name} = %{version}-%{release}

%description epsg
This package contains additional EPSG dataset.

%prep
%setup -q
%patch0 -p0
%patch1 -p0

# disable internal libtool to avoid hardcoded r-path
for makefile in `find . -type f -name 'Makefile.in'`; do
sed -i 's|@LIBTOOL@|%{_bindir}/libtool|g' $makefile
done

# Prepare nad
cd nad
unzip %{SOURCE1}
cd ..
# fix shebag header of scripts
for script in `find nad/ -type f -perm -a+x`; do
sed -i -e '1,1s|:|#!/bin/bash|' $script
done

%build

# fix version info to respect new ABI
sed -i -e 's|5\:4\:5|6\:4\:6|' src/Makefile*

%configure
make OPTIMIZE="$RPM_OPT_FLAGS" %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
%makeinstall
install -p -m 0644 nad/pj_out27.dist nad/pj_out83.dist nad/td_out.dist $RPM_BUILD_ROOT%{_datadir}/%{name}
install -p -m 0755 nad/test27 nad/test83 nad/testvarious $RPM_BUILD_ROOT%{_datadir}/%{name}
install -p -m 0644 nad/epsg $RPM_BUILD_ROOT%{_datadir}/%{name}

# Install projects.h manually, per #830496:
install -p -m 0644 src/projects.h $RPM_BUILD_ROOT%{_includedir}/

%check
pushd nad
# set test enviroment for porj
export PROJ_LIB=$RPM_BUILD_ROOT%{_datadir}/%{name}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH%{buildroot}%{_libdir}
# run tests for proj
./test27      $RPM_BUILD_ROOT%{_bindir}/%{name} || exit 0
./test83      $RPM_BUILD_ROOT%{_bindir}/%{name} || exit 0
./testIGNF    $RPM_BUILD_ROOT%{_bindir}/%{name} || exit 0
./testntv2    $RPM_BUILD_ROOT%{_bindir}/%{name} || exit 0
./testvarious $RPM_BUILD_ROOT%{_bindir}/%{name} || exit 0
popd

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%doc NEWS AUTHORS COPYING README ChangeLog
%{_bindir}/*
%{_mandir}/man1/*.1*
%{_libdir}/*.so.*

%files devel
%defattr(-,root,root,-)
%{_mandir}/man3/*.3*
%{_includedir}/*.h
%{_libdir}/*.so
%exclude %{_libdir}/*.a
%exclude %{_libdir}/libproj.la

%files static
%defattr(-,root,root,-)
%{_libdir}/*.a
%{_libdir}/libproj.la


%files nad
%defattr(-,root,root,-)
%doc nad/README
%attr(0755,root,root) %{_datadir}/%{name}/test27
%attr(0755,root,root) %{_datadir}/%{name}/test83
%attr(0755,root,root) %{_datadir}/%{name}/testvarious
%attr(0755,root,root) %{_libdir}/pkgconfig/%{name}.pc
%exclude %{_datadir}/%{name}/epsg
%{_datadir}/%{name}

%files epsg
%defattr(-,root,root,-)
%doc nad/README
%attr(0644,root,root) %{_datadir}/%{name}/epsg

%changelog
* Thu Aug 16 2012 Devrim GÜNDÜZ <devrim@gunduz.org> 4.8.0-3
- Install projects.h manually, per #830496.

* Sat Jul 21 2012 Fedora Release Engineering <rel-eng@lists.fedoraproject.org> - 4.8.0-2
- Rebuilt for https://fedoraproject.org/wiki/Fedora_18_Mass_Rebuild

* Fri Apr 20 2012 Devrim GÜNDÜZ <devrim@gunduz.org> 4.8.0-1
- Update to 4.8.0, per bz #814851

* Sat Jan 14 2012 Fedora Release Engineering <rel-eng@lists.fedoraproject.org> - 4.7.0-5
- Rebuilt for https://fedoraproject.org/wiki/Fedora_17_Mass_Rebuild

* Wed Feb 09 2011 Fedora Release Engineering <rel-eng@lists.fedoraproject.org> - 4.7.0-4
- Rebuilt for https://fedoraproject.org/wiki/Fedora_15_Mass_Rebuild

* Thu Mar 18 2010 Balint Cristian <cristian.balint@gmail.com> - 4.7.0-3
- fix for bz#562671

* Thu Mar 18 2010 Balint Cristian <cristian.balint@gmail.com> - 4.7.0-2
- fix for bz#556091

* Fri Dec 4 2009 Devrim GÜNDÜZ <devrim@gunduz.org> 4.7.0-1
- Update to 4.7.0
- Update to new datumgrid (1.5)
