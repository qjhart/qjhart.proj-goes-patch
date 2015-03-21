# Introduction

This project supplies a patch for the inclusion of NOAA GOES
satellites into the PROJ.4 http://trac.osgeo.org/proj/ library and
binary.  This allows projections to and from GOES using the standard
PROJ.4 utilities (e.g. cs2cs), and also allows for the inclusion of
other software that use the libproj library for projections, like
GRASS.

Debian and RedHat comes with a standard proj package.  This is the package that
is used for projecting between coordinate systems in grass, among
other software packages.  A patched version of the software needs to
applied to the systems, for the GOES projection.

# Software

This patch adds a new projection +proj=goes, corresponding to NOAA
GOES satellites.  The following flags are associated with the the
projection:
* +proj=goes
* +goes=15
  * (or +goes=1506 for summer flip ~MAR20-SEP20)
  * (or +goes=1512 for winter flip ~ SEP20-MAR20 )
* +nsinc
* +ewinc

if you have *gvar_inspector* installed, you can get the projection
information form a GVAR file, with the *--proj* command-line.

# Compiling

## debian

The debian patch uses quilting,
http://wiki.debian.org/UsingQuilt#Using_quilt_with_Debian_source_packages
as the methodology for adding the patch.  Note that below the
qjhart.proj-goes-patch/proj-4.7.0/ directory contains some files in
the src/ directory, but those are for info only, their source is
included in the patch, and the patch is the cannonical format.

```
#!bash
# Checkout the proj patch, get location.
gh=~/qjhart.proj-goes-patch
# Get the proj source
apt-get source proj
# Copy apply the patch and build
(cd proj-4.7.0; quilt import $gh/add-goes-4.7; dpkg-buildpackage)
# install these sources
sudo dpkg --install *proj*.deb
```

## RedHat

We can use the same patch for the redhat patch.
http://bradthemad.org/tech/notes/patching_rpms.php provides some insight on
how to do this.

```
#!bash
# Checkout the proj patch
git clone https://github.com/qjhart/qjhart.proj-goes-patch.git
# First get the proj source and build environment
mkdir -p ~/rpmbuild/{BUILD,RPMS,SOURCES,SPECS,SRPMS}
echo '%_topdir %(echo $HOME)/rpmbuild' > ~/.rpmmacros
#Setup your build environment ( this uses home directory)
sudo yum install rpm-build redhat-rpm-config
# Get the proj source
yumdownloader --source proj
rpm -ivh proj-4.7.0-2_0.el6.src.rpm
# Copy the patch to the rpmbuild directory
cp ~/qjhart.proj-goes-patch/rpmbuild/add-goes-4.7 ~/rpmbuild/SOURCES
```

Now you need to add the patch to the ~/rpmbuild/SPEC/proj.spec file.
You need to add two lines, see below for those diffs.  These are the
changes to the proj.spec files

```{bash}
--- /home/qhart/rpmbuild/SPECS/proj.spec        2011-11-08 11:47:16.000000000 -0800
+++ /home/qhart/rpmbuild/SPECS/proj.spec+       2015-03-20 18:53:38.460888519 -0700
@@ -9,6 +9,7 @@
 Source0:   http://download.osgeo.org/proj/proj-%{version}.tar.gz
 Source1:   http://download.osgeo.org/proj/proj-datumgrid-1.5.zip
 BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
+Patch0: add-goes-4.7

 BuildRequires: libtool
 Requires: %{name}-epsg = %{version}-%{release}
@@ -44,6 +45,7 @@

 %prep
 %setup -q
+%patch0 -p1

 # disable internal libtool to avoid hardcoded r-path
 for makefile in `find . -type f -name 'Makefile.in'`; do
```

Then, build the package (if this fails you made need to install more dependencies)

```{bash}
rpmbuild -ba rpmbuild/SPECS/proj.spec
# install these sources
sudo rpm -Uvh --force ~/rpmbuild/RPMS/x86_64/proj-*.rpm
```

Responds with:
```
Preparing...                ########################################### [100%]
   1:proj-epsg              ########################################### [ 20%]
   2:proj                   ########################################### [ 40%]
   3:proj-devel             ########################################### [ 60%]
   4:proj-nad               ########################################### [ 80%]
   5:proj-debuginfo         ########################################### [100%]
```

You now have a set of packages that can be installed in any matching
redhat distribution.  You can check with 

```
proj +proj=goes +goes=15
```

Should not fail
 
# History

Starting around 2003 the CSTARS program began modeling incoming solar
radiation for California from GOES satellite data.  We manipulate the
rasters using grass.  In order to project the GOES images into our
projection scheme, we wrote a patch to proj, adding a goes projection
In 2012, with the launch of GOES15, we contacted NOAA, and they kindly
sent us their C code from projecting GOES data, though it was more of
a standalone program.  We, have incorporated their code for GOES
instruments beyond GOES-11

# Maintenance

http://pkg-perl.alioth.debian.org/howto/quilt.html

# Acknowledgments

I'd like to that *NOAA* for providing a listing of their GOES projection software routines.
