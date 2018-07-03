# Install script for directory: /home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local/HDF_Group/HDF5/1.8.18")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "headers")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/hdf5.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5api_adpt.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5public.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5version.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5overflow.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Apkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Apublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5ACpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5ACpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5B2pkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5B2public.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Bpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Bpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Dpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Dpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Edefin.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Einit.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Epkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Epubgen.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Epublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Eterm.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Fpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Fpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDcore.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDdirect.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDfamily.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDlog.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDmpi.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDmpio.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDmulti.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDsec2.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FDstdio.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FSpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5FSpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Gpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Gpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5HFpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5HFpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5HGpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5HGpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5HLpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5HLpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5MPpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Opkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Opublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Oshared.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Ppkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Ppublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5PLextern.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5PLpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Rpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Rpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Spkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Spublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5SMpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Tpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Tpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Zpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Zpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Cpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Cpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Ipkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Ipublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Lpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Lpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5MMpublic.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Rpkg.h;/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5/H5Rpublic.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/annika/COCO/build/lib/gatb-core/gatb-core/include/hdf5" TYPE FILE FILES
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/hdf5.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5api_adpt.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5public.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5version.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5overflow.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Apkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Apublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5ACpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5ACpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5B2pkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5B2public.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Bpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Bpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Dpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Dpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Edefin.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Einit.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Epkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Epubgen.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Epublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Eterm.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Fpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Fpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDcore.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDdirect.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDfamily.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDlog.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDmpi.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDmpio.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDmulti.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDsec2.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FDstdio.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FSpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5FSpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Gpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Gpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5HFpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5HFpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5HGpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5HGpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5HLpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5HLpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5MPpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Opkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Opublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Oshared.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Ppkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Ppublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5PLextern.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5PLpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Rpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Rpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Spkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Spublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5SMpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Tpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Tpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Zpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Zpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Cpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Cpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Ipkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Ipublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Lpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Lpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5MMpublic.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Rpkg.h"
    "/home/annika/COCO/lib/gatb-core/gatb-core/thirdparty/hdf5/src/H5Rpublic.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "libraries")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/annika/COCO/build/lib/gatb-core/gatb-core/lib/libhdf5.a")
endif()

