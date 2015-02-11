# File automatically generated: DO NOT EDIT.

# Get the exported informations about the targets
get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
#include(${_dir}/FCCSWExports.cmake)

# Set useful properties
get_filename_component(_dir "${_dir}" PATH)
set(FCCSW_INCLUDE_DIRS ${_dir}/include)
set(FCCSW_LIBRARY_DIRS ${_dir}/lib)

set(FCCSW_BINARY_PATH ${_dir}/bin ${_dir}/scripts)
set(FCCSW_PYTHON_PATH ${_dir}/python)

set(FCCSW_COMPONENT_LIBRARIES )
set(FCCSW_LINKER_LIBRARIES )

set(FCCSW_ENVIRONMENT PREPEND;PATH;/usr/bin;PREPEND;PATH;/afs/cern.ch/exp/fcc/sw/0.2/Python/2.7.6/x86_64-slc6-gcc48-dbg/bin;PREPEND;PATH;\${.}/scripts;PREPEND;PATH;\${.}/bin;PREPEND;LD_LIBRARY_PATH;\${.}/lib;PREPEND;ROOT_INCLUDE_PATH;\${.}/include;PREPEND;PYTHONPATH;\${.}/python;PREPEND;PYTHONPATH;\${.}/python/lib-dynload;SET;FCCSW_PROJECT_ROOT;\${.}/../)

set(FCCSW_EXPORTED_SUBDIRS)
foreach(p )
  get_filename_component(pn ${p} NAME)
  if(EXISTS ${_dir}/cmake/${pn}Export.cmake)
    set(FCCSW_EXPORTED_SUBDIRS ${FCCSW_EXPORTED_SUBDIRS} ${p})
  endif()
endforeach()

set(FCCSW_OVERRIDDEN_SUBDIRS )
