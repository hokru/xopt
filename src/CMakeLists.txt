# handle target and source files


# NOTE: This helper function assumes no generator expressions are used
#       for the source files
function(target_sources_local target)
  if(POLICY CMP0076)
    # New behavior is available, so just forward to it by ensuring
    # that we have the policy set to request the new behavior, but
    # don't change the policy setting for the calling scope
    cmake_policy(PUSH)
    cmake_policy(SET CMP0076 NEW)
    target_sources(${target} ${ARGN})
    cmake_policy(POP)
    return()
  endif()

  # Must be using CMake 3.12 or earlier, so simulate the new behavior
  message('cmake <=3.12 detected')
  unset(_srcList)
  get_target_property(_targetSourceDir ${target} SOURCE_DIR)

  foreach(src ${ARGN})
    if(NOT src STREQUAL "PRIVATE" AND
       NOT src STREQUAL "PUBLIC" AND
       NOT src STREQUAL "INTERFACE" AND
       NOT IS_ABSOLUTE "${src}")
      # Relative path to source, prepend relative to where target was defined
      file(RELATIVE_PATH src "${_targetSourceDir}" "${CMAKE_CURRENT_LIST_DIR}/${src}")
    endif()
    list(APPEND _srcList ${src})
  endforeach()
  target_sources(${target} ${_srcList})
endfunction()



target_sources_local(xopt
    PRIVATE
 string.f90
 Hmass.f90
 TRproj.f90
 anc.f90
 ancopt.f90
 ciopt.f90
 conjgrad.f90
 constrains.f90
 control.f90
 copt.f90
 defaults.f90
 eval_opt.f90
 gdiis.f90
 getHess.f90
 getgrad.F90
 hessup.f90
 intcoords.f90
 intopt.f90
 io.f90
 irc_opt.f90
 lambda.f90
 math.f90
 md.f90
 molecule.f90
 oniom.f90
 printmat.f90
 progs.f90
 random.f90
 tools.f90
 hbonds.f90
 modules.F90 main.F90 print.F90
)

