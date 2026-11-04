# Options renamed in PCMS 0.4. For one release the old names are still
# honored: a value given under the old name is copied to the new name, with a
# deprecation message. Remove this file in the release after 0.4.

# pcms_renamed_option(<old> <new>)
#   If <old> is set (cache or -D), store its value in the cache under <new>
#   and drop <old> from the cache. If both are set to different values, <new>
#   wins and a warning is printed. Must be called before option(<new> ...).
function(pcms_renamed_option old new)
  if(NOT DEFINED ${old})
    return()
  endif()
  if(DEFINED ${new} AND NOT "${${new}}" STREQUAL "${${old}}")
    message(WARNING "Both ${old} (deprecated) and ${new} are set; "
                    "using ${new}=${${new}}.")
  else()
    message(DEPRECATION "${old} is deprecated and will be removed in the next "
                        "release; use ${new} instead.")
    set(${new} "${${old}}" CACHE BOOL "" FORCE)
  endif()
  unset(${old} CACHE)
endfunction()

pcms_renamed_option(PCMS_ENABLE_Python PCMS_ENABLE_PYTHON)
pcms_renamed_option(PCMS_ENABLE_Fortran PCMS_ENABLE_FORTRAN)
pcms_renamed_option(PETSC_LINK_STATIC PCMS_PETSC_LINK_STATIC)
