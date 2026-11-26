# Options renamed in PCMS 0.4. For one release the old names are still
# honored: a value given under the old name is copied to the new name, with a
# deprecation message. Remove this file in the release after 0.4.

# pcms_renamed_option(<old> <new>)
#   If <old> is set (cache or -D), store its value in the cache under <new>
#   and drop <old> from the cache. Because <old> is dropped every time, it is
#   only set when it was just given (-D on this configure), so it also wins
#   over a <new> value cached by an earlier configure. Must be called before
#   option(<new> ...).
function(pcms_renamed_option old new)
  if(NOT DEFINED ${old})
    return()
  endif()
  message(DEPRECATION "${old} is deprecated and will be removed in the next "
                      "release; use ${new} instead.")
  if(DEFINED ${new} AND NOT "${${new}}" STREQUAL "${${old}}")
    message(WARNING "${old}=${${old}} replaces ${new}=${${new}}.")
  endif()
  set(${new} "${${old}}" CACHE BOOL "" FORCE)
  unset(${old} CACHE)
endfunction()

pcms_renamed_option(PCMS_ENABLE_Python PCMS_ENABLE_PYTHON)
pcms_renamed_option(PCMS_ENABLE_Fortran PCMS_ENABLE_FORTRAN)
pcms_renamed_option(PETSC_LINK_STATIC PCMS_PETSC_LINK_STATIC)
pcms_renamed_option(PCMS_ENABLE_MESHFIELDS PCMS_ENABLE_OMEGA_H_MESHFIELDS)
# the projections are always built since 0.5; what the flag added on top of
# OMEGA_H_MESHFIELDS (which it required) is the PETSc solve
pcms_renamed_option(PCMS_ENABLE_CONSERVATIVE_PROJECTION PCMS_ENABLE_PETSC)

# pcms_removed_option(<old> <reason>)
#   Options that no longer exist: if <old> is set, say that it has no effect
#   (for one release) and drop it from the cache.
function(pcms_removed_option old reason)
  if(DEFINED ${old})
    message(DEPRECATION "${old} is removed and has no effect (${reason}).")
    unset(${old} CACHE)
  endif()
endfunction()
