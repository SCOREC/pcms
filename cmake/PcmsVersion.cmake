# Writes pcms/version.h with the version derived from git. Run in script mode
# (cmake -P) at build time, so the version follows the checked out commit
# without reconfiguring. configure_file only rewrites the header when its
# contents change, so files including it are not rebuilt on every build.
#
# Inputs: SOURCE_DIR, GIT_EXECUTABLE, PROJECT_VERSION, INPUT, OUTPUT, and
# optionally VERSION_OVERRIDE
#
# The version follows PEP 440 with the setuptools-scm conventions:
#   at tag vX.Y.Z:         X.Y.Z
#   N commits past it:     X.Y.(Z+1).devN+g<hash>
#   uncommitted changes:   X.Y.(Z+1).devN+g<hash>.dYYYYMMDD
# When no vX.Y.Z tag is reachable (e.g. a shallow clone) X.Y.Z is the
# project() version and N is 0. Without git history (e.g. a source tarball) it
# falls back to the project() version and the revision is "unknown". VERSION_OVERRIDE, when set, replaces
# the version; pip builds pass the version setuptools-scm computed for the
# package metadata, so the two always agree (also when building an sdist).

string(REPLACE "." ";" components "${PROJECT_VERSION}")
list(GET components 0 PROJECT_VERSION_MAJOR)
list(GET components 1 PROJECT_VERSION_MINOR)
list(GET components 2 PROJECT_VERSION_PATCH)

set(PCMS_VERSION_STRING "${PROJECT_VERSION}")
set(PCMS_GIT_REVISION "unknown")

if(GIT_EXECUTABLE)
  # Only use git when SOURCE_DIR is the top of the repository, so a tarball
  # unpacked inside some other repository does not take that one's version.
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --show-toplevel
    WORKING_DIRECTORY "${SOURCE_DIR}"
    OUTPUT_VARIABLE toplevel
    RESULT_VARIABLE result
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  get_filename_component(toplevel "${toplevel}" REALPATH)
  get_filename_component(source_dir "${SOURCE_DIR}" REALPATH)
  if(result EQUAL 0 AND toplevel STREQUAL source_dir)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
      WORKING_DIRECTORY "${SOURCE_DIR}"
      OUTPUT_VARIABLE revision
      RESULT_VARIABLE result
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(result EQUAL 0)
      set(PCMS_GIT_REVISION "${revision}")
    endif()
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" describe --tags --long --dirty
              --match "v[0-9]*"
      WORKING_DIRECTORY "${SOURCE_DIR}"
      OUTPUT_VARIABLE describe
      RESULT_VARIABLE result
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    # e.g. v0.3.0-511-g8a9fa572 or v0.3.0-511-g8a9fa572-dirty
    if(result EQUAL 0 AND describe MATCHES
       "^v([0-9]+)\\.([0-9]+)\\.([0-9]+)-([0-9]+)-g([0-9a-f]+)(-dirty)?$")
      set(major "${CMAKE_MATCH_1}")
      set(minor "${CMAKE_MATCH_2}")
      set(patch "${CMAKE_MATCH_3}")
      set(distance "${CMAKE_MATCH_4}")
      set(hash "${CMAKE_MATCH_5}")
      set(dirty "${CMAKE_MATCH_6}")
      set(tagged TRUE)
    else()
      # No release tag is reachable, e.g. in a shallow clone. Count from the
      # project() version with an unknown distance (0), so the build still
      # reports a development version rather than the release.
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --always --dirty
        WORKING_DIRECTORY "${SOURCE_DIR}"
        OUTPUT_VARIABLE describe
        RESULT_VARIABLE result
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if(result EQUAL 0 AND describe MATCHES "^([0-9a-f]+)(-dirty)?$")
        set(major "${PROJECT_VERSION_MAJOR}")
        set(minor "${PROJECT_VERSION_MINOR}")
        set(patch "${PROJECT_VERSION_PATCH}")
        set(distance 0)
        set(hash "${CMAKE_MATCH_1}")
        set(dirty "${CMAKE_MATCH_2}")
        set(tagged FALSE)
      endif()
    endif()
    if(DEFINED hash)
      if(tagged AND distance EQUAL 0 AND NOT dirty)
        set(PCMS_VERSION_STRING "${major}.${minor}.${patch}")
      else()
        math(EXPR patch "${patch} + 1")
        set(PCMS_VERSION_STRING
            "${major}.${minor}.${patch}.dev${distance}+g${hash}")
        if(dirty)
          string(TIMESTAMP today "%Y%m%d")
          string(APPEND PCMS_VERSION_STRING ".d${today}")
        endif()
      endif()
    endif()
  endif()
endif()

if(VERSION_OVERRIDE)
  set(PCMS_VERSION_STRING "${VERSION_OVERRIDE}")
endif()

configure_file("${INPUT}" "${OUTPUT}" @ONLY)
