cmake_minimum_required(VERSION 3.15.0...3.21.0)
file(
  GLOB bpDirs
  LIST_DIRECTORIES TRUE
  *.bp*)
message(STATUS "before bpDirs ${bpDirs}")
if(bpDirs)
  file(REMOVE_RECURSE ${bpDirs})
  file(
    GLOB bpDirs
    LIST_DIRECTORIES TRUE
    *.bp*)
  message(STATUS "after bpDirs ${bpDirs}")
endif()
file(
  GLOB sstDirs
  LIST_DIRECTORIES TRUE
  *.sst*)
message(STATUS "before sstDirs ${sstDirs}")
if(sstDirs)
  file(REMOVE_RECURSE ${sstDirs})
  file(
    GLOB sstDirs
    LIST_DIRECTORIES TRUE
    *.sst*)
  message(STATUS "after sstDirs ${sstDirs}")
endif()
