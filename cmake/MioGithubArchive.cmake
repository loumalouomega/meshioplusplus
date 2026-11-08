# mio_github_archive(<out_var> REPO <owner/name> TAG <tag> PREFIX <dir/> SHA256 <hash>)
#
# Sets <out_var> to a local copy of https://github.com/<repo>/archive/refs/tags/<tag>.tar.gz
# whose SHA-256 is <hash>, for a FetchContent `URL`. GitHub's archive downloads
# have gone down on their own for hours (503/504, November 2026) while clones
# kept working. GitHub builds a tag archive with `git archive` piped through
# `gzip -cn`, so when the download fails the same command on a shallow clone of
# the tag rebuilds the archive's exact bytes. Either copy is refused unless its
# SHA-256 is <hash>. <dir/> is the archive's top directory, GitHub's
# `<name>-<tag without a leading v>/`.

function(mio_github_archive out_var)
  cmake_parse_arguments(PARSE_ARGV 1 arg "" "REPO;TAG;PREFIX;SHA256" "")
  string(REPLACE "/" "-" stem "${arg_REPO}-${arg_TAG}")
  set(dir "${CMAKE_BINARY_DIR}/_deps/mio-archives")
  set(archive "${dir}/${stem}.tar.gz")
  set(url "https://github.com/${arg_REPO}/archive/refs/tags/${arg_TAG}.tar.gz")

  if(EXISTS "${archive}")
    file(SHA256 "${archive}" got)
    if(got STREQUAL arg_SHA256)
      set(${out_var} "${archive}" PARENT_SCOPE)
      return()
    endif()
    file(REMOVE "${archive}")
  endif()
  file(MAKE_DIRECTORY "${dir}")

  file(DOWNLOAD "${url}" "${archive}.part" STATUS status TLS_VERIFY ON)
  list(GET status 0 code)
  if(code EQUAL 0)
    file(SHA256 "${archive}.part" got)
    if(got STREQUAL arg_SHA256)
      file(RENAME "${archive}.part" "${archive}")
      set(${out_var} "${archive}" PARENT_SCOPE)
      return()
    endif()
    message(WARNING "meshio++: ${url} has SHA-256 ${got}, expected ${arg_SHA256}")
  else()
    list(GET status 1 reason)
    message(STATUS "meshio++: downloading ${url} failed (${reason}); rebuilding it from a git clone")
  endif()
  file(REMOVE "${archive}.part")

  find_package(Git QUIET)
  find_program(MIO_GZIP_EXECUTABLE gzip)
  if(NOT GIT_FOUND OR NOT MIO_GZIP_EXECUTABLE)
    message(FATAL_ERROR "meshio++: cannot download ${url}, and rebuilding it needs git and gzip")
  endif()
  set(clone "${dir}/${stem}.git")
  file(REMOVE_RECURSE "${clone}")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" clone --quiet --bare --depth 1 --branch "${arg_TAG}"
            "https://github.com/${arg_REPO}.git" "${clone}"
    RESULT_VARIABLE rc)
  if(NOT rc EQUAL 0)
    file(REMOVE_RECURSE "${clone}")
    message(FATAL_ERROR "meshio++: cannot download ${url} or clone ${arg_REPO}")
  endif()
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${clone}" -c "tar.tar.gz.command=${MIO_GZIP_EXECUTABLE} -cn"
            archive --format=tar.gz "--prefix=${arg_PREFIX}" "${arg_TAG}"
    OUTPUT_FILE "${archive}.part"
    RESULT_VARIABLE rc)
  file(REMOVE_RECURSE "${clone}")
  if(rc EQUAL 0)
    file(SHA256 "${archive}.part" got)
  endif()
  if(NOT rc EQUAL 0 OR NOT got STREQUAL arg_SHA256)
    file(REMOVE "${archive}.part")
    message(FATAL_ERROR "meshio++: the archive of ${arg_REPO} ${arg_TAG} rebuilt from git does not match SHA-256 ${arg_SHA256}")
  endif()
  file(RENAME "${archive}.part" "${archive}")
  message(STATUS "meshio++: rebuilt ${url} from a git clone (SHA-256 matches)")
  set(${out_var} "${archive}" PARENT_SCOPE)
endfunction()
