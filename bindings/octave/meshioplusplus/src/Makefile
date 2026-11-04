# Octave `pkg install` build of the MEX gateway (meshioplusplus_mex).
#
# The C library is located like the R binding's configure script does:
#   1. MESHIOPLUSPLUS_CFLAGS / MESHIOPLUSPLUS_LIBS  -- fully manual
#   2. MESHIOPLUSPLUS_HOME                          -- an install prefix
#   3. pkg-config meshioplusplus                    -- the .pc the C API installs
# MATLAB (and Octave without pkg) builds the same sources with inst/mio_build.m.

MKOCTFILE ?= mkoctfile

ifneq ($(strip $(MESHIOPLUSPLUS_CFLAGS)$(MESHIOPLUSPLUS_LIBS)),)
  MIO_CFLAGS := $(MESHIOPLUSPLUS_CFLAGS)
  MIO_LIBS := $(MESHIOPLUSPLUS_LIBS)
else ifneq ($(strip $(MESHIOPLUSPLUS_HOME)),)
  MIO_CFLAGS := -I$(MESHIOPLUSPLUS_HOME)/include
  MIO_LIBS := -L$(MESHIOPLUSPLUS_HOME)/lib -lmeshioplusplus
else
  MIO_CFLAGS := $(shell pkg-config --cflags meshioplusplus 2>/dev/null)
  MIO_LIBS := $(shell pkg-config --libs meshioplusplus 2>/dev/null)
endif

SRC := $(wildcard *.c)

all: meshioplusplus_mex.mex

meshioplusplus_mex.mex: $(SRC) mio_oct.h
	@if [ -z "$(MIO_LIBS)" ]; then \
	  echo "meshioplusplus: the meshio++ C library was not found." >&2; \
	  echo "  Build and install it (build/configure.sh --c-api --build;" >&2; \
	  echo "  cmake --install <tree> --prefix <p>), then set MESHIOPLUSPLUS_HOME=<p>" >&2; \
	  echo "  or put <p>/lib/pkgconfig on PKG_CONFIG_PATH." >&2; \
	  exit 1; \
	fi
	$(MKOCTFILE) --mex $(MIO_CFLAGS) $(SRC) $(MIO_LIBS) -o $@

clean:
	rm -f *.o *.mex

.PHONY: all clean
