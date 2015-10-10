include ../Makefile.common

TARGET=VolFilters.$(DEXT)

LIBS=-lfftw3f -lfftw3
ifeq "$(PLATFORM)" "linux"
LIBS+=-lfftw3f_threads -lfftw3_threads
endif

SRCS=volfilters.cpp convolution.cpp
OBJS=volfilters.o
ifeq "$(PLATFORM)" "mingw32"
OBJS+=res.o
endif

include ../Makefile.target

volfilters.o: $(SRCS) $(CORE)

res.o: VolFilters.rc resource.h
	$(WRC) -I. -o $@ $<

