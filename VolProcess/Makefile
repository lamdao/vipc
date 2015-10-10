include ../Makefile.common

TARGET=VolProcess.$(DEXT)
SRCS=volprocess.cpp
OBJS=volprocess.o
ifeq "$(PLATFORM)" "mingw32"
OBJS+=res.o
endif

include ../Makefile.target

volprocess.o: $(SRCS) $(CORE)

res.o: VolProcess.rc resource.h
	$(WRC) -I. -o $@ $<

