include ../Makefile.common

TARGET=VolTensors.$(DEXT)

SRCS=voltensors.cpp eig3.cpp
OBJS=voltensors.o eig3.o
ifeq "$(PLATFORM)" "mingw32"
OBJS+=res.o
endif

include ../Makefile.target

voltensors.o: $(SRCS) $(CORE)

res.o: VolTensors.rc resource.h
	$(WRC) -I. -o $@ $<

