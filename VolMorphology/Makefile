include ../Makefile.common

TARGET=VolMorphology.$(DEXT)
SRCS=morphology.cpp
OBJS=morphology.o
ifeq "$(PLATFORM)" "mingw32"
OBJS+=res.o
endif

include ../Makefile.target

morphology.o: $(SRCS) $(CORE)

res.o: VolMorphology.rc resource.h
	$(WRC) -I. -o $@ $<

