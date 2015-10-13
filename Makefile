TARGET=VolFilters VolMorphology VolProcess VolTensors

all:
	@for t in $(TARGET); do	\
		make -C $$t;		\
	done

clean:
	@for t in $(TARGET); do	\
		make -C $$t clean;	\
	done

distclean:
	@for t in $(TARGET); do		\
		make -C $$t distclean;	\
	done

