LDLIBS = -lm
CFLAGS = -Wall -g -O3

ALL: model validate

test_harmonic: harmonics.f90 test_harmonic.f90
	gfortran $^ -o test


model: model.o harmonics.o
validate: validate.o harmonics.o
model.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f model validate *.o