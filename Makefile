LDLIBS = -lm -lgmp
CFLAGS = -Wall -g -O3

ALL: model model_ejovo validate

model: model.o harmonics.o
model_ejovo: model_ejovo.o harmonics.o

validate: validate.o harmonics.o 
model.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f model validate *.o