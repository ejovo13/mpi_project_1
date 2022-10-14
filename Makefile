INC = -I/usr/local/include/ejovo -I/usr/local/include/ejovo/matrix
LDLIBS = -lm -lgmp -lejovo
CFLAGS = -Wall -g -O3 -Wno-unused-variable $(INC) -Werror


ALL: model model_ejovo validate

model: model.o harmonics.o
model_ejovo: model_ejovo.o harmonics.o

validate: validate.o harmonics.o 
model.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

# model_ejovo.o: model_ejovo.c
# 	gcc -c model_ejovo.c $(CFLAGS) $(LDLIBS)

.PHONY: clean

clean:
	rm -f model validate *.o