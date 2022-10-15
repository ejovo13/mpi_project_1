INC = -I/usr/local/include/ejovo -I/usr/local/include/ejovo/matrix
LDLIBS = -lm -lejovo
CFLAGS = -Wall -g -O3 -Wno-unused-variable $(INC) -Werror -Wextra


ALL: model model_ejovo validate reduce

model: model.o harmonics.o
reduce: reduce.o 
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