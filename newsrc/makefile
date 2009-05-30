#======================== GRFFDE ===========================#
#
#
#CFLAGS = -pedantic -O0
#CC = gcc3
#LIBS = -lm
#
#
# prof doesn't seem to produce much in the way of improvement
CC = icc
#CFLAGS = -O0
#CFLAGS = -O3 -ipo #-unroll -w2 -axiMKW #-prof_use #-prof_gen # -prof_use
#
#LIBS = -lsvml -limf -lm
#CFLAGS = -c99 -g  -O0 -mp -fp -prec_div -DDEBUG -DNUMERICAL
#CFLAGS = -O3 -ipo -w2 -Wall -ansi -c99 -axW -tpp7 -Xc -DGRIDOUTPUT -DDEBUG

#CFLAGS = -O3 -ipo -w2 -Wall -ansi -c99 -axW -tpp7 -Xc -DGRIDOUTPUT
CFLAGS = -O3 -ipo -w2 -Wall

#CFLAGS = -pg
#
#
#CC = gcc
#CFLAGS = -O0 -pedantic -std=gnu99 -Wall \
# -W -Wfloat-equal -Wtraditional -Wundef -Wshadow -Wpointer-arith \
# -Wbad-function-cast -Wcast-qual -Wcast-align -Wwrite-strings \
# -Wconversion -Wsign-compare -Waggregate-return -Wstrict-prototypes \
# -Wmissing-prototypes -Wmissing-declarations -Wmissing-noreturn \
# -Wmissing-format-attribute -Wpacked -Wpadded -Wredundant-decls \
# -Wnested-externs -Wunreachable-code -Winline -Wlong-long \
# -Wdisabled-optimization \
# -DDEBUG
#LIBS = -lm

#
#CC = gcc
#LIBS = -lm
#CFLAGS = -pg

SRCS = \
bdump.c bounds.c coord.c diag.c dump.c ener.c \
init.c interp.c main.c metric.c phys.c \
step_ch.c timestep.c

OBJS = \
bdump.o bounds.o coord.o diag.o dump.o ener.o \
init.o interp.o main.o metric.o phys.o \
step_ch.o timestep.o
 
gr: $(OBJS) makefile
	$(CC) $(CFLAGS) -o gr $(OBJS) $(LIBS)

# dependencies
$(OBJS) : defs.h decs.h makefile

clean:
	rm *.o *.il

dump:
	rm dump??? ener.out

tild:
	rm *.c~ *.h~ *.m~ makefile~

ps :     $(SRCS) defs.h decs.h makefile
	a2ps --media=letter --portrait -sduplex --columns=2 --rows=1 \
	--major=columns --file-align=virtual -o griffey.ps \
	$(SRCS) defs.h decs.h makefile

tar :    $(SRCS) defs.h decs.h makefile
	tar cvf griffey.tar $(SRCS) defs.h decs.h makefile ; \
	gzip --best griffey.tar
