include common.mk

#OPT += -DUSE_FGETS

target=sort_forests
SRC1=sort_forests.c utils.c progressbar.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = utils.h progressbar.h 

all: $(target) $(SRC1) $(INCL) Makefile 

$(target): $(OBJS1) $(INCL) Makefile common.mk
	$(CC) $(OBJS1) $(CLINK) $(CFLAGS) -o $@

%.o: %.c $(INCL) Makefile common.mk
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan 

clean:
	$(RM) $(target) $(OBJS1) 

clena: clean

celan: celan



