

SOURCES=cgranges.c
OBJECTS = $(SOURCES:.c=.o)
LIBFILE = libcgranges.a

$(LIBFILE): $(OBJECTS)
	ar -rc $(LIBFILE)  $(OBJECTS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@


clean:
	${RM} *.o *.a
