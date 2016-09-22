.SUFFIXES: .cpp .c .o
.cpp.o:
	$(GCC) $(COMPFLAGS) $*.cpp
.c.o:
	$(GCC) $(COMPFLAGS) $*.c

GPP = g++
COMPFLAGS = -std=gnu++11 
LINKFLAGS = 

# assume it's linux
3DOBJ = serial.exe
LDFLAGS = -lGL -lGLU -lglut -lm

# check if it's a windows machine
ifeq "$(OS)" "Windows_NT"
	LDFLAGS = -lopeng132 -lglu32 -lglut32

#if it's not windows, check if it's mac 
else
	UNAME_S := $(shell uname -s)
	ifeq "$(UNAME_S)" "Darwin"
		LDFLAGS = -framework OpenGL -framework GLUT
	endif
endif

$(3DOBJ) : serial.cpp
	$(GPP) -o $@ $< $(COMPFLAGS) $(LDFLAGS) 
