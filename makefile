.SUFFIXES: .cu .o
	.cu.o:
		$(GCC) $(COMPFLAGS) $*.cu
GPP = nvcc
COMPFLAGS = -std=c++11
LINKFLAGS = 

3DOBJ = master.exe
LDFLAGS = -lGL -lGLU -lglut

$(3DOBJ) : main.cu
	$(GPP) -o $@ $< $(COMPFLAGS) $(LDFLAGS)
