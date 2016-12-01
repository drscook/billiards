.SUFFIXES: .cu .o
	.cu.o:
		$(GCC) $(COMPFLAGS) $*.cu
GPP = nvcc
COMPFLAGS = -std=c++11
LINKFLAGS = 
INCD = -I"/use/local/cuda/samples/common/inc"

3DOBJ = master.exe
LDFLAGS = -lGL -lGLU -lglut

$(3DOBJ) : main.cu
	$(GPP) -o $@ $< $(COMPFLAGS) $(INCD) $(LDFLAGS)
