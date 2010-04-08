#
# Mac OSX makefile for Perple_X 07
#
# To compile the Perple_X programs with this file first edit the compiler
# variables so that they are consistent with the fortran and C installation on 
# your system. Then type the command
#
#                             make <target>
#
# where target can be one or more of the programs you want compiled, e.g.,
# "make all" - make all the programs
# "make build vertex pssect" - make only these three programs
#
# To start over, you can type
#
#							make clean
#
# which will remove all programs and compiled code from the directory, leaving only the source code.
#
# NOTE: file conversion from windows to unix may result in a missing carraige
# return/line feed at the end of Perple_X source files and some compilers may 
# report this as an error, in which case edit the source file and add a blank line
# at the end of the file.
#
# JADC, Jan 19, 2008    
# 
##################### COMPILER VARIABLES ####################################  

# Unused options to consider:
# -static-libgfortran -arch x86_64 -O0 -g -pedantic -Wunused -arch ppc -arch i386  -arch x86_64 
# Notes: 
# -march=core2 = Intel Core2 CPU with 64-bit extensions, MMX, SSE, SSE2, SSE3 and SSSE3 instruction set support. 
# -mtune=generic = Produce code optimized for the most common IA32/AMD64/EM64T processors. If you know the CPU on which your code will run, then you should use the corresponding -mtune option instead of -mtune=generic. But, if you do not know exactly what CPU users of your application will have, then you should use this option. 

COMPILE_OPTIONS_32 = -m32 -march=i386 -mtune=core2 -msse2 -maccumulate-outgoing-args 

COMPILE_OPTIONS_64 = -m64 -march=core2 -mtune=core2 -msse2 -maccumulate-outgoing-args 

COMPILE_OPTIONS_PPC = -m32 -mcpu=G4

COMPILE_OPTIONS = -v -O3 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4  -ffast-math -funroll-loops --param max-unroll-times=4 -ftree-vectorize
                   
COMP77 = /sw/bin/gfortran 		# the name of the local fortran77-compatible compiler
COMP77_PPC = /ppc/bin/gfortran

export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/local/lib:/ppc/lib

BASE_FORTRAN_FLAGS = -pedantic -Wunused 	#desired additional compile options for fortran
FFLAGSPPC = $(BASE_FORTRAN_FLAGS) $(COMPILE_OPTIONS) $(COMPILE_OPTIONS_PPC) 
FFLAGS32 = $(BASE_FORTRAN_FLAGS) $(COMPILE_OPTIONS) $(COMPILE_OPTIONS_32) 
FFLAGS64 = $(BASE_FORTRAN_FLAGS) $(COMPILE_OPTIONS) $(COMPILE_OPTIONS_64) 

CC = gcc-4		# name of the local C compiler
CC_PPC = /ppc/bin/gcc	#location of the PPC C compiler

BASE_C_FLAGS = -pedantic -Wunused 		#desired additional compile options for C
CFLAGSPPC = $(BASE_C_FLAGS) $(COMPILE_OPTIONS) $(COMPILE_OPTIONS_PPC)
CFLAGS32 = $(BASE_C_FLAGS) $(COMPILE_OPTIONS) $(COMPILE_OPTIONS_32)
CFLAGS64 = $(BASE_C_FLAGS) $(COMPILE_OPTIONS) $(COMPILE_OPTIONS_64)

#LDFLAGS = $(COMPILE_OPTIONS_64)					#desired additional linker options


CAIRO_INCLUDE = /Users/dave/BuildCairo/include/cairo
CAIRO_LIB = /Users/dave/BuildCairo/lib
CAIRO_LIB_FLAGS = -L$(CAIRO_LIB) -lcairo
LIBRARYINSTALLATIONTARGET=/usr/local/lib/perplex

MYOBJ = actcor build cohsrk ctransf frendly hptover htog meemum pscontor pscontor_rel pspts psvdraw pssect pscontor_new pscontor_rel_new pspts_new psvdraw_new pssect_new pt2curv species vertex werami 

all:  $(MYOBJ)
	mkdir -p ./bin
	mv $(MYOBJ) ./bin

FORT_SOURCE = $(MYOBJ) tlib rlib flib clib dlib olib cont_lib nlib sclib vrlib resub

ALLPROG := $(foreach name,$(MYOBJ),$(name)_32) $(foreach name,$(MYOBJ),$(name)_64) $(foreach name,$(MYOBJ),$(name)_ppc)

clean: 
	rm -f *.o $(ALLPROG) $(MYOBJ)

cleanbin: 
	rm -f $(ALLPROG) $(MYOBJ)
	
test:
	echo $(ALLPROG)

FINISHLINK = lipo -create $@_32 $@_64 $@_ppc -output $@; install_name_tool -change /sw/lib/gcc4.4/lib/i386/libgfortran.3.dylib $(LIBRARYINSTALLATIONTARGET)/libgfortran.3.dylib $@; install_name_tool -change /sw/lib/gcc4.4/lib/libgcc_s.1.dylib $(LIBRARYINSTALLATIONTARGET)/libgcc_s.1.dylib $@;rm $@_32 $@_64 $@_ppc

FINISHLINKCAIRO = $(FINISHLINK);install_name_tool -change $(CAIRO_LIB)/libcairo.2.dylib $(LIBRARYINSTALLATIONTARGET)/libcairo.2.dylib $@


###################### TARGETS FOR FORTRAN PROGRAMS #########################   
# 
actcor: actcor_32 actcor_64 actcor_ppc
	$(FINISHLINK)
actcor_32: actcor_32.o tlib_32.o 
	$(COMP77) $(FFLAGS32) actcor_32.o tlib_32.o -o $@
actcor_64: actcor_64.o tlib_64.o 
	$(COMP77) $(FFLAGS64) actcor_64.o tlib_64.o -o $@
actcor_ppc: actcor_ppc.o tlib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) actcor_ppc.o tlib_ppc.o -o $@


build: build_32 build_64 build_ppc
	$(FINISHLINK)
build_32: build_32.o tlib_32.o rlib_32.o flib_32.o 
	$(COMP77) $(FFLAGS32) build_32.o tlib_32.o rlib_32.o flib_32.o -o $@
build_64: build_64.o tlib_64.o rlib_64.o flib_64.o 
	$(COMP77) $(FFLAGS64) build_64.o tlib_64.o rlib_64.o flib_64.o -o $@
build_ppc: build_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o 
	$(COMP77_PPC) $(FFLAGSPPC) build_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o -o $@

cohsrk: cohsrk_32 cohsrk_64 cohsrk_ppc
	$(FINISHLINK)
cohsrk_32: cohsrk_32.o tlib_32.o flib_32.o 
	$(COMP77) $(FFLAGS32) cohsrk_32.o tlib_32.o flib_32.o -o $@
cohsrk_64: cohsrk_64.o tlib_64.o flib_64.o 
	$(COMP77) $(FFLAGS64) cohsrk_64.o tlib_64.o flib_64.o -o $@
cohsrk_ppc: cohsrk_ppc.o tlib_ppc.o flib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) cohsrk_ppc.o tlib_ppc.o flib_ppc.o -o $@


ctransf: ctransf_32 ctransf_64 ctransf_ppc
	$(FINISHLINK)
ctransf_32: ctransf_32.o tlib_32.o 
	$(COMP77) $(FFLAGS32) ctransf_32.o tlib_32.o -o $@
ctransf_64: ctransf_64.o tlib_64.o 
	$(COMP77) $(FFLAGS64) ctransf_64.o tlib_64.o -o $@
ctransf_ppc: ctransf_ppc.o tlib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) ctransf_ppc.o tlib_ppc.o -o $@


frendly: frendly_32 frendly_64 frendly_ppc
	$(FINISHLINK)
frendly_32: frendly_32.o tlib_32.o rlib_32.o flib_32.o 
	$(COMP77) $(FFLAGS32) frendly_32.o tlib_32.o rlib_32.o flib_32.o -o $@
frendly_64: frendly_64.o tlib_64.o rlib_64.o flib_64.o 
	$(COMP77) $(FFLAGS64) frendly_64.o tlib_64.o rlib_64.o flib_64.o -o $@
frendly_ppc: frendly_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) frendly_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o -o $@


hptover: hptover_32 hptover_64 hptover_ppc
	$(FINISHLINK)
hptover_32: hptover_32.o
	$(COMP77) $(FFLAGS32) hptover_32.o -o $@
hptover_64: hptover_64.o
	$(COMP77) $(FFLAGS64) hptover_64.o -o $@
hptover_ppc: hptover_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) hptover_ppc.o -o $@


htog: htog_32 htog_64 htog_ppc
	$(FINISHLINK)
htog_32: htog_32.o
	$(COMP77) $(FFLAGS32) htog_32.o -o $@
htog_64: htog_64.o
	$(COMP77) $(FFLAGS64) htog_64.o -o $@
htog_ppc: htog_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) htog_ppc.o -o $@


meemum: meemum_32 meemum_64 meemum_ppc
	$(FINISHLINK)
meemum_32: meemum_32.o tlib_32.o rlib_32.o flib_32.o nlib_32.o clib_32.o resub_32.o olib_32.o
	$(COMP77) $(FFLAGS32) meemum_32.o tlib_32.o rlib_32.o flib_32.o nlib_32.o clib_32.o resub_32.o olib_32.o -o $@
meemum_64: meemum_64.o tlib_64.o rlib_64.o flib_64.o nlib_64.o clib_64.o resub_64.o olib_64.o
	$(COMP77) $(FFLAGS64) meemum_64.o tlib_64.o rlib_64.o flib_64.o nlib_64.o clib_64.o resub_64.o olib_64.o -o $@
meemum_ppc: meemum_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o nlib_ppc.o clib_ppc.o resub_ppc.o olib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) meemum_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o nlib_ppc.o clib_ppc.o resub_ppc.o olib_ppc.o -o $@
   

pscontor: pscontor_32 pscontor_64 pscontor_ppc
	$(FINISHLINK)
pscontor_32: pscontor_32.o pscom_32.o pslib_32.o tlib_32.o cont_lib_32.o 
	$(COMP77) $(FFLAGS32) pscontor_32.o pscom_32.o pslib_32.o tlib_32.o cont_lib_32.o -o $@
pscontor_64: pscontor_64.o pscom_64.o pslib_64.o tlib_64.o cont_lib_64.o 
	$(COMP77) $(FFLAGS64) pscontor_64.o pscom_64.o pslib_64.o tlib_64.o cont_lib_64.o -o $@
pscontor_ppc: pscontor_ppc.o pscom_ppc.o pslib_ppc.o tlib_ppc.o cont_lib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) pscontor_ppc.o pscom_ppc.o pslib_ppc.o tlib_ppc.o cont_lib_ppc.o -o $@

pscontor_rel: pscontor_rel_32 pscontor_rel_64 pscontor_rel_ppc
	$(FINISHLINK)
pscontor_rel_32: pscontor_rel_32.o pslib_32.o pscom_32.o tlib_32.o cont_lib_32.o
	$(COMP77) $(FFLAGS32) pscontor_rel_32.o pslib_32.o pscom_32.o tlib_32.o cont_lib_32.o -o $@
pscontor_rel_64: pscontor_rel_64.o pslib_64.o pscom_64.o tlib_64.o cont_lib_64.o
	$(COMP77) $(FFLAGS64) pscontor_rel_64.o pslib_64.o pscom_64.o tlib_64.o cont_lib_64.o -o $@
pscontor_rel_ppc: pscontor_rel_ppc.o pslib_ppc.o pscom_ppc.o tlib_ppc.o cont_lib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) pscontor_rel_ppc.o pslib_ppc.o pscom_ppc.o tlib_ppc.o cont_lib_ppc.o -o $@

pspts: pspts_32 pspts_64 pspts_ppc
	$(FINISHLINK)
pspts_32: pspts_32.o pslib_32.o tlib_32.o pscom_32.o 
	$(COMP77) $(FFLAGS32) pspts_32.o pslib_32.o tlib_32.o pscom_32.o -o $@
pspts_64: pspts_64.o pslib_64.o tlib_64.o pscom_64.o 
	$(COMP77) $(FFLAGS64) pspts_64.o pslib_64.o tlib_64.o pscom_64.o -o $@
pspts_ppc: pspts_ppc.o pslib_ppc.o tlib_ppc.o pscom_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) pspts_ppc.o pslib_ppc.o tlib_ppc.o pscom_ppc.o -o $@

psvdraw: psvdraw_32 psvdraw_64 psvdraw_ppc
	$(FINISHLINK)
psvdraw_32: psvdraw_32.o pslib_32.o tlib_32.o pscom_32.o
	$(COMP77) $(FFLAGS32) psvdraw_32.o pslib_32.o tlib_32.o pscom_32.o -o $@
psvdraw_64: psvdraw_64.o pslib_64.o tlib_64.o pscom_64.o
	$(COMP77) $(FFLAGS64) psvdraw_64.o pslib_64.o tlib_64.o pscom_64.o -o $@
psvdraw_ppc: psvdraw_ppc.o pslib_ppc.o tlib_ppc.o pscom_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) psvdraw_ppc.o pslib_ppc.o tlib_ppc.o pscom_ppc.o -o $@

pssect: pssect_32 pssect_64 pssect_ppc
	$(FINISHLINK)
pssect_32: pssect_32.o pslib_32.o pscom_32.o tlib_32.o rlib_32.o flib_32.o clib_32.o dlib_32.o 
	$(COMP77) $(FFLAGS32) pssect_32.o pslib_32.o pscom_32.o tlib_32.o rlib_32.o flib_32.o clib_32.o dlib_32.o -o $@
pssect_64: pssect_64.o pslib_64.o pscom_64.o tlib_64.o rlib_64.o flib_64.o clib_64.o dlib_64.o 
	$(COMP77) $(FFLAGS64) pssect_64.o pslib_64.o pscom_64.o tlib_64.o rlib_64.o flib_64.o clib_64.o dlib_64.o -o $@
pssect_ppc: pssect_ppc.o pslib_ppc.o pscom_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o clib_ppc.o dlib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) pssect_ppc.o pslib_ppc.o pscom_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o clib_ppc.o dlib_ppc.o -o $@

## NEW VERSIONS ##
pscontor_new: pscontor_new_32 pscontor_new_64 pscontor_new_ppc
	$(FINISHLINKCAIRO)
pscontor_new_32: pscontor_new_32.o pscom_new_32.o graphic_output_32.o tlib_32.o cont_lib_32.o 
	$(COMP77) $(FFLAGS32) $(CAIRO_LIB_FLAGS) pscontor_new_32.o pscom_new_32.o graphic_output_32.o tlib_32.o cont_lib_32.o -o $@
pscontor_new_64: pscontor_new_64.o pscom_new_64.o graphic_output_64.o tlib_64.o cont_lib_64.o 
	$(COMP77) $(FFLAGS64) $(CAIRO_LIB_FLAGS) pscontor_new_64.o pscom_new_64.o graphic_output_64.o tlib_64.o cont_lib_64.o -o $@
pscontor_new_ppc: pscontor_new_ppc.o pscom_new_ppc.o graphic_output_ppc.o tlib_ppc.o cont_lib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) $(CAIRO_LIB_FLAGS) pscontor_new_ppc.o pscom_new_ppc.o graphic_output_ppc.o tlib_ppc.o cont_lib_ppc.o -o $@

pscontor_rel_new: pscontor_rel_new_32 pscontor_rel_new_64 pscontor_rel_new_ppc
	$(FINISHLINKCAIRO)
pscontor_rel_new_32: pscontor_rel_new_32.o graphic_output_32.o pscom_new_32.o tlib_32.o cont_lib_32.o
	$(COMP77) $(FFLAGS32) $(CAIRO_LIB_FLAGS) pscontor_rel_new_32.o graphic_output_32.o pscom_new_32.o tlib_32.o cont_lib_32.o -o $@
pscontor_rel_new_64: pscontor_rel_new_64.o graphic_output_64.o pscom_new_64.o tlib_64.o cont_lib_64.o
	$(COMP77) $(FFLAGS64) $(CAIRO_LIB_FLAGS) pscontor_rel_new_64.o graphic_output_64.o pscom_new_64.o tlib_64.o cont_lib_64.o -o $@
pscontor_rel_new_ppc: pscontor_rel_new_ppc.o graphic_output_ppc.o pscom_new_ppc.o tlib_ppc.o cont_lib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) $(CAIRO_LIB_FLAGS) pscontor_rel_new_ppc.o graphic_output_ppc.o pscom_new_ppc.o tlib_ppc.o cont_lib_ppc.o -o $@

pspts_new: pspts_new_32 pspts_new_64 pspts_new_ppc
	$(FINISHLINKCAIRO)
pspts_new_32: pspts_new_32.o graphic_output_32.o tlib_32.o pscom_new_32.o 
	$(COMP77) $(FFLAGS32) $(CAIRO_LIB_FLAGS) pspts_new_32.o graphic_output_32.o tlib_32.o pscom_new_32.o -o $@
pspts_new_64: pspts_new_64.o graphic_output_64.o tlib_64.o pscom_new_64.o 
	$(COMP77) $(FFLAGS64) $(CAIRO_LIB_FLAGS) pspts_new_64.o graphic_output_64.o tlib_64.o pscom_new_64.o -o $@
pspts_new_ppc: pspts_new_ppc.o graphic_output_ppc.o tlib_ppc.o pscom_new_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) $(CAIRO_LIB_FLAGS) pspts_new_ppc.o graphic_output_ppc.o tlib_ppc.o pscom_new_ppc.o -o $@

psvdraw_new: psvdraw_new_32 psvdraw_new_64 psvdraw_new_ppc
	$(FINISHLINKCAIRO)
psvdraw_new_32: psvdraw_new_32.o graphic_output_32.o tlib_32.o pscom_new_32.o
	$(COMP77) $(FFLAGS32) $(CAIRO_LIB_FLAGS) psvdraw_new_32.o graphic_output_32.o tlib_32.o pscom_new_32.o -o $@
psvdraw_new_64: psvdraw_new_64.o graphic_output_64.o tlib_64.o pscom_new_64.o
	$(COMP77) $(FFLAGS64) $(CAIRO_LIB_FLAGS) psvdraw_new_64.o graphic_output_64.o tlib_64.o pscom_new_64.o -o $@
psvdraw_new_ppc: psvdraw_new_ppc.o graphic_output_ppc.o tlib_ppc.o pscom_new_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) $(CAIRO_LIB_FLAGS) psvdraw_new_ppc.o graphic_output_ppc.o tlib_ppc.o pscom_new_ppc.o -o $@

pssect_new: pssect_new_32 pssect_new_64 pssect_new_ppc
	$(FINISHLINKCAIRO)
pssect_new_32: pssect_new_32.o graphic_output_32.o pscom_new_32.o tlib_32.o rlib_32.o flib_32.o clib_32.o dlib_32.o 
	$(COMP77) $(FFLAGS32) $(CAIRO_LIB_FLAGS) pssect_new_32.o graphic_output_32.o pscom_new_32.o tlib_32.o rlib_32.o flib_32.o clib_32.o dlib_32.o -o $@
pssect_new_64: pssect_new_64.o graphic_output_64.o pscom_new_64.o tlib_64.o rlib_64.o flib_64.o clib_64.o dlib_64.o 
	$(COMP77) $(FFLAGS64) $(CAIRO_LIB_FLAGS) pssect_new_64.o graphic_output_64.o pscom_new_64.o tlib_64.o rlib_64.o flib_64.o clib_64.o dlib_64.o -o $@
pssect_new_ppc: pssect_new_ppc.o graphic_output_ppc.o pscom_new_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o clib_ppc.o dlib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) $(CAIRO_LIB_FLAGS) pssect_new_ppc.o graphic_output_ppc.o pscom_new_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o clib_ppc.o dlib_ppc.o -o $@


pt2curv: pt2curv_32 pt2curv_64 pt2curv_ppc
	$(FINISHLINK)
pt2curv_32: pt2curv_32.o
	$(COMP77) $(FFLAGS32) pt2curv_32.o -o $@
pt2curv_64: pt2curv_64.o
	$(COMP77) $(FFLAGS64) pt2curv_64.o -o $@
pt2curv_ppc: pt2curv_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) pt2curv_ppc.o -o $@


species: species_32 species_64 species_ppc
	$(FINISHLINK)
species_32: species_32.o flib_32.o tlib_32.o
	$(COMP77) $(FFLAGS32) species_32.o flib_32.o tlib_32.o -o $@
species_64: species_64.o flib_64.o tlib_64.o
	$(COMP77) $(FFLAGS64) species_64.o flib_64.o tlib_64.o -o $@
species_ppc: species_ppc.o flib_ppc.o tlib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) species_ppc.o flib_ppc.o tlib_ppc.o -o $@

         
vertex: vertex_32 vertex_64 vertex_ppc
	$(FINISHLINK)
vertex_32: vertex_32.o tlib_32.o rlib_32.o flib_32.o nlib_32.o clib_32.o resub_32.o
	$(COMP77) $(FFLAGS32) vertex_32.o tlib_32.o rlib_32.o flib_32.o nlib_32.o clib_32.o resub_32.o -o $@
vertex_64: vertex_64.o tlib_64.o rlib_64.o flib_64.o nlib_64.o clib_64.o resub_64.o
	$(COMP77) $(FFLAGS64) vertex_64.o tlib_64.o rlib_64.o flib_64.o nlib_64.o clib_64.o resub_64.o -o $@
vertex_ppc: vertex_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o nlib_ppc.o clib_ppc.o resub_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) vertex_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o nlib_ppc.o clib_ppc.o resub_ppc.o -o $@

werami: werami_32 werami_64 werami_ppc
	$(FINISHLINK)
werami_32: reader_32.o tlib_32.o rlib_32.o flib_32.o clib_32.o dlib_32.o olib_32.o
	$(COMP77) $(FFLAGS32) reader_32.o tlib_32.o rlib_32.o flib_32.o clib_32.o dlib_32.o olib_32.o -o $@
werami_64: reader_64.o tlib_64.o rlib_64.o flib_64.o clib_64.o dlib_64.o olib_64.o
	$(COMP77) $(FFLAGS64) reader_64.o tlib_64.o rlib_64.o flib_64.o clib_64.o dlib_64.o olib_64.o -o $@
werami_ppc: reader_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o clib_ppc.o dlib_ppc.o olib_ppc.o
	$(COMP77_PPC) $(FFLAGSPPC) reader_ppc.o tlib_ppc.o rlib_ppc.o flib_ppc.o clib_ppc.o dlib_ppc.o olib_ppc.o -o $@


# targets missing from '07:
#rk: rk.o flib.o tlib.o
#	$(COMP77) $(FFLAGS) $@.o tlib.o flib.o -o $@
#ps_p_contor: ps_p_contor.o pslib.o
#	$(COMP77) $(FFLAGS) $@.o pslib.o -o $@
#ge0pt: ge0pt.o
#	$(COMP77) $(FFLAGS) $@.o -o $@
#satsurf: satsurf.o flib.o tlib.o
#	$(COMP77) $(FFLAGS) $@.o flib.o tlib.o -o $@
#rewrite,sox,gox,cohscont 

#################################################################################

actcor_32.o: actcor.f
	$(COMP77) $(FFLAGS32) -c actcor.f -o $@
build_32.o: build.f
	$(COMP77) $(FFLAGS32) -c build.f -o $@
cohsrk_32.o: cohsrk.f
	$(COMP77) $(FFLAGS32) -c cohsrk.f -o $@
cont_lib_32.o: cont_lib.f
	$(COMP77) $(FFLAGS32) -c cont_lib.f -o $@
ctransf_32.o: ctransf.f
	$(COMP77) $(FFLAGS32) -c ctransf.f -o $@
frendly_32.o: frendly.f
	$(COMP77) $(FFLAGS32) -c frendly.f -o $@
ge0pt_32.o: ge0pt.f
	$(COMP77) $(FFLAGS32) -c ge0pt.f -o $@
hptover_32.o: hptover.f
	$(COMP77) $(FFLAGS32) -c hptover.f -o $@
htog_32.o: htog.f
	$(COMP77) $(FFLAGS32) -c htog.f -o $@
pssect_32.o: pssect.f
	$(COMP77) $(FFLAGS32) -c pssect.f -o $@
psvdraw_32.o: psvdraw.f
	$(COMP77) $(FFLAGS32) -c psvdraw.f -o $@
pscom_32.o: pscom.f
	$(COMP77) $(FFLAGS32) -c pscom.f -o $@
pspts_32.o: pspts.f
	$(COMP77) $(FFLAGS32) -c pspts.f -o $@
pscontor_32.o: pscontor.f
	$(COMP77) $(FFLAGS32) -c pscontor.f -o $@
pscontor_rel_32.o: pscontor_rel.f
	$(COMP77) $(FFLAGS32) -c pscontor_rel.f -o $@

### NEW VERSIONS:
psvdraw_new_32.o: psvdraw_new.f
	$(COMP77) $(FFLAGS32) -c psvdraw_new.f -o $@
pscom_new_32.o: pscom_new.f
	$(COMP77) $(FFLAGS32) -c pscom_new.f -o $@
pssect_new_32.o: pssect_new.f
	$(COMP77) $(FFLAGS32) -c pssect_new.f -o $@
pspts_new_32.o: pspts_new.f
	$(COMP77) $(FFLAGS32) -c pspts_new.f -o $@
pscontor_new_32.o: pscontor_new.f
	$(COMP77) $(FFLAGS32) -c pscontor_new.f -o $@
pscontor_rel_new_32.o: pscontor_rel_new.f
	$(COMP77) $(FFLAGS32) -c pscontor_rel_new.f -o $@

pt2curv_32.o: pt2curv.f
	$(COMP77) $(FFLAGS32) -c pt2curv.f -o $@
satsurf_32.o: satsurf.f
	$(COMP77) $(FFLAGS32) -c satsurf.f -o $@
species_32.o: species.f
	$(COMP77) $(FFLAGS32) -c species.f -o $@
rk_32.o: rk.f
	$(COMP77) $(FFLAGS32) -c rk.f -o $@
reader_32.o: reader.f
	$(COMP77) $(FFLAGS32) -c reader.f -o $@
resub_32.o: resub.f
	$(COMP77) $(FFLAGS32) -c resub.f -o $@
meemum_32.o: meemum.f
	$(COMP77) $(FFLAGS32) -c meemum.f -o $@
vertex_32.o: vertex.f
	$(COMP77) $(FFLAGS32) -c vertex.f -o $@
clib_32.o: clib.f
	$(COMP77) $(FFLAGS32) -c clib.f -o $@
dlib_32.o: dlib.f
	$(COMP77) $(FFLAGS32) -c dlib.f -o $@
flib_32.o: flib.f
	$(COMP77) $(FFLAGS32) -c flib.f -o $@
pslib_32.o: pslib.f
	$(COMP77) $(FFLAGS32) -c pslib.f -o $@
nlib_32.o: nlib.f
	$(COMP77) $(FFLAGS32) -c nlib.f -o $@
rlib_32.o: rlib.f
	$(COMP77) $(FFLAGS32) -c rlib.f -o $@
tlib_32.o: tlib.f
	$(COMP77) $(FFLAGS32) -c tlib.f -o $@
olib_32.o: olib.f
	$(COMP77) $(FFLAGS32) -c olib.f -o $@
biglib_32.o: biglib.f
	$(COMP77) $(FFLAGS32) -c biglib.f -o $@

# DMH 2010Feb02
graphic_output_32.o: graphic_output.c
	$(CC) $(CFLAGS32) $(INCL) -c -I$(CAIRO_INCLUDE) graphic_output.c -o $@


actcor_64.o: actcor.f
	$(COMP77) $(FFLAGS64) -c actcor.f -o $@
build_64.o: build.f
	$(COMP77) $(FFLAGS64) -c build.f -o $@
cohsrk_64.o: cohsrk.f
	$(COMP77) $(FFLAGS64) -c cohsrk.f -o $@
cont_lib_64.o: cont_lib.f
	$(COMP77) $(FFLAGS64) -c cont_lib.f -o $@
ctransf_64.o: ctransf.f
	$(COMP77) $(FFLAGS64) -c ctransf.f -o $@
frendly_64.o: frendly.f
	$(COMP77) $(FFLAGS64) -c frendly.f -o $@
ge0pt_64.o: ge0pt.f
	$(COMP77) $(FFLAGS64) -c ge0pt.f -o $@
hptover_64.o: hptover.f
	$(COMP77) $(FFLAGS64) -c hptover.f -o $@
htog_64.o: htog.f
	$(COMP77) $(FFLAGS64) -c htog.f -o $@
pssect_64.o: pssect.f
	$(COMP77) $(FFLAGS64) -c pssect.f -o $@
psvdraw_64.o: psvdraw.f
	$(COMP77) $(FFLAGS64) -c psvdraw.f -o $@
pscom_64.o: pscom.f
	$(COMP77) $(FFLAGS64) -c pscom.f -o $@
pspts_64.o: pspts.f
	$(COMP77) $(FFLAGS64) -c pspts.f -o $@
pscontor_64.o: pscontor.f
	$(COMP77) $(FFLAGS64) -c pscontor.f -o $@
pscontor_rel_64.o: pscontor_rel.f
	$(COMP77) $(FFLAGS64) -c pscontor_rel.f -o $@

### NEW VERSIONS:
psvdraw_new_64.o: psvdraw_new.f
	$(COMP77) $(FFLAGS64) -c psvdraw_new.f -o $@
pscom_new_64.o: pscom_new.f
	$(COMP77) $(FFLAGS64) -c pscom_new.f -o $@
pssect_new_64.o: pssect_new.f
	$(COMP77) $(FFLAGS64) -c pssect_new.f -o $@
pspts_new_64.o: pspts_new.f
	$(COMP77) $(FFLAGS64) -c pspts_new.f -o $@
pscontor_new_64.o: pscontor_new.f
	$(COMP77) $(FFLAGS64) -c pscontor_new.f -o $@
pscontor_rel_new_64.o: pscontor_rel_new.f
	$(COMP77) $(FFLAGS64) -c pscontor_rel_new.f -o $@

pt2curv_64.o: pt2curv.f
	$(COMP77) $(FFLAGS64) -c pt2curv.f -o $@
satsurf_64.o: satsurf.f
	$(COMP77) $(FFLAGS64) -c satsurf.f -o $@
species_64.o: species.f
	$(COMP77) $(FFLAGS64) -c species.f -o $@
rk_64.o: rk.f
	$(COMP77) $(FFLAGS64) -c rk.f -o $@
reader_64.o: reader.f
	$(COMP77) $(FFLAGS64) -c reader.f -o $@
resub_64.o: resub.f
	$(COMP77) $(FFLAGS64) -c resub.f -o $@
meemum_64.o: meemum.f
	$(COMP77) $(FFLAGS64) -c meemum.f -o $@
vertex_64.o: vertex.f
	$(COMP77) $(FFLAGS64) -c vertex.f -o $@
clib_64.o: clib.f
	$(COMP77) $(FFLAGS64) -c clib.f -o $@
dlib_64.o: dlib.f
	$(COMP77) $(FFLAGS64) -c dlib.f -o $@
flib_64.o: flib.f
	$(COMP77) $(FFLAGS64) -c flib.f -o $@
pslib_64.o: pslib.f
	$(COMP77) $(FFLAGS64) -c pslib.f -o $@
nlib_64.o: nlib.f
	$(COMP77) $(FFLAGS64) -c nlib.f -o $@
rlib_64.o: rlib.f
	$(COMP77) $(FFLAGS64) -c rlib.f -o $@
tlib_64.o: tlib.f
	$(COMP77) $(FFLAGS64) -c tlib.f -o $@
olib_64.o: olib.f
	$(COMP77) $(FFLAGS64) -c olib.f -o $@
biglib_64.o: biglib.f
	$(COMP77) $(FFLAGS64) -c biglib.f -o $@

# DMH 2010Feb02
graphic_output_64.o: graphic_output.c
	$(CC) $(CFLAGS64) $(INCL) -c -I$(CAIRO_INCLUDE) graphic_output.c -o $@
	


actcor_ppc.o: actcor.f
	$(COMP77_PPC) $(FFLAGSPPC) -c actcor.f -o $@
build_ppc.o: build.f
	$(COMP77_PPC) $(FFLAGSPPC) -c build.f -o $@
cohsrk_ppc.o: cohsrk.f
	$(COMP77_PPC) $(FFLAGSPPC) -c cohsrk.f -o $@
cont_lib_ppc.o: cont_lib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c cont_lib.f -o $@
ctransf_ppc.o: ctransf.f
	$(COMP77_PPC) $(FFLAGSPPC) -c ctransf.f -o $@
frendly_ppc.o: frendly.f
	$(COMP77_PPC) $(FFLAGSPPC) -c frendly.f -o $@
ge0pt_ppc.o: ge0pt.f
	$(COMP77_PPC) $(FFLAGSPPC) -c ge0pt.f -o $@
hptover_ppc.o: hptover.f
	$(COMP77_PPC) $(FFLAGSPPC) -c hptover.f -o $@
htog_ppc.o: htog.f
	$(COMP77_PPC) $(FFLAGSPPC) -c htog.f -o $@
pssect_ppc.o: pssect.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pssect.f -o $@
psvdraw_ppc.o: psvdraw.f
	$(COMP77_PPC) $(FFLAGSPPC) -c psvdraw.f -o $@
pscom_ppc.o: pscom.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pscom.f -o $@
pspts_ppc.o: pspts.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pspts.f -o $@
pscontor_ppc.o: pscontor.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pscontor.f -o $@
pscontor_rel_ppc.o: pscontor_rel.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pscontor_rel.f -o $@

### NEW VERSIONS:
psvdraw_new_ppc.o: psvdraw_new.f
	$(COMP77_PPC) $(FFLAGSPPC) -c psvdraw_new.f -o $@
pscom_new_ppc.o: pscom_new.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pscom_new.f -o $@
pssect_new_ppc.o: pssect_new.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pssect_new.f -o $@
pspts_new_ppc.o: pspts_new.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pspts_new.f -o $@
pscontor_new_ppc.o: pscontor_new.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pscontor_new.f -o $@
pscontor_rel_new_ppc.o: pscontor_rel_new.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pscontor_rel_new.f -o $@

pt2curv_ppc.o: pt2curv.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pt2curv.f -o $@
satsurf_ppc.o: satsurf.f
	$(COMP77_PPC) $(FFLAGSPPC) -c satsurf.f -o $@
species_ppc.o: species.f
	$(COMP77_PPC) $(FFLAGSPPC) -c species.f -o $@
rk_ppc.o: rk.f
	$(COMP77_PPC) $(FFLAGSPPC) -c rk.f -o $@
reader_ppc.o: reader.f
	$(COMP77_PPC) $(FFLAGSPPC) -c reader.f -o $@
resub_ppc.o: resub.f
	$(COMP77_PPC) $(FFLAGSPPC) -c resub.f -o $@
meemum_ppc.o: meemum.f
	$(COMP77_PPC) $(FFLAGSPPC) -c meemum.f -o $@
vertex_ppc.o: vertex.f
	$(COMP77_PPC) $(FFLAGSPPC) -c vertex.f -o $@
clib_ppc.o: clib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c clib.f -o $@
dlib_ppc.o: dlib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c dlib.f -o $@
flib_ppc.o: flib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c flib.f -o $@
pslib_ppc.o: pslib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c pslib.f -o $@
nlib_ppc.o: nlib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c nlib.f -o $@
rlib_ppc.o: rlib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c rlib.f -o $@
tlib_ppc.o: tlib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c tlib.f -o $@
olib_ppc.o: olib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c olib.f -o $@
biglib_ppc.o: biglib.f
	$(COMP77_PPC) $(FFLAGSPPC) -c biglib.f -o $@

# DMH 2010Feb02
graphic_output_ppc.o: graphic_output.c
	$(CC_PPC) $(CFLAGSPPC) $(INCL) -c -I$(CAIRO_INCLUDE) graphic_output.c -o $@