name = main
#compile_thread = -j8

#if you put the exteranl libaray in one directory

oplevel0 = -O0 -g  -DBOOST_LOG_DYN_LINK -fopenmp
oplevel3 = -O3 -ffast-math -DBOOST_LOG_DYN_LINK -fopenmp

FC  = mpiifort ${oplevel3}
CC  = mpiicc   ${oplevel3}
CXX = mpiicpc  ${oplevel3}

METIS_LIB = /THL8/home/zhangyi/software/metis-5.1.0/libmetis/

FCFLAGS  = -J${CURDIR}/ -fno-trapping-math -fno-signaling-nans -ffree-line-length-0 

CXXFLAGS = -fPIC --std=c++0x -Werror=return-type \
		-fno-trapping-math -fno-signaling-nans

CFLAGS   = -Werror=return-type -fno-trapping-math \
		-fno-signaling-nans

LDFLAGS = -L${METIS_LIB} 

LIBS =  -lstdc++ ${EXT_LIB} ${EXT_LIB_LINK} 

SRC_DIR = ../src/
TEST_DIR = ../test/
OBJ_DIR = ../obj/


OBJS =  $(addsuffix .o,  $(basename $(wildcard $(SRC_DIR)/*.cpp))) \
	$(addsuffix .o,  $(basename $(wildcard $(SRC_DIR)/*.F90)))

OBJS_TEST = $(addsuffix .o,  $(basename $(wildcard $(TEST_DIR)/*.cpp))) \
	   $(addsuffix .o,  $(basename $(wildcard $(TEST_DIR)/*.F90))) 


.DEFAULT_GOAL := all

MAKEFILE = makefile.linux


%.o: %.cpp
	$(CXX) ${EXT_INCLUDE} -c $(CXXFLAGS) $< -o $@

%.o: %.c
	$(CXX) ${EXT_INCLUDE} -c $(CFLAGS) $< -o $@

%.o: %.F90
	$(FC) ${FCFLAGS} -c $< -o $@


%.d: %.cpp
	@set -e; rm -f $@; $(CXX) ${EXT_INCLUDE} -MM $< $(INCLUDEFLAGS) > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.d: %.c
	@set -e; rm -f $@; $(CC) ${EXT_INCLUDE} -MM $< $(INCLUDEFLAGS) > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

%.d: %.F90 
	@set -e; rm -f $@; $(FC) ${FCFLAGS} -MM $< $(INCLUDEFLAGS) > $@.$$$$ 2>/dev/null; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$


-include $(OBJS:.o=.d) 
-include $(OBJS_TEST:.o=.d) 

debug:
	@echo ${OBJS}   

test : ${OBJS_TEST} lib
	-${CXX} -o main ${OBJS_TEST}  	${LDFLAGS} \
	-lgfortran -lmetis -lgrist_lib -lmpi_mpifh -L$(CURDIR)/ \
	@./main

lib: ${OBJS}
	@ar rcs libgrist_lib.a ${OBJS}

small:
	@make all -f ${MAKEFILE}
	@mpirun -n 4 ./main 4 3 2

clean:
	@find ../ -name "*.o" \
	-or -name "*.d*" \
	-or -name "*.log" \
	-or -name "*.mod" \
	-or -name "main" \
	-or -name "*.a" \
	-or -name "testall_main" | \
	xargs rm -f  &>/dev/null
