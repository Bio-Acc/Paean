CXXFLAGS= -O3 -fPIC -std=c++11

#HIP_PATH=`hipconfig --path`
#HIP_PATH=/public/software/compiler/rocm/rocm-3.3.0
HIP_PATH=/opt/rocm

#HIPCC= /opt/rocm/bin/hipcc
#HIPCC= ${HIP_PATH}/bin/hipcc
HIPCC= OMPI_CXX=hipcc mpic++

#rocminclude=-I/opt/rocm/include/hip/hip/
rocminclude=-I${HIP_PATH}/include/hip/hip

#hiplib=-L/opt/rocm/hip/lib -lhip_hcc -lhiprtc
hiplib=-L${HIP_PATH}/hip/lib -lhip_hcc -lhiprtc

#hipcubinc=-I/opt/rocm/hipcub/include/
#hipcublib=-L/opt/rocm/hipcub/
hipcubinc=-I${HIP_PATH}/hipcub/include
hipcublib=-L${HIP_PATH}/hipcub

gff_source= $(wildcard ./src/gff/*.cpp)
fmt_source= $(wildcard ./src/fmt/*.cc)

htsinclude=-I/public/home/zhoucb/wangsq/htslib-develop/include
htslib=-L/public/home/zhoucb/wangsq/htslib-develop -lhts
htslibpath=/public/home/zhoucb/wangsq/htslib-develop/libhts.so.3

projectinc=-I./include

OBJS= src/libgff.so src/libfmt.so src/libparse.so src/libfusion.so 
GOBJS= src/libcub_sort.o
CPULIB= src/cpu.a
GPULIB = src/gpu.a
EXE = src/paean

TOP_DIR=.
INCLUDE_DIRS= -I$(TOP_DIR)/include 

ALL:$(CPULIB) $(GPULIB) $(EXE)
$(CPULIB):$(OBJS)
$(GPULIB):$(GOBJS)
src/libgff.so:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) $(gff_source) -shared -o src/libgff.so
src/libparse.so:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) $(htsinclude) src/parse.cpp -shared -o src/libparse.so
src/libfusion.so:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) src/fusion.cpp -shared -o src/libfusion.so
src/libfmt.so:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) $(fmt_source) -shared -o src/libfmt.so
src/libcub_sort.o:
	$(HIPCC) -c $(rocminclude) $(projectinc) $(CXXFLAGS) \
	src/cub_sort.cu
	mv cub_sort.o src/
src/paean:
	$(HIPCC) -c src/gene.cu $(rocminclude) $(projectinc) $(CXXFLAGS) 
	mv gene.o src/
	cp $(htslibpath) src/
	$(HIPCC) -o src/paean src/gene.o src/cub_sort.o $(rocminclude) $(projectinc) $(CXXFLAGS) -L./src -lz -lm -pthread -lgff -lfmt -lparse -lfusion $(htslib) $(hiplib)
	

.PHONY :clean
clean:
	-rm  src/*.o
	-rm  src/*.so
	-rm  src/paean

