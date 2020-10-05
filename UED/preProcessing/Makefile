BASETOOLS=/reg/neh/home/khegazy/baseTools/tools
BASECODE=/reg/neh/home/khegazy/baseTools/UEDanalysis/preProcessing
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wno-unused-local-typedefs $(shell root-config --cflags) -I/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel7-gcc48-opt/include -I/reg/g/psdm/sw/external/opencv2/2.4.11/x86_64-rhel7-gcc48-opt/include -I/reg/g/psdm/sw/external/boost/1.55.0-python2.7/x86_64-rhel7-gcc48-opt/include -I/reg/neh/home/khegazy/packages/Eigen/include
LINKLIBS=$(shell root-config --libs) -L/reg/g/psdm/sw/external/fftw/3.3.4/x86_64-rhel7-gcc48-opt/lib -lfftw3 -L/reg/g/psdm/sw/external/opencv2/2.4.11/x86_64-rhel7-gcc48-opt/lib -lopencv_core -lopencv_highgui -lopencv_contrib -lopencv_imgproc -lopencv_objdetect -L/reg/g/psdm/sw/external/boost/1.55.0-python2.7/x86_64-rhel7-gcc48-opt/lib -lboost_system -lboost_math_c99 -lboost_math_tr1 

DEP_HPP_FILES := $(wildcard $(BASETOOLS)/*.h) 
DEP_CPP_FILES := $(wildcard $(BASETOOLS)/*.cpp)
LOCAL_DEP_CPP_FILES := ../../parameters.cpp preProcFuncts.cpp
DEPS := $(DEP_HPP_FILES) $(DEP_CPP_FILES)
LOCAL_DEPS := $(LOCAL_DEP_CPP_FILES)
CPP_FILES := $(wildcard $(BASECODE)/*.cpp)
HPP_FILES := $(wildcard $(BASECODE)/*.h)
OBJ_FILES := $(notdir $(CPP_FILES:%.cpp=%.o))
DEPS := $(DEP_HPP_FILES) $(DEP_CPP_FILES)
DEP_OBJ_FILES := $(DEP_CPP_FILES:%.cpp=%.o)
LOCAL_DEP_OBJ_FILES := $(LOCAL_DEP_CPP_FILES:%.cpp=%.o)
EXE_FILES := $(OBJ_FILES:%.o=%.exe)
STD := -std=c++0x

$(info ${OBJ_FILES})
$(info ${DEPS})

vpath %.cpp $(sort $(dir $(CPP_FILES)))

all: $(EXE_FILES)

.PRECIOUS: %.o

%.o: %.cpp $(LOCAL_DEPS) 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp $(DEPS) $(LOCAL_DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.exe: %.o $(OBJ_FILES) $(DEP_OBJ_FILES) $(LOCAL_DEP_OBJ_FILES)
	$(CXX) -o $@ $^ $(LINKLIBS)

clean:
	rm -f *.o *.exe ../../parameters.o ${BASETOOLS}/*.o
