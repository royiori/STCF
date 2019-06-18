# Add ROOT
ROOTETC = $(shell root-config --etcdir)/Makefile.arch
include $(ROOTETC)

# Compiler flags
CXXFLAGS  += -Wall
CXXFLAGS  += -I. -I./inc

LDFLAGS += -lGui

#------------------------------------------------------------------------------

SRCS = $(wildcard src/*.cpp)
HHDS = $(SRCS:src/%.cpp=inc/%.h)
OBJS = $(SRCS:src/%.cpp=obj/%.o)

LibSrc = lib/RICHDict.cpp
LibObj = obj/RICHDict.o
LibOut = lib/libRICHDict.so

#RICHS = $(wildcard *.cpp)
#RICHO = $(patsubst %.cpp, %.o,$(notdir ${RICHS}))

#------------------------------------------------------------------------------
MAIN = RICH
all: $(MAIN)


$(OBJS) : obj/%.o:src/%.cpp inc/%.h

$(LibSrc): ./inc/MyRICHDetector.h ./inc/LinkDef.h
		@echo "Generating dictionary $@..."
		$(ROOTCLING) -f $@ -c $^

$(LibOut): $(LibObj)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)


$(MAIN): %:%.o ${OBJS} ${LibObj}
		$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)
		@cp ./lib/RICHDict_rdict.pcm ./
		@echo "----> $@ done!"

#RICH: $(RICHO) 
#		$(LD) $(LDFLAGS) $^ $(GLIBS) -lEG -lHtml  $(OutPutOpt)$@
#		@echo "----> $@ done!"

clean:
		@rm -f $(OBJS)
		@rm -f RICHDict*
		@rm -f ./src/RICHDict*
		@rm -f *.o 
#
#%.o:%.C
#`		$(CXX) $(CXXFLAGS) -c $< -o $@

$(MAIN).o : $(MAIN).cpp $(HHDS)
	$(CXX) -o $@ $(CXXFLAGS) -c $<

obj/%.o : src/%.cpp
		$(CXX) -o $@ $(CXXFLAGS) -c $<

obj/%.o : lib/%.cpp
		$(CXX) -o $@ $(CXXFLAGS) -c $<