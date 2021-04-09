# Add ROOT
ROOTETC = $(shell root-config --etcdir)/Makefile.arch
include $(ROOTETC)

# Compiler flags
CXXFLAGS  += -Wall
CXXFLAGS  += -I. -I./inc

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
lib: $(LibOut)


$(OBJS) : obj/%.o:src/%.cpp inc/%.h

$(LibSrc): ./inc/MyGuiMainFrame.h ./inc/LinkDef.h
		@echo "Generating dictionary $@..."
		$(ROOTCLING) -f $@ -c $^

$(LibOut): $(LibObj)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)


$(MAIN): %:%.o ${OBJS} ${LibObj}
		$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) -lGui -lGeom -lRGL
		@cp ./lib/*.pcm ./
		@echo "----> $@ done!"

#RICH: $(RICHO) 
#		$(LD) $(LDFLAGS) $^ $(GLIBS) -lEG -lHtml  $(OutPutOpt)$@
#		@echo "----> $@ done!"

clean:
		@rm -f $(OBJS)
		@rm -f *Dict*
		@rm -f ./src/*Dict*
		@rm -f ./lib/*
		@rm -f *.o 
		@rm -f ./obj/*.o
#
#%.o:%.C
#`		$(CXX) $(CXXFLAGS) -c $< -o $@

$(MAIN).o : $(MAIN).cpp $(HHDS)
	$(CXX) -o $@ $(CXXFLAGS) -c $<

obj/%.o : src/%.cpp
		$(CXX) -o $@ $(CXXFLAGS) -c $<

obj/%.o : lib/%.cpp
		$(CXX) -o $@ $(CXXFLAGS) -c $<