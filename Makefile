#object_files (testGene.o)
OBJECTS := $(wildcard *.o)

#root_stuff (root libraries and needed root options)
ROOTLIBS := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore
ROOTCINT := $(shell which rootcint)

#exe_files
EXECUTABLE6  := saveEffiGeneReFit
CLASS1       := RooGenePdf
CLASSDICT1   := $(CLASS1)Dictionary.cxx
CLASS2       := RooBernsteinEffi
CLASSDICT2   := $(CLASS2)Dictionary.cxx
CLASS3       := RooProdGenePdfBernsteinEffi
CLASSDICT3   := $(CLASS3)Dictionary.cxx

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS1 := $(CLASS1).cxx  $(CLASSDICT1)
LIBS2 := $(CLASS2).cxx  $(CLASSDICT2)
LIBS3 := $(CLASS3).cxx  $(CLASSDICT3)

	
all: $(CLASSDICT1) $(CLASSDICT2) $(CLASSDICT3) $(EXECUTABLE6)

dict: $(CLASSDICT1) $(CLASSDICT2) $(CLASSDICT3) 

dict1: $(CLASSDICT1)

dict2: $(CLASSDICT2)

dict3: $(CLASSDICT3)

dict4: $(CLASSDICT4)

plot: $(EXMAKEPLOT)

list: $(EXREADLIST)

$(CLASSDICT1): $(CLASS1).h $(CLASS1)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^

$(CLASSDICT2): $(CLASS2).h $(CLASS2)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^

$(CLASSDICT3): $(CLASS3).h $(CLASS3)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^

$(CLASSDICT4): $(CLASS4).h $(CLASS4)LinkDef.h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^

	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS2) $(LIBS3) $(ROOTLIBS) $(ROOTFLAGS) -I.

$(EXECUTABLE6): $(EXECUTABLE6).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS1) $(LIBS2) $(LIBS3) $(ROOTLIBS) $(ROOTFLAGS) -I.



#cleaning options
.PHONY: clean cleanall
clean:
	rm -f $(OBJECTS) && rm -f $(EXECUTABLE6) *Dictionary.cxx *Dictionary_rdict.pcm
cleanplot:
	rm -f  $(EXMAKEPLOT)

