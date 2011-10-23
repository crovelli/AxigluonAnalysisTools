ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lMinuit2
GLIBS          = $(filter-out -lNew, $(NGLIBS))

# For MT2:
# MT2LIB_INC_DIR := mT2
# MT2LIB_LIB_DIR := mT2/src
# MT2LIB_CPPFLAGS := -I $(MT2LIB_INC_DIR)
# MT2LIB_LDFLAGS  := $(MT2LIB_LIB_DIR)/libMt2.so
# CXXFLAGS += $(MT2LIB_CPPFLAGS)
# LDFLAGS  += $(MT2LIB_LDFLAGS)

INCLUDEDIR       = ./
INCLUDEDIRCOMMON = ../
CXX	         += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/
OUTLIBCOMMON     = $(INCLUDEDIRCOMMON)/CommonTools/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)AxigluonBase.o: $(INCLUDEDIR)/src/AxigluonBase.C \
	$(INCLUDEDIR)/src/Axigluon.cc \
	$(INCLUDEDIR)/src/AxigluonSelection.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)AxigluonBase.o $<
$(OUTLIBCOMMON)Conditions.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Conditions.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Conditions.o $<
$(OUTLIBCOMMON)PUWeight.o: $(INCLUDEDIRCOMMON)/CommonTools/src/PUWeight.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)PUWeight.o $<
$(OUTLIBCOMMON)LumiReWeighting.o: $(INCLUDEDIRCOMMON)/CommonTools/src/LumiReWeighting.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)LumiReWeighting.o $<
$(OUTLIBCOMMON)Utils.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Utils.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIBCOMMON)Utils.o $<
$(OUTLIBCOMMON)Counters.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Counters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Counters.o $<
$(OUTLIBCOMMON)Selection.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Selection.o $<
$(OUTLIBCOMMON)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EfficiencyEvaluator.o $<
$(OUTLIBCOMMON)Monitor.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Monitor.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Monitor.o $<
$(OUTLIBCOMMON)SprDataFiller.o: $(INCLUDEDIRCOMMON)/CommonTools/src/SprDataFiller.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)SprDataFiller.o $<
$(OUTLIBCOMMON)TriggerMask.o: $(INCLUDEDIRCOMMON)/CommonTools/src/TriggerMask.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)TriggerMask.o $<
$(OUTLIB)EcalCleaner.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/EcalCleaner.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)EcalCleaner.o $<
$(OUTLIB)CutBasedEleIDSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)CutBasedEleIDSelector.o $<
$(OUTLIB)CiCBasedEleSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CiCBasedEleSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)CiCBasedEleSelector.o $<
$(OUTLIB)Axigluon.o: $(INCLUDEDIR)/src/Axigluon.cc $(OUTLIB)JetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Axigluon.o $<
$(OUTLIB)JetCorrectorParameters.o: $(INCLUDEDIRCOMMON)/HiggsAnalysisTools/src/JetCorrectorParameters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)JetCorrectorParameters.o $<
$(OUTLIB)SimpleJetCorrectionUncertainty.o: $(INCLUDEDIRCOMMON)/HiggsAnalysisTools/src/SimpleJetCorrectionUncertainty.cc \
	$(OUTLIB)JetCorrectorParameters.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleJetCorrectionUncertainty.o $<
$(OUTLIB)JetCorrectionUncertainty.o: $(INCLUDEDIRCOMMON)/HiggsAnalysisTools/src/JetCorrectionUncertainty.cc \
	$(OUTLIB)SimpleJetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectionUncertainty.o $<
$(OUTLIB)AxigluonSelection.o: $(INCLUDEDIR)/src/AxigluonSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)AxigluonSelection.o $<
$(OUTLIB)RedAxiTree.o: $(INCLUDEDIR)/src/RedAxiTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedAxiTree.o $<
$(OUTLIB)CutBasedAxiSelector.o: $(INCLUDEDIR)/src/CutBasedAxiSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutBasedAxiSelector.o $<
$(OUTLIB)QGLikelihoodCalculator.o: $(INCLUDEDIR)/src/QGLikelihoodCalculator.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)QGLikelihoodCalculator.o $<

# ==================== HiggsApp =============================================
AxigluonApp: $(INCLUDEDIR)/src/AxigluonApp.C \
	$(OUTLIB)AxigluonBase.o \
	$(OUTLIB)Axigluon.o \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)PUWeight.o \
	$(OUTLIBCOMMON)LumiReWeighting.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)Monitor.o \
	$(OUTLIBCOMMON)SprDataFiller.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIB)CutBasedEleIDSelector.o \
	$(OUTLIB)QGLikelihoodCalculator.o \
	$(OUTLIB)CiCBasedEleSelector.o \
	$(OUTLIB)RedAxiTree.o \
	$(OUTLIB)EcalCleaner.o \
	$(OUTLIB)CutBasedAxiSelector.o 
	$(CXX) $(CXXFLAGS) -ldl -o AxigluonApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(GLIBS) $(LDFLAGS) $ $<
AxigluonApp.clean:
	rm -f AxigluonApp

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f AxigluonApp

all:  AxigluonApp
