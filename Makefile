CXX            = g++
INCLUDES       = -I.
STRICT         = -pedantic -Wall -Wno-variadic-macros
EXTRA_DEBUG    = -g
EXTRA_CXXFLAGS = -c -Wall $(EXTRA_DEBUG) $(STRICT) $(INCLUDES)
EXTRA_LDFLAGS  = $(EXTRA_DEBUG) -lnetpbm -ljpeg -lopencv_core -lopencv_highgui -lopencv_imgproc

ifndef DEBUG
EXTRA_DEBUG   += -DNODEBUG
endif

ifeq ($(shell uname), Linux)
EXTRA_LDFLAGS +=  -static
endif

HEADERS        = error.h image.h measure.h 
SOURCE         = error.cc image.cc measure.cc main.cc ssim.cc
OBJECTS        = $(SOURCE:.cc=.o)
EXECUTABLE     = compare

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(HEADERS)
	$(CXX) $(LDFLAGS) $(EXTRA_LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -c $<

.cc :
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $< -o $@

clean:
	rm -fr $(EXECUTABLE).dSYM
	rm -f *~ *# *jpeg $(EXECUTABLE)
	rm -f *aux *bbl *blg *dvi *log *pdf *toc *snm *nav *out
	rm -f $(OBJECTS)

rebuild: clean all
