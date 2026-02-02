CXX = g++					# Specify which C++ compiler.
CXXFLAGS = -std=c++20 -O2	# Specify C++ 20, O2 optimization.

HTSLIB_PATH = /cbi/dabseq_v2/software/htslib-1.22.1
ZLIB_ROOT = /cbi/dabseq_v2/software/zlib

CXXFLAGS += -I$(HTSLIB_PATH)/include # Allows compiler to find `#include <htslib/file.h>
LDFLAGS = -L$(HTSLIB_PATH)/lib			# Allows linker to find `libhts.so` and `libhts.a` when you use `-lhts`.

LDLIBS = -lhts 							# Tells linker to link against hts library (installed in `/cbi/dabseq_v2/software`).

TARGET = main #main_orig
SRCS = fastq_reader.cpp barcode_index.cpp dabseq_utilities.cpp main.cpp #main_orig.cpp #main.cpp
OBJS = $(SRCS:.cpp=.o)

# $< outputs the first prerequisite
# $@ outputs target name
# $^ outputs all prerequisites

all: $(TARGET)

# Link step: executable from .o file.
# g++ <link_flags> -o <target_executable> <.o files> 
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# g++ <compile_flags> -c <.cpp> <.o file>
# Note, % is a general rule saying for .o targets, we have .cpp sources.
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: run clean


run: $(TARGET)
	@echo "Start: $$(date +'%Y-%m-%d %H:%M:%S %Z')"
	./$(TARGET) data/MB11_TS11_L004_R1_001.fastq.gz data/MB11_TS11_L004_R2_001.fastq.gz data/mb_cell_barcodes_v2.csv data/ab_barcodes.45plex_8.csv
	@echo "End:   $$(date +'%Y-%m-%d %H:%M:%S %Z')"

clean:
	rm -f $(TARGET) $(OBJS)
