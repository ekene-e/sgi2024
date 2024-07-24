MARCHING_TRIANGLES_TARGET=MarchingTriangles
MARCHING_TRIANGLES_SOURCE=MarchingTriangles/MarchingTriangles.cpp
MULTI_MARCHING_TRIANGLES_TARGET=MultiMarchingTriangles
MULTI_MARCHING_TRIANGLES_SOURCE=MultiMarchingTriangles/MultiMarchingTriangles.cpp
CURVE_TO_TUBE_TARGET=CurveToTube
CURVE_TO_TUBE_SOURCE=CurveToTube/CurveToTube.cpp
SMOOTH_CURVE_TARGET=SmoothCurve
SMOOTH_CURVE_SOURCE=SmoothCurve/SmoothCurve.cpp
PROCESS_GRID_TARGET=ProcessGrid
PROCESS_GRID_SOURCE=ProcessGrid/ProcessGrid.cpp
JITTER_GRID_TARGET=Jitter
JITTER_GRID_SOURCE=Jitter/Jitter.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++17 -Wno-invalid-offsetof -g
	LFLAGS += -lgomp -lstdc++ -lpthread
	CC=gcc
	CXX=g++
else
	CFLAGS += -Wno-deprecated -std=c++17 -Wno-invalid-offsetof -Wno-dangling-else -Wno-null-dereference -g
	LFLAGS += -lstdc++
	CC=clang
	CXX=clang++
endif

CFLAGS += -O3 -DRELEASE -funroll-loops -ffast-math -g
##LFLAGS += -O3 -g -lqhullstatic
LFLAGS += -O3 -g

BIN = Bin/Linux/
BIN_O = Obj/Linux/
INCLUDE = . -I /mnt/c/Research/Libraries/Include


MD=mkdir

MARCHING_TRIANGLES_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(MARCHING_TRIANGLES_SOURCE))))
MARCHING_TRIANGLES_OBJECT_DIR=$(dir $(MARCHING_TRIANGLES_OBJECTS))
MULTI_MARCHING_TRIANGLES_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(MULTI_MARCHING_TRIANGLES_SOURCE))))
MULTI_MARCHING_TRIANGLES_OBJECT_DIR=$(dir $(MULTI_MARCHING_TRIANGLES_OBJECTS))
CURVE_TO_TUBE_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(CURVE_TO_TUBE_SOURCE))))
CURVE_TO_TUBE_OBJECT_DIR=$(dir $(CURVE_TO_TUBE_OBJECTS))
SMOOTH_CURVE_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SMOOTH_CURVE_SOURCE))))
SMOOTH_CURVE_OBJECT_DIR=$(dir $(SMOOTH_CURVE_OBJECTS))
PROCESS_GRID_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(PROCESS_GRID_SOURCE))))
PROCESS_GRID_OBJECT_DIR=$(dir $(PROCESS_GRID_OBJECTS))
JITTER_GRID_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(JITTER_GRID_SOURCE))))
JITTER_GRID_OBJECT_DIR=$(dir $(JITTER_GRID_OBJECTS))

all: make_dirs
all: $(BIN)$(MARCHING_TRIANGLES_TARGET)
all: $(BIN)$(MULTI_MARCHING_TRIANGLES_TARGET)
all: $(BIN)$(CURVE_TO_TUBE_TARGET)
all: $(BIN)$(SMOOTH_CURVE_TARGET)
all: $(BIN)$(PROCESS_GRID_TARGET)
all: $(BIN)$(JITTER_GRID_TARGET)

MarchingTriangles: make_dirs
MarchingTriangles: $(BIN)$(MARCHING_TRIANGLES_TARGET)

MultiMarchingTriangles: make_dirs
MultiMarchingTriangles: $(BIN)$(MULTI_MARCHING_TRIANGLES_TARGET)

CurveToTube: make_dirs
CurveToTube: $(BIN)$(CURVE_TO_TUBE_TARGET)

SmoothCurve: make_dirs
SmoothCurve: $(BIN)$(SMOOTH_CURVE_TARGET)

ProcessGrid: make_dirs
ProcessGrid: $(BIN)$(PROCESS_GRID_TARGET)

JitterGrid: make_dirs
JitterGrid: $(BIN)$(JITTER_GRID_TARGET)

clean:
	rm -rf $(BIN)$(MARCHING_TRIANGLES_TARGET)
	rm -rf $(BIN)$(MULTI_MARCHING_TRIANGLES_TARGET)
	rm -rf $(BIN)$(CURVE_TO_TUBE_TARGET)
	rm -rf $(BIN)$(SMOOTH_CURVE_TARGET)
	rm -rf $(BIN)$(PROCESS_GRID_TARGET)
	rm -rf $(BIN)$(JITTER_GRID_TARGET)
	rm -rf $(BIN_O)

make_dirs: FORCE
	$(MD) -p $(BIN)
	$(MD) -p $(BIN_O)
	$(MD) -p $(MARCHING_TRIANGLES_OBJECT_DIR)
	$(MD) -p $(MULTI_MARCHING_TRIANGLES_OBJECT_DIR)
	$(MD) -p $(CURVE_TO_TUBE_OBJECT_DIR)
	$(MD) -p $(SMOOTH_CURVE_OBJECT_DIR)
	$(MD) -p $(PROCESS_GRID_OBJECT_DIR)
	$(MD) -p $(JITTER_GRID_OBJECT_DIR)

$(BIN)$(MARCHING_TRIANGLES_TARGET): $(MARCHING_TRIANGLES_OBJECTS)
	$(CXX) -o $@ $(MARCHING_TRIANGLES_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(MULTI_MARCHING_TRIANGLES_TARGET): $(MULTI_MARCHING_TRIANGLES_OBJECTS)
	$(CXX) -o $@ $(MULTI_MARCHING_TRIANGLES_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(CURVE_TO_TUBE_TARGET): $(CURVE_TO_TUBE_OBJECTS)
	$(CXX) -o $@ $(CURVE_TO_TUBE_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(SMOOTH_CURVE_TARGET): $(SMOOTH_CURVE_OBJECTS)
	$(CXX) -o $@ $(SMOOTH_CURVE_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(PROCESS_GRID_TARGET): $(PROCESS_GRID_OBJECTS)
	$(CXX) -o $@ $(PROCESS_GRID_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(JITTER_GRID_TARGET): $(JITTER_GRID_OBJECTS)
	$(CXX) -o $@ $(JITTER_GRID_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE: