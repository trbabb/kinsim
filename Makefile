CC = clang++
LD = clang++
AR = ar

# compile
INCLUDES = /usr/local/boost src /Users/tbabb/code/gltoy/src
CFLAGS   = -g -std=c++11 -O3 -Wall -c -fmessage-length=0 -Wno-unused -Wno-unused-local-typedefs
IFLAGS   = $(addprefix -I, $(INCLUDES))

# link
LIBDIRS  = /System/Library/Frameworks/OpenGL.framework/Libraries \
	   /Users/tbabb/code/gltoy/lib
LIBS     = geomc GL png z GLU gltoy
LDFLAGS  = $(addprefix -l, $(LIBS)) \
           $(addprefix -L, $(LIBDIRS)) \
           -framework GLUT -framework OpenGL

# sources
MODULES  = imu
SRC      = $(wildcard src/*.cpp) \
           $(foreach m, $(MODULES), $(wildcard src/$(m)/*.cpp))
OBJ      = $(patsubst src/%.cpp, build/%.o, $(SRC))


all: sim ksim qsim

clean:
	rm -rf ./build/*

## binaries

sensors: build/SerialSensor.o
	$(CC) $(LDFLAGS) build/SerialSensor.o -o bin/sensors

rotsim: build/RigidBody.o  build/rotsim.o
	$(CC) $(LDFLAGS) build/RigidBody.o build/rotsim.o -o bin/rotsim

qsim:  build/qsim.o
	$(CC) $(LDFLAGS) build/qsim.o -o bin/qsim

ksim: build/KinematicSimulator.o build/KinematicSolver.o build/KinematicSimulator.o build/KinSim.o
	$(CC) $(LDFLAGS) build/KinematicSolver.o build/KinematicSimulator.o build/KinSim.o -o bin/ksim

kaltest: build/kaltest.o
	$(CC) $(LDFLAGS) build/kaltest.o -o bin/kaltest

sim: build/RigidBody.o build/KeyIntegrator.o build/sim.o
	$(CC) $(LDFLAGS) build/RigidBody.o build/KeyIntegrator.o build/sim.o -o bin/sim

build/%.o : src/%.cpp
	mkdir -p $(patsubst src/%, build/%, $(dir $<))
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

