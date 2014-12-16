CC = g++
LD = g++
AR = ar

# compile
INCLUDES = /usr/local/boost src /Users/tbabb/code/gltoy/src
CFLAGS   = -std=c++11 -O3 -Wall -c -fmessage-length=0 -Wno-unused -Wno-unused-local-typedefs
IFLAGS   = $(addprefix -I, $(INCLUDES))

# link
LIBDIRS  = /usr/local/boost/stage/lib \
           /System/Library/Frameworks/OpenGL.framework/Libraries \
	   /Users/tbabb/code/gltoy/lib
LIBS     = geomc GL png SDL_mixer z GLU boost_system gltoy
LDFLAGS  = $(addprefix -l, $(LIBS)) \
           $(addprefix -L, $(LIBDIRS)) \
           -framework GLUT -framework OpenGL \
           `/opt/local/bin/sdl-config --libs --cflags --static-libs` \

# sources
MODULES  = imu
SRC      = $(wildcard src/*.cpp) \
           $(foreach m, $(MODULES), $(wildcard src/$(m)/*.cpp))
OBJ      = $(patsubst src/%.cpp, build/%.o, $(SRC))


all: sim ksim qsim

clean:
	rm -rf ./build/*

## binaries

qsim:  build/qsim.o
	$(CC) $(LDFLAGS) build/qsim.o -o bin/qsim

ksim: build/KinematicSimulator.o build/KinematicSolver.o build/KinematicSimulator.o build/KinSim.o
	$(CC) $(LDFLAGS) build/KinematicSolver.o build/KinematicSimulator.o build/KinSim.o -o bin/ksim

sim: $(OBJ)
	$(CC) $(LDFLAGS) build/RigidBody.o build/KeyIntegrator.o build/sim.o -o bin/sim

build/%.o : src/%.cpp
	mkdir -p $(patsubst src/%, build/%, $(dir $<))
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

