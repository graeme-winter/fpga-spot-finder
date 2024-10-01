RM = rm -f
RMDIR = rm -rf

CXXFLAGS += -std=c++1y -I$(XILINX_XRT)/include -I$(XILINX_HLS)/include \
  -DVITIS_PLATFORM=$(PLATFORM) -O2 -Wall -c -fmessage-length=0
LDFLAGS += -L$(XILINX_XRT)/lib -luuid -lxrt_coreutil -lxilinxopencl \
  -lpthread -lrt -ldl -lcrypt -lstdc++

HOST_OBJECTS = spot_finder.o

HOST_EXE = spot_finder
HOST_CXX = g++

.PHONY: all clean
all:  $(HOST_EXE)

clean:
	-$(RM) $(HOST_EXE) $(HOST_OBJECTS)

.PHONY: incremental
incremental: all

nothing:

check_env:
ifndef PLATFORM
	$(error PLATFORM not provided)
endif

host.o: host.cc check_env
	-@mkdir -p $(@D)
	$(HOST_CXX) $(CXXFLAGS) -o "$@" "$<"

$(HOST_EXE): $(HOST_OBJECTS)
	$(HOST_CXX) -o "$@" $(+) $(LDFLAGS)
