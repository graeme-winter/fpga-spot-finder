VPP_LINKER ?= ${XILINX_VITIS}/bin/v++
TARGET ?= hw

RM = rm -f
RMDIR = rm -rf

VPP_OPTS = --target $(TARGET)

BINARY_CONTAINERS += spot_finder.xclbin

OBJECTS += ../kernel/spot_finder.xo

ALL_MESSAGE_FILES = $(subst .xclbin,.mdb,$(BINARY_CONTAINERS))

.PHONY: all clean
all: $(BINARY_CONTAINERS)

clean:
	-$(RM) $(BINARY_CONTAINERS) $(ALL_MESSAGE_FILES) 
	-$(RM) *.xclbin.sh *.xclbin.info *.xclbin.link_summary*
	-$(RMDIR) $(BUILD_SUBDIRS)
	-$(RMDIR) .Xil

ifeq ($(MAKECMDGOALS),incremental)
.PHONY: spot_finder.xclbin
.PHONY: ../kernel/spot_finder.xo
endif

.PHONY: incremental
incremental: all

nothing:

check_env:
ifndef PLATFORM
	$(error PLATFORM not provided)
endif

spot_finder.xclbin: $(OBJECTS) check_env
	$(VPP_LINKER) $(VPP_OPTS) --link --platform $(PLATFORM) -o"$@" $(OBJECTS)
