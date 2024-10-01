VPP ?= ${XILINX_VITIS}/bin/v++
TARGET ?= hw

RM = rm -f
RMDIR = rm -rf

VPP_OPTS = --target $(TARGET)

HW_KERNEL_OBJS += spot_finder.xo
SRC_KERNEL_OBJS += spot_finder.xo

ALL_MESSAGE_FILES = $(subst .xo,.o,.mdb,$(SRC_KERNEL_OBJS))

.PHONY: all clean
all: $(HW_KERNEL_OBJS)

clean:
	-$(RM) $(SRC_KERNEL_OBJS) $(ALL_MESSAGE_FILES) 
	-$(RMDIR) $(BUILD_SUBDIRS)
	-$(RMDIR) .Xil

.PHONY: incremental
incremental: all

nothing:

check_env:
ifndef PLATFORM
	$(error PLATFORM not provided)
endif

spot_finder.xo: spot_finder.cc check_env
	-@mkdir -p $(@D)
	-@$(RM) $@
	$(VPP) $(VPP_OPTS) --compile -I"$(<D)" --platform $(PLATFORM) --kernel spot_finder -o"$@" "$<"

