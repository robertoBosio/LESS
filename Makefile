#
# Copyright 2019-2020 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

MK_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CUR_DIR := $(patsubst %/,%,$(dir $(MK_PATH)))
XF_PROJ_ROOT ?= $(shell bash -c 'export MK_PATH=$(MK_PATH); echo $${MK_PATH%L1/*}')
XPART := "xck26-sfvc784-2LV-c"

.PHONY: help

help::
	@echo ""
	@echo "Makefile Usage:"
	@echo ""
	@echo "  make runhls CSIM=1 CSYNTH=1 COSIM=1"
	@echo "      Command to run the selected tasks for specified device."
	@echo ""
	@echo "      Valid tasks are CSIM, CSYNTH, COSIM, EXPORT, IMPL"
	@echo ""
	@echo "  make clean "
	@echo "      Command to remove the generated files."
	@echo ""

.PHONY: check_vivado
check_vivado:
ifeq (,$(wildcard $(XILINX_VIVADO)/bin/vivado))
	@echo "Cannot locate Vivado installation. Please set XILINX_VIVADO variable." && false
endif

.PHONY: check_hls
check_hls:
ifeq (,$(wildcard $(XILINX_HLS)/bin/vitis_hls))
	@echo "Cannot locate Vitis_hls installation. Please set XILINX_HLS variable." && false
endif

export PATH := $(XILINX_VIVADO)/bin:$(PATH)
export PATH := $(XILINX_VITIS)/bin:$(XILINX_XRT)/bin:$(PATH)
.PHONY: run setup runhls clean cleanall check

# Alias to run, for legacy test script
check: run

CSIM ?= 0
CSYNTH ?= 0
COSIM ?= 0
VIVADO_SYN ?= 0
VIVADO_IMPL ?= 0
EXPORT ?= 0
IMPL ?= 0
QOR_CHECK ?= 0

# at least RTL synthesis before check QoR
ifeq (1,$(QOR_CHECK))
ifeq (0,$(VIVADO_IMPL))
override VIVADO_SYN := 1
endif
endif

# need synthesis before cosim or vivado
ifeq (1,$(VIVADO_IMPL))
override CSYNTH := 1
endif

ifeq (1,$(VIVADO_SYN))
override CSYNTH := 1
endif

ifeq (1,$(COSIM))
override CSYNTH := 1
endif

ifeq (1,$(EXPORT))
override CSYNTH := 1
endif

ifeq (1,$(IMPL))
override CSYNTH := 1
override EXPORT := 1
endif

run: setup runhls runimpl

setup:
	@rm -f ./scripts/settings.tcl
	@if [ -n "$$CLKP" ]; then echo 'set CLKP $(CLKP)' >> ./scripts/settings.tcl ; fi
	@echo 'set XPART $(XPART)' >> ./scripts/settings.tcl
	@echo 'set CSIM $(CSIM)' >> ./scripts/settings.tcl
	@echo 'set CSYNTH $(CSYNTH)' >> ./scripts/settings.tcl
	@echo 'set COSIM $(COSIM)' >> ./scripts/settings.tcl
	@echo 'set VIVADO_SYN $(VIVADO_SYN)' >> ./scripts/settings.tcl
	@echo 'set VIVADO_IMPL $(VIVADO_IMPL)' >> ./scripts/settings.tcl
	@echo 'set EXPORT $(EXPORT)' >> ./scripts/settings.tcl
	@echo 'set XF_PROJ_ROOT "$(XF_PROJ_ROOT)"' >> ./scripts/settings.tcl
	@echo 'set CUR_DIR "$(CUR_DIR)"' >> ./scripts/settings.tcl
	@echo "Configured: scripts/settings.tcl"
	@echo "----"
	@cat ./scripts/settings.tcl
	@echo "----"

HLS ?= vitis_hls
runhls: setup | check_hls
	$(HLS) -f scripts/run_hls.tcl;

runimpl: setup | check_vivado
	vivado -mode tcl -source scripts/script_vivado.tcl

XSA := subgraph.xsa
MK_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
COMMON_REPO ?= $(shell bash -c 'export MK_PATH=$(MK_PATH); echo $${MK_PATH%/*}')
PWD = $(shell readlink -f .)

TARGET ?= hw
TOP_NAME = subgraphIsomorphism
TEMP_DIR := $(CUR_DIR)/_x.$(TARGET).$(XSA)
BUILD_DIR := $(CUR_DIR)/build_dir.$(TARGET).$(XSA)
LINK_OUTPUT := $(BUILD_DIR)/$(TOP_NAME).link.xclbin
PACKAGE_OUT = ./package.$(TARGET)

VPP_PFLAGS :=
VPP_LDFLAGS :=

############################## Setting up Host Variables ##############################
CMD_ARGS = -x $(BUILD_DIR)/$(TOP_NAME).xclbin

#Include Required Host Source Files
CXXFLAGS += -I$(CUR_DIR)/source
CXXFLAGS += -I$(XILINX_XRT)/include
CXXFLAGS += -I$(XILINX_VIVADO)/include
CXXFLAGS += -I$(XILINX_HLS)/include
CXXFLAGS += -Wno-attributes
CXXFLAGS += -Wno-unknown-pragmas
CXXFLAGS += -O3
CXXFLAGS += -g
# CXXFLAGS += -std=c++17
CXXFLAGS += -fmessage-length=0
CXXFLAGS += -DXILINX_XRT

HOST_SRCS := $(CUR_DIR)/source/cmdlineparser.cpp
HOST_SRCS += $(CUR_DIR)/source/logger.cpp
HOST_SRCS += $(CUR_DIR)/source/subgraphIsomorphism.cpp

# Host compiler global settings
LDFLAGS += -lrt -lstdc++
LDFLAGS += -luuid -lxrt_coreutil -pthread -lopencv_imgproc -lopencv_core -lopencv_imgcodecs
LDFLAGS += -L$(XILINX_XRT)/lib -pthread -lOpenCL

EXECUTABLE = $(CUR_DIR)/host.o
EMCONFIG_DIR = $(TEMP_DIR)

############################## Kernel Source Files Repository##########################
emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(PLATFORM) --od $(EMCONFIG_DIR)

.PHONY: test
test: $(EXECUTABLE) emconfig check-target check-top
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	cd $(CUR_DIR) && cp $(EMCONFIG_DIR)/emconfig.json .
	cd $(CUR_DIR) && XCL_EMULATION_MODE=$(TARGET) $(EXECUTABLE) $(CMD_ARGS)
else
	$(EXECUTABLE) $(CMD_ARGS)
endif

.PHONY: host
host: $(EXECUTABLE)

############################## Setting up Kernel Variables ##############################
# Kernel compiler global settings
VPP_FLAGS += --save-temps
VPP_FLAGS += --config subiso_hls.cfg
# VPP_FLAGS += --profile_kernel data:all:all:all

VPP_LDFLAGS += --vivado.synth.jobs 4
VPP_LDFLAGS += --vivado.impl.jobs 4
VPP_LDFLAGS += --config subiso_link.cfg
# VPP_LDFLAGS += --profile_kernel data:all:all:all

$(TEMP_DIR)/$(TOP_NAME).xo: $(CUR_DIR)/source/subgraphIsomorphism.cpp
	mkdir -p $(TEMP_DIR)
	v++ -g -c $(VPP_FLAGS) -t $(TARGET) --platform $(PLATFORM) -k $(TOP_NAME) --temp_dir $(TEMP_DIR) -o'$@' '$<'

$(BUILD_DIR)/$(TOP_NAME).xclbin: $(TEMP_DIR)/$(TOP_NAME).xo
	mkdir -p $(BUILD_DIR)
	v++ -g -l  $(VPP_FLAGS) $(VPP_LDFLAGS) -t $(TARGET) --platform $(PLATFORM) --temp_dir $(TEMP_DIR) -o $(LINK_OUTPUT) $(+)
	v++ -p $(LINK_OUTPUT) $(VPP_FLAGS) -t $(TARGET) --platform $(PLATFORM) --package.out_dir $(PACKAGE_OUT) -o $(BUILD_DIR)/$(TOP_NAME).xclbin

$(EXECUTABLE): $(HOST_SRCS)
	g++ -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -rf settings.tcl *_hls.log vivado*.log vivado*.jou vivado*.str

cleanall: clean
	rm -rf less_hls less_bd