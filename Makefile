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
	@rm -f ./settings.tcl
	@if [ -n "$$CLKP" ]; then echo 'set CLKP $(CLKP)' >> ./settings.tcl ; fi
	@echo 'set XPART $(XPART)' >> ./settings.tcl
	@echo 'set CSIM $(CSIM)' >> ./settings.tcl
	@echo 'set CSYNTH $(CSYNTH)' >> ./settings.tcl
	@echo 'set COSIM $(COSIM)' >> ./settings.tcl
	@echo 'set VIVADO_SYN $(VIVADO_SYN)' >> ./settings.tcl
	@echo 'set VIVADO_IMPL $(VIVADO_IMPL)' >> ./settings.tcl
	@echo 'set EXPORT $(EXPORT)' >> ./settings.tcl
	@echo 'set XF_PROJ_ROOT "$(XF_PROJ_ROOT)"' >> ./settings.tcl
	@echo 'set CUR_DIR "$(CUR_DIR)"' >> ./settings.tcl
	@echo "Configured: settings.tcl"
	@echo "----"
	@cat ./settings.tcl
	@echo "----"

HLS ?= vitis_hls
runhls: setup | check_hls
	$(HLS) -f run_hls.tcl;

runimpl: setup | check_vivado
	vivado -mode tcl -source script_vivado.tcl

clean:
	rm -rf settings.tcl *_hls.log vivado*.log vivado*.jou vivado*.str

cleanall: clean
