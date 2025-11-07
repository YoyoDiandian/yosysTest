#pragma once

#include "kernel/yosys.h"

YOSYS_NAMESPACE_BEGIN

// Public entry point for LUT6D mapping (implemented in lut6d_map.cc)
void LUT6DMapper(Module *module);

YOSYS_NAMESPACE_END
