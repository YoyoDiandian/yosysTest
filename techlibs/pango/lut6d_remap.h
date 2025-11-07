#ifndef LUT6D_REMAP_H
#define LUT6D_REMAP_H

#include "kernel/yosys.h"
#include "kernel/sigtools.h"

YOSYS_NAMESPACE_BEGIN

// Exposed globals from lut6d_remap.cc
// extern dict<Cell *, size_t> lut2levels;
// extern dict<size_t, pool<Cell *>> level2luts;

// // Public APIs implemented in lut6d_remap.cc
// void GetTopoSortedLUTs(Module *module, vector<Cell *> &luts);

// // Build a GTP_LUT6D cell (dual-output LUT), used by both remap and map flows
// RTLIL::Cell *addLut6d(Module *module, const string &name,
//                       const vector<SigBit> &inputs,
//                       const vector<SigBit> &outputs,
//                       Const INIT = Const(0, 64));

// // Convenience IO helpers for LUT-style cells
// void GetCellOutputsVector(Cell *cell, vector<SigBit> &outputs);

// // Check if two single-output LUTs can be merged into a LUT6D with given I5
// std::pair<Const, vector<SigBit>> checkCombinable(Cell *lut1, Cell *lut2, SigBit I5);

// // Update local data structures after a merge
// void updateDataStructuresAfterMerge(Module *module, Cell *old_lut1, Cell *old_lut2, Cell *new_lut6d);

// // Perform a concrete merge attempt with a candidate and common inputs
// Cell *performLUT6DMerge(Module *module, Cell *lut1, Cell *lut2, pool<SigBit> &common_inputs);

// // Per-root remap driver and top-level remapper
// void LUT6DRemapMain(Module *module, Cell *root);
void LUT6DRemap(Module *module);

YOSYS_NAMESPACE_END

#endif