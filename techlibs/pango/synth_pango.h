#pragma once

#include "kernel/yosys.h"
#include "kernel/sigtools.h"
#include "techlibs/pango/pango_common.h"  // Inline hot-path functions

YOSYS_NAMESPACE_BEGIN

// Extern globals provided by synth_pango.cc (or related synthesis helpers)
extern dict<Cell*, vector<SigBit>> cell2bits;
extern dict<SigBit, Cell*> bit2driver;
extern dict<SigBit, vector<Cell*>> bit2reader;
extern dict<Cell*, dict<pool<SigBit>, pool<Cell*>>> cell2cuts;

// Some passes/utilities also rely on a shared SigMap
extern SigMap sigmap;

// Mapper configuration flags and controls
extern bool using_internel_lut_type;
extern size_t MAX_INTERATIONS;

// Mapper main driver
void MapperMain(Module *module);

// Helper APIs used across LUT6D remap/map implementation
void GetPrimeInputOuput(Module *module, pool<SigBit> &prime_inputs, pool<SigBit> &prime_outputs);

bool IsCombinationalCell(Cell *cell);
bool IsGTP_LUT6D(Cell *cell);

// Inline functions are now in pango_common.h:
// - IsCombinationalGate, IsAND, IsOR, IsNOT, IsMUX, IsXOR, IsGTP

// Generic cell IO helpers
RTLIL::SigBit GetCellOutput(Cell *cell);
void GetCellInputsVector(Cell *cell, vector<SigBit> &inputs);
void GetCellInputsSet(Cell *cell, pool<SigBit> &inputs);

// Flow helpers
void MapperInit(Module *module);
void CheckCellWidth(Module *module);
void GetTopoSortedGates(Module *module, vector<Cell*> &gates);
void GenerateCuts(Module *module);
void GenerateCuts(Cell *cell);

YOSYS_NAMESPACE_END
