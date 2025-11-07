/*
 *  yosys -- Yosys Open SYnthesis Suite
 *
 *  Copyright (C) 2025  Shenzhen Pango Microsystems Co., Ltd.
 *
 *  Common inline functions for Pango FPGA synthesis
 */

#ifndef PANGO_COMMON_H
#define PANGO_COMMON_H

#include "kernel/yosys.h"

YOSYS_NAMESPACE_BEGIN

// Inline hot-path functions for performance
static inline bool IsNOT(RTLIL::Cell *cell)
{
	return (cell->type == ID($not) || cell->type == ID($_NOT_));
}

static inline bool IsAND(RTLIL::Cell *cell)
{
	return (cell->type == ID($and) || cell->type == ID($_AND_));
}

static inline bool IsOR(RTLIL::Cell *cell)
{
	return (cell->type == ID($or) || cell->type == ID($_OR_));
}

static inline bool IsXOR(RTLIL::Cell *cell)
{
	return (cell->type == ID($xor) || cell->type == ID($_XOR_));
}

static inline bool IsMUX(RTLIL::Cell *cell)
{
	return (cell->type == ID($mux) || cell->type == ID($_MUX_));
}

static inline bool IsGTP(RTLIL::Cell *cell) 
{ 
	return cell->type.begins_with("\\GTP_"); 
}

static inline bool IsCombinationalGate(RTLIL::Cell *cell) 
{ 
	return IsAND(cell) || IsOR(cell) || IsNOT(cell) || IsMUX(cell) || IsXOR(cell); 
}

YOSYS_NAMESPACE_END

#endif // PANGO_COMMON_H
