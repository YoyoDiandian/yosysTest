/*
 * LUT6D Merge Implementation
 * Based on the algorithm for merging single-output LUTs into dual-output LUT6D
 */

#include "kernel/yosys.h"
#include "kernel/sigtools.h"
#include "kernel/celltypes.h"
#include "kernel/modtools.h"
#include "kernel/satgen.h"
#include "techlibs/pango/synth_pango.h"
// #include "techlibs/pango/lut6d_mapping.h"
#include <queue>
#include <algorithm>

USING_YOSYS_NAMESPACE
using namespace std;
YOSYS_NAMESPACE_BEGIN

// Use different name to avoid conflict with synth_pango.cc's bit2States
dict<SigBit, vector<State>> g_lut6d_bit2States; // for combination pre-check in LUT6D mapping
// Track visited by name to avoid holding stale pointers across clean()
// pool<IdString> visited;
// Precomputed original combinational fanout count per cell (reader cells count)
static dict<Cell*, int> g_orig_fanout;
static bool g_fanout_ready = false;

// Dynamic parameters that will be adjusted based on netlist size
static int MAX_CUTS_PER_GATE = 8;        // how many cuts per gate to consider at most
static int MAX_CANDS_PER_CUT = 16;       // how many candidates per cut to consider at most
static int MAX_CAND_CUTS_PER_CAND = 8;   // how many cand_cuts per candidate to consider at most
static int MAX_TFO_DEPTH = 8;            // limit TFO exploration depth from cut bits
static int MAX_TFO_CELLS = 256;          // max cells collected in TFO per cut
static int MAX_CONE_CELLS = 64;          // skip cand_cut whose cone has too many cells
static int PATTERN_STAGE1 = 32;          // staged precheck patterns before full scan
static int MAX_LEVEL_DELTA = 8;          // |level(gate)-level(cand)| window


// Unified number of random patterns for partial simulation
static int NUM_PATTERNS = 128;           // number of random patterns for simulation

// ULTRA FAST MODE: Skip expensive SAT verification, rely on pattern matching
// This trades some accuracy for ~10-15% speedup
// static bool ENABLE_SAT_VERIFICATION = false;  // Set to true for highest quality, false for speed

// Store bit states indexed by signal name to survive clean operations
static dict<string, vector<State>> g_bit_states_cache;

// Current module pointer to validate cell pointers after clean/removal
static Module *g_current_module = nullptr;

// Forward declaration
pool<SigBit> CombinationPreCheck(const SigBit &Z1, const SigBit &Z2, const pool<SigBit> &cut1, const pool<SigBit> &cut2);

// Dynamic parameter adjustment based on netlist size
// More fine-grained adjustment: 5 tiers from tiny to huge netlists
// Principle: Bigger netlists need more aggressive pruning for acceptable runtime
static void AdjustParametersBasedOnNetlistSize(int num_combinational_gates) {
    log("  Netlist has %d combinational gates\n", num_combinational_gates);
    
    if (num_combinational_gates < 3000) {
        // Tiny netlist (<3K): maximum quality, original aggressive parameters
        log("  TINY netlist (<3K gates) - using MAXIMUM QUALITY parameters\n");
        MAX_CUTS_PER_GATE = 16;
        MAX_CANDS_PER_CUT = 32;
        MAX_CAND_CUTS_PER_CAND = 16;
        MAX_TFO_DEPTH = 12;
        MAX_TFO_CELLS = 2048;
        MAX_CONE_CELLS = 128;
        PATTERN_STAGE1 = 64;
        MAX_LEVEL_DELTA = 12;
        NUM_PATTERNS = 256;
    } else if (num_combinational_gates < 10000) {
        // Small netlist (3K-10K): high quality with slight optimization
        log("  SMALL netlist (3K-10K gates) - using HIGH QUALITY parameters\n");
        MAX_CUTS_PER_GATE = 14;
        MAX_CANDS_PER_CUT = 28;
        MAX_CAND_CUTS_PER_CAND = 14;
        MAX_TFO_DEPTH = 11;
        MAX_TFO_CELLS = 1536;
        MAX_CONE_CELLS = 112;
        PATTERN_STAGE1 = 56;
        MAX_LEVEL_DELTA = 11;
        NUM_PATTERNS = 224;
    } else if (num_combinational_gates < 30000) {
        // Medium netlist (10K-30K): balanced quality vs speed
        log("  MEDIUM netlist (10K-30K gates) - using BALANCED parameters\n");
        MAX_CUTS_PER_GATE = 12;
        MAX_CANDS_PER_CUT = 24;
        MAX_CAND_CUTS_PER_CAND = 12;
        MAX_TFO_DEPTH = 10;
        MAX_TFO_CELLS = 1024;
        MAX_CONE_CELLS = 96;
        PATTERN_STAGE1 = 48;
        MAX_LEVEL_DELTA = 10;
        NUM_PATTERNS = 192;
    } else if (num_combinational_gates < 80000) {
        // Large netlist (30K-80K): favor speed over quality
        log("  LARGE netlist (30K-80K gates) - using SPEED-FAVORED parameters\n");
        MAX_CUTS_PER_GATE = 10;
        MAX_CANDS_PER_CUT = 20;
        MAX_CAND_CUTS_PER_CAND = 10;
        MAX_TFO_DEPTH = 9;
        MAX_TFO_CELLS = 512;
        MAX_CONE_CELLS = 80;
        PATTERN_STAGE1 = 40;
        MAX_LEVEL_DELTA = 9;
        NUM_PATTERNS = 160;
    } else {
        // Huge netlist (>80K): maximum speed optimization
        log("  HUGE netlist (>80K gates) - using MAXIMUM SPEED parameters\n");
        MAX_CUTS_PER_GATE = 8;
        MAX_CANDS_PER_CUT = 16;
        MAX_CAND_CUTS_PER_CAND = 8;
        MAX_TFO_DEPTH = 8;
        MAX_TFO_CELLS = 256;
        MAX_CONE_CELLS = 64;
        PATTERN_STAGE1 = 32;
        MAX_LEVEL_DELTA = 8;
        NUM_PATTERNS = 128;
    }
    
    log("  Parameters: MAX_CUTS_PER_GATE=%d, MAX_CANDS_PER_CUT=%d, MAX_CAND_CUTS_PER_CAND=%d\n",
        MAX_CUTS_PER_GATE, MAX_CANDS_PER_CUT, MAX_CAND_CUTS_PER_CAND);
    log("  Parameters: MAX_TFO_DEPTH=%d, MAX_TFO_CELLS=%d, MAX_CONE_CELLS=%d\n",
        MAX_TFO_DEPTH, MAX_TFO_CELLS, MAX_CONE_CELLS);
    log("  Parameters: PATTERN_STAGE1=%d, MAX_LEVEL_DELTA=%d, NUM_PATTERNS=%d\n",
        PATTERN_STAGE1, MAX_LEVEL_DELTA, NUM_PATTERNS);
}

// Track cells that have been cleaned/removed to avoid processing
static pool<IdString> g_removed_cell_names;

// Cache for TFO checks to avoid repeated BFS traversals
static dict<pair<SigBit, SigBit>, bool> g_tfo_cache;

// Cache for fanin cone computations
static dict<SigBit, pair<pool<Cell*>, pool<SigBit>>> g_fanin_cone_cache;

// Cache for cell output - critical hot path optimization
static dict<Cell*, SigBit> g_cell_output_cache;

// Cache for SAT verification results
struct SatCacheKey {
    Cell* z5;
    Cell* z;
    SigBit I5;
    bool operator==(const SatCacheKey& other) const {
        return z5 == other.z5 && z == other.z && I5 == other.I5;
    }
    unsigned int hash() const {
        // Combine hashes of pointers and SigBit
        unsigned int h = 0;
        h = ((uintptr_t)z5) ^ (((uintptr_t)z) << 7);
        // For SigBit, combine wire pointer and offset
        if (I5.wire) {
            h = h ^ (((uintptr_t)I5.wire) << 14) ^ (I5.offset << 21);
        }
        return h;
    }
    // Required by Yosys hashlib - work with Hasher object
    template<typename Hasher>
    Hasher hash_into(Hasher h) const {
        h.hash32((uintptr_t)z5);
        h.hash32((uintptr_t)z);
        if (I5.wire) {
            h.hash32((uintptr_t)I5.wire);
            h.hash32(I5.offset);
        }
        return h;
    }
};
static dict<SatCacheKey, bool> g_sat_cache;

// 是否会被移除的节点（它的output不可用）
pool<SigBit> removed_sigbit;

// 已经被合并
pool<Cell *> visited;

// Helper function to check if a cut contains any removed signal bits
static inline bool CutHasRemovedBits(const pool<SigBit> &cut) {
    // Early exit if removed_sigbit is empty (common case)
    if (removed_sigbit.empty()) return false;
    
    // Use smaller set for iteration
    if (removed_sigbit.size() < cut.size()) {
        for (const auto &bit : removed_sigbit) {
            if (cut.count(bit)) return true;
        }
    } else {
        for (const auto &bit : cut) {
            if (removed_sigbit.count(bit)) return true;
        }
    }
    return false;
}

// Fast cached version of GetCellOutput - CRITICAL HOT PATH
// This function is called 2.5+ billion times, so every cycle counts
static inline SigBit GetCellOutputCached(Cell *cell) {
    // Use cache lookup first - much faster than cell2bits lookup
    auto it = g_cell_output_cache.find(cell);
    if (it != g_cell_output_cache.end()) {
        return it->second;
    }
    // Fallback to uncached version (should rarely happen)
    SigBit out = GetCellOutputCached(cell);
    g_cell_output_cache[cell] = out;
    return out;
}

// simple level map for combinational cells (topo depth)
static dict<Cell*, int> g_cell_level;
static void BuildCellLevels(Module *module)
{
    g_cell_level.clear();
    vector<Cell*> order;
    GetTopoSortedGates(module, order);
    std::reverse(order.begin(), order.end()); // inputs -> outputs
    for (Cell *c : order) {
        if (!IsCombinationalGate(c)) continue;
        int max_in = 0;
        vector<SigBit> inbits;
        GetCellInputsVector(c, inbits);
        for (auto &b : inbits) {
            Cell *drv = bit2driver.count(b) ? bit2driver[b] : nullptr;
            if (drv && IsCombinationalGate(drv)) {
                if (g_cell_level.count(drv)) max_in = std::max(max_in, g_cell_level[drv]);
            }
        }
        g_cell_level[c] = max_in + 1;
    }
}

// Wrap the ordered candidate structure for readability and helpers.
// vector< pair< cut, vector< pair< cand, vector< pair< cand_cut, cone > > > > > >
// list: [cut: [(cand, [(cand_cut, cone), ...]), ...]]
struct NestedCandMap {
    using cut = pool<SigBit>;
    using cand_cut = pool<SigBit>;
    using potentialI5 = pool<SigBit>;
    using CandCuts = vector<pair<cand_cut, potentialI5>>;
    using CandEntry = pair<Cell*, CandCuts>;
    using OuterEntry = pair<cut, vector<CandEntry>>;

    vector<OuterEntry> data;

    NestedCandMap() = default;
    explicit NestedCandMap(const vector<OuterEntry>& v) : data(v) {}

    size_t size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    void clear() { data.clear(); }
    void push_back(const OuterEntry& x) { data.push_back(x); }

    OuterEntry& operator[](size_t i) { return data[i]; }
    const OuterEntry& operator[](size_t i) const { return data[i]; }

    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }
};

// Structure to hold all merge-related information
struct MergeInfo {
    Cell *gate;
    Cell *cand_merge;
    vector<SigBit> LUT6Dinput;
    vector<bool> LUT6Dinit;
    SigBit Z;
    SigBit Z5;
    // pool<Cell*> cells_to_remove;  // Cells that will be removed after merge

    MergeInfo() : gate(nullptr), cand_merge(nullptr), Z(), Z5() {}
    
    MergeInfo(Cell *g, Cell *c, const vector<SigBit> &inputs, const vector<bool> &init, const SigBit &z, const SigBit &z5)
        : gate(g), cand_merge(c), LUT6Dinput(inputs), LUT6Dinit(init), Z(z), Z5(z5) {}
};

// Soft evaluation: identical logic but never logs errors. Returns Sx on any unknown/unsupported.
State StateEvalSoft(dict<SigBit, State> &bit_map, SigBit out)
{
	// Handle constants directly
	if (!out.is_wire()) {
		State const_state = out.data;
		bit_map[out] = const_state;
		return const_state;
	}
	
	if (bit_map.count(out)) {
        return bit_map[out];
    }
	if (!bit2driver.count(out)) return State::Sx;
	Cell *cell = bit2driver[out];
	if (!cell) return State::Sx;
    if (!cell2bits.count(cell)) return State::Sx;
	// Handle LUT6D cells
	if (cell->type == ID(GTP_LUT6D)) {
		// Get INIT parameter
		if (!cell->hasParam(ID(INIT))) return State::Sx;
		Const init_const = cell->getParam(ID(INIT));
		if (init_const.size() != 64) return State::Sx;
		
		// Get input signals
		vector<SigBit> inputs(6);
		inputs[0] = sigmap(cell->getPort(ID(I0)))[0];
		inputs[1] = sigmap(cell->getPort(ID(I1)))[0];
		inputs[2] = sigmap(cell->getPort(ID(I2)))[0];
		inputs[3] = sigmap(cell->getPort(ID(I3)))[0];
		inputs[4] = sigmap(cell->getPort(ID(I4)))[0];
		inputs[5] = sigmap(cell->getPort(ID(I5)))[0];
		
		// Evaluate input states
		int addr = 0;
		bool has_unknown = false;
		for (int i = 0; i < 6; i++) {
			State in_state = StateEvalSoft(bit_map, inputs[i]);
			if (in_state == State::Sx) {
				has_unknown = true;
				break;
			}
			if (in_state == State::S1) addr |= (1 << i);
		}
		
		if (has_unknown) return bit_map[out] = State::Sx;
		
		// Look up INIT value at computed address
		// State result = init_const[addr];
		
		// Determine which output we're evaluating (Z or Z5)
		SigBit z_bit = sigmap(cell->getPort(ID(Z)))[0];
		SigBit z5_bit = sigmap(cell->getPort(ID(Z5)))[0];
		
		if (out == z_bit) {
			// Z output (when I5=1 in standard LUT6D configuration)
			// For now, use the INIT value directly
			bit_map[out] = init_const[addr];
			return bit_map[out];
		} else if (out == z5_bit) {
			// Z5 output (when I5=0)
			// For now, use the INIT value directly
            if (addr > 31) {
                bit_map[out] = init_const[addr - 32];
            }
            else {
                bit_map[out] = init_const[addr];
            }
			return bit_map[out];
		}
		
		return bit_map[out] = State::Sx;
	}

	vector<SigBit> bits = cell2bits[cell];
	
	if (IsAND(cell)) {
		for (size_t i = 1; i < bits.size(); i++) {
            // log("\t%s is not in bitmap, it navigates to AND %s\n", log_signal(out), log_id(cell->name));
			State tmp = StateEvalSoft(bit_map, bits[i]);
			if (tmp == State::S0) return bit_map[bits[0]] = State::S0;
			if (tmp == State::Sx) return bit_map[bits[0]] = State::Sx;
		}
		return bit_map[bits[0]] = State::S1;
	} else if (IsOR(cell)) {
		for (size_t i = 1; i < bits.size(); i++) {
            // log("\t%s is not in bitmap, it navigates to OR %s\n", log_signal(out), log_id(cell->name));
			State tmp = StateEvalSoft(bit_map, bits[i]);
			if (tmp == State::S1) return bit_map[bits[0]] = State::S1;
			if (tmp == State::Sx) return bit_map[bits[0]] = State::Sx;
		}
		return bit_map[bits[0]] = State::S0;
	} else if (IsNOT(cell)) {
        // log("\t%s is not in bitmap, it navigates to NOT %s\n", log_signal(out), log_id(cell->name));
		State tmp = StateEvalSoft(bit_map, bits[1]);
		if (tmp == State::S0) return bit_map[bits[0]] = State::S1;
		if (tmp == State::S1) return bit_map[bits[0]] = State::S0;
		return bit_map[bits[0]] = State::Sx;
	} else if (IsMUX(cell)) {
		SigSpec sig = cell->getPort(ID(S)); sig = sigmap(sig); if (!sig.size()) return State::Sx;
		SigBit sel_bit = sig[0];
		SigBit i0_bit = sigmap(cell->getPort(ID(A)))[0];
		SigBit i1_bit = sigmap(cell->getPort(ID(B)))[0];
        // log("\t%s is not in bitmap, it navigates to MUX %s\n", log_signal(out), log_id(cell->name));
		State sel_state = StateEvalSoft(bit_map, sel_bit);
		if (sel_state == State::S0) {
            // log("\t%s is not in bitmap, it navigates to MUX i0 %s\n", log_signal(out), log_id(cell->name));
            return bit_map[bits[0]] = StateEvalSoft(bit_map, i0_bit);
        }
		if (sel_state == State::S1) {
            // log("\t%s is not in bitmap, it navigates to  MUX i1 %s\n", log_signal(out), log_id(cell->name));
            return bit_map[bits[0]] = StateEvalSoft(bit_map, i1_bit);
        }
		return bit_map[bits[0]] = State::Sx;
	} else if (IsXOR(cell)) {
        // log("\t%s is not in bitmap, it navigates to XOR %s\n", log_signal(out), log_id(cell->name));
		State a = StateEvalSoft(bit_map, bits[1]);
        // log("\t%s is not in bitmap, it navigates to XOR %s\n", log_signal(out), log_id(cell->name));
		State b = StateEvalSoft(bit_map, bits[2]);
		if (a == State::Sx || b == State::Sx) return bit_map[bits[0]] = State::Sx;
		if (a == b) return bit_map[bits[0]] = State::S0;
		return bit_map[bits[0]] = State::S1;
	} else {
		// unsupported types in soft 
		log("  for output %s, unhandled cell %s \n", log_signal(out), cell->type.c_str());
		return State::Sx;
	}
}

void partialSimulator(const pool<SigBit> &seed_inputs) {
    for (SigBit bit : seed_inputs) {
        vector<State> states;
        states.reserve(NUM_PATTERNS);
        for (int i = 0; i < NUM_PATTERNS; i++) {
            State random_state = (rand() % 2 == 0) ? State::S0 : State::S1;
            states.push_back(random_state);
        }
        g_lut6d_bit2States[bit] = states;
    }
}

void GenerateBitState(Module *module, const pool<SigBit> &prime_inputs, bool incremental) {
    // Only clear if not incremental update
    if (!incremental) {
        g_lut6d_bit2States.clear();
        g_bit_states_cache.clear();
    }
    // 先基于简单门区域的“边界信号”生成输入模式：这些比特的驱动不是可处理的组合门（例如 PI、寄存器输出、黑盒等）
    // pool<SigBit> seeds = prime_inputs;
    // vector<Cell *> gates_tmp;
    // GetTopoSortedGates(module, gates_tmp);
    // for (Cell *cell : gates_tmp) {
    //     if (!IsCombinationalGate(cell)) continue;
    //     pool<SigBit> inputs;
    //     GetCellInputsSet(cell, inputs);
    //     for (auto &ib : inputs) {
    //         Cell *drv = bit2driver.count(ib) ? bit2driver[ib] : nullptr;
    //         if (!drv || !IsCombinationalGate(drv)) {
    //             seeds.insert(ib);
    //         }
    //     }
    // }
    partialSimulator(prime_inputs);
	
	// 得到部分真值表
    vector<Cell *> gates;
    GetTopoSortedGates(module, gates);

	for (Cell *cell : gates) {
		if (!IsCombinationalGate(cell)) {
			continue;
		}
		
		SigBit output = GetCellOutputCached(cell);
        vector<State> output_states;
        output_states.reserve(NUM_PATTERNS);

		pool<SigBit> inputs;
		GetCellInputsSet(cell, inputs);
		
		// int num_of_zero = 0;
		// int num_of_one = 0;

        for (int i = 0; i < NUM_PATTERNS; i++) {
			dict<SigBit, State> bit_map;
			
			for (SigBit input_bit : inputs) {
				if (g_lut6d_bit2States.count(input_bit)) {
					bit_map[input_bit] = g_lut6d_bit2States[input_bit][i];
				}
			}
			
            State output_state = StateEvalSoft(bit_map, output);
			output_states.push_back(output_state);
			
			// if (output_state == State::S0) {
			// 	num_of_zero++;
			// } else if (output_state == State::S1) {
			// 	num_of_one++;
			// }
		}
		
		// 检查0和1的数量是否都大于5
		// if (num_of_zero > 5 && num_of_one > 5) {
		// 	log("Cell %s output has balanced states: %d zeros, %d ones\n", 
		// 		cell->name.c_str(), num_of_zero, num_of_one);
		// }
		// else {
		// 	log("Cell %s output has unbalanced states: %d zeros, %d ones\n", 
		// 		cell->name.c_str(), num_of_zero, num_of_one);
		// }
        g_lut6d_bit2States[output] = output_states;
        
        // Cache the states using signal name for persistence across clean operations
        string sig_name = log_signal(output);
        g_bit_states_cache[sig_name] = output_states;
	}
}

// Slow baseline (fallback)
void ComputeMFFC(Cell *root, pool<Cell*> &mffc)
{
    mffc.clear();
    if (!root || !IsCombinationalGate(root)) return;
    // 构建 fanout 计数
    dict<Cell*, int> fanout_cnt;
    for (auto &kv : cell2bits) {
        Cell *c = kv.first;
        if (!IsCombinationalGate(c)) continue;
        vector<SigBit> inbits; 
        GetCellInputsVector(c, inbits);
        for (auto &b : inbits) {
            Cell *drv = bit2driver.count(b) ? bit2driver[b] : nullptr;
            if (drv && IsCombinationalGate(drv)) fanout_cnt[drv]++;
        }
    }
    // DFS 栈式扩展
    vector<Cell*> stack; stack.push_back(root);
    while (!stack.empty()) {
        Cell *c = stack.back(); stack.pop_back();
        if (!c || !IsCombinationalGate(c) || mffc.count(c)) continue;
        mffc.insert(c);
        vector<SigBit> inbits;
        GetCellInputsVector(c, inbits);
        for (auto &b : inbits) {
            Cell *drv = bit2driver.count(b) ? bit2driver[b] : nullptr;
            if (!drv || !IsCombinationalGate(drv) || mffc.count(drv)) continue;
            if (fanout_cnt[drv] == 1) {
                vector<SigBit> pp; 
                GetCellInputsVector(drv, pp);
                for (auto &bb : pp) {
                    Cell *pdrv = bit2driver.count(bb) ? bit2driver[bb] : nullptr;
                    if (pdrv && IsCombinationalGate(pdrv)) fanout_cnt[pdrv]--;
                }
                stack.push_back(drv);
            }
        }
    }
}

// Build original combinational fanout map once (per module init)
static void BuildOrigFanout()
{
    g_orig_fanout.clear();
    for (auto &kv : cell2bits) {
        Cell *c = kv.first;
        if (!IsCombinationalGate(c)) continue;
        vector<SigBit> inbits; GetCellInputsVector(c, inbits);
        for (auto &b : inbits) {
            Cell *drv = bit2driver.count(b) ? bit2driver[b] : nullptr;
            if (drv && IsCombinationalGate(drv)) g_orig_fanout[drv]++;
        }
    }
    g_fanout_ready = true;
}

// Fast MFFC using precomputed original fanout and local deltas
static void ComputeMFFCFast(Cell *root, pool<Cell*> &mffc)
{
    mffc.clear();
    if (!root || !IsCombinationalGate(root)) return;
    if (!g_fanout_ready) BuildOrigFanout();

    dict<Cell*, int> local_cnt; // only touched nodes
    auto get_cnt = [&](Cell *c)->int {
        if (!c || !IsCombinationalGate(c)) return 0;
        if (local_cnt.count(c)) return local_cnt[c];
        return g_orig_fanout.count(c) ? g_orig_fanout[c] : 0;
    };
    auto dec_cnt = [&](Cell *c){
        if (!c || !IsCombinationalGate(c)) return;
        int v = get_cnt(c);
        local_cnt[c] = v - 1;
    };

    vector<Cell*> stack; stack.push_back(root);
    while (!stack.empty()) {
        Cell *c = stack.back(); stack.pop_back();
        if (!c || !IsCombinationalGate(c) || mffc.count(c)) continue;
        mffc.insert(c);
        vector<SigBit> inbits;
        GetCellInputsVector(c, inbits);
        for (auto &b : inbits) {
            Cell *drv = bit2driver.count(b) ? bit2driver[b] : nullptr;
            if (!drv || !IsCombinationalGate(drv) || mffc.count(drv)) continue;
            if (get_cnt(drv) == 1) {
                // include drv and propagate decrement to its inputs' drivers
                vector<SigBit> pp;
                GetCellInputsVector(drv, pp);
                for (auto &bb : pp) {
                    Cell *pdrv = bit2driver.count(bb) ? bit2driver[bb] : nullptr;
                    if (pdrv && IsCombinationalGate(pdrv)) dec_cnt(pdrv);
                }
                stack.push_back(drv);
            }
        }
    }
}

// 计算一组比特的组合传递扇出 TFO（到达的组合 Cell 集合）
static void ComputeTFO(const pool<SigBit> &src_bits, pool<Cell*> &tfo_cells, int max_depth = -1)
{
    tfo_cells.clear();
    // 从比特出发，沿 bit -> reader cell -> cell output bit -> 下一层 递推
    deque<pair<SigBit,int>> qbits; // (bit, depth)
    pool<SigBit> seen_bits;
    for (auto b : src_bits) {
        qbits.push_back({b, 0});
        seen_bits.insert(b);
    }
    while (!qbits.empty()) {
        auto [b, d] = qbits.front(); qbits.pop_front();
        if (!bit2reader.count(b)) continue;
        for (Cell *c : bit2reader[b]) {
            if (!c || !IsCombinationalGate(c)) continue;
            if (!tfo_cells.count(c)) {
                if ((int)tfo_cells.size() >= MAX_TFO_CELLS) {
                    // log("ComputeTFO达到最大TFO单元限制，停止进一步处理。\n");
                    continue; // cap TFO size
                }
                tfo_cells.insert(c);

            }
            // 有可能已经被合并了，cell2bits中的这个cell被删除了。
            if (cell2bits.count(c) == 0) continue;
            // 推进到它的输出 bit
            SigBit ob = GetCellOutputCached(c);
            if (!seen_bits.count(ob)) {
                seen_bits.insert(ob);
                if (max_depth < 0 || d+1 <= max_depth)
                    qbits.push_back({ob, d+1});
                else {
                    // 达到深度限制，停止推进
                    // log("ComputeTFO达到深度限制，停止推进Cell %s\n", c->name.c_str());
                }
            }
        }
    }
}

// 统计两个 pool 的交并信息（只需要大小）
static inline pair<int,int> InterUnionSize(const pool<SigBit> &A, const pool<SigBit> &B)
{
    int inter = 0;
    // 遍历较小集合以加速
    if (A.size() <= B.size()) {
        for (auto &x : A) if (B.count(x)) inter++;
    } else {
        for (auto &x : B) if (A.count(x)) inter++;
    }
    int uni = int(A.size() + B.size() - inter);
    return {inter, uni};
}

// 计算输出比特 out 的完整组合 fan-in 锥（带缓存）：
// - cone 返回锥内所有组合 cell（包含 out 的驱动 cell 若为组合门）
// - leaves 返回锥的叶子比特（无驱动或驱动为非组合单元/PI/常量/寄存器输出等）
// - max_depth: 限制BFS深度，防止在大型CPU设计中追溯过深（默认30层）
// - max_cells: 限制cone中cell数量，防止SAT公式膨胀（默认300个cell）
static void GetFullFanInCone(SigBit out, pool<Cell*> &cone, pool<SigBit> &leaves, int max_depth = 30, int max_cells = 300)
{
    // Check cache first
    if (g_fanin_cone_cache.count(out)) {
        const auto &cached = g_fanin_cone_cache[out];
        cone = cached.first;
        leaves = cached.second;
        return;
    }
    
    cone.clear();
    leaves.clear();

    // 以比特为节点进行后向遍历：bit -> driver cell -> input bits -> ...
    // 使用pair<SigBit, depth>来跟踪深度
    deque<pair<SigBit, int>> work;
    pool<SigBit> seen_bits;
    work.push_back({out, 0});
    seen_bits.insert(out);

    while (!work.empty()) {
        auto [b, depth] = work.front();
        work.pop_front();

        Cell *drv = bit2driver.count(b) ? bit2driver[b] : nullptr;
        // 没有驱动或驱动不是组合门：视为叶子
        if (!drv || !IsCombinationalGate(drv)) {
            leaves.insert(b);
            continue;
        }

        // 深度限制：达到最大深度时，将当前bit视为叶子
        if (depth >= max_depth) {
            leaves.insert(b);
            continue;
        }

        // Cell数量限制：防止cone过大导致SAT爆炸
        if (cone.size() >= (size_t)max_cells) {
            leaves.insert(b);
            continue;
        }

        // 加入锥体并扩展其输入
        if (!cone.count(drv)) {
            cone.insert(drv);
            vector<SigBit> inbits;
            GetCellInputsVector(drv, inbits);
            for (auto &ib : inbits) {
                if (!seen_bits.count(ib)) {
                    seen_bits.insert(ib);
                    work.push_back({ib, depth + 1});
                }
            }
        }
    }
    
    // Store in cache (aggressive caching for large designs)
    if (g_fanin_cone_cache.size() < 100000) {  // increase from 50000 -> 100000
        g_fanin_cone_cache[out] = {cone, leaves};
    }
}

// Check whether two cells have shared fanin cone (either shared logic cells or shared input signals)
// Returns true if they share at least one combinational cell in their fanin or share at least one leaf signal
// 增强版本：如果cone重叠过多（>80%），也返回true以拒绝合并（防止SAT爆炸）
// static bool HasSharedFaninCone(Cell *gate, Cell *cand)
// {
//     if (!gate || !cand) return false;
    
//     pool<Cell*> gate_cone, cand_cone;
//     pool<SigBit> gate_leaves, cand_leaves;
    
//     // Get full fanin cone for both cells (with limits)
//     GetFullFanInCone(GetCellOutputCached(gate), gate_cone, gate_leaves, 20, 200);
//     GetFullFanInCone(GetCellOutputCached(cand), cand_cone, cand_leaves, 20, 200);
    
//     // Check if cones have intersection (shared intermediate logic)
//     size_t shared_cells = 0;
//     for (auto *c : gate_cone) {
//         if (cand_cone.count(c)) shared_cells++;
//     }
    
//     // 优化：如果cone重叠度过高（>80%），拒绝合并，避免SAT复杂度爆炸
//     // 这种情况通常出现在CPU的ALU/寄存器文件等全局共享逻辑中
//     size_t min_cone_size = std::min(gate_cone.size(), cand_cone.size());
//     if (min_cone_size > 0 && shared_cells > 0) {
//         float overlap_ratio = (float)shared_cells / (float)min_cone_size;
//         if (overlap_ratio > 0.8f) {
//             // log("  Rejecting merge due to high cone overlap: %zu/%zu (%.1f%%)\n", 
//             //     shared_cells, min_cone_size, overlap_ratio * 100);
//             return true; // 拒绝：视为"有共享"而过滤
//         }
//         return true; // 有共享逻辑
//     }
    
//     // Check if leaves have intersection (shared input signals)
//     for (auto &b : gate_leaves) {
//         if (cand_leaves.count(b)) return true;
//     }
    
//     return false;
// }

// cand 列表：vector< pair<cand, vector< pair<cand_cut, cone_cells> >> >，并按规则有序
void getPotentialCandidatesForCut(Cell *gate, const pool<SigBit> &cut_bits,
                                  vector<pair<Cell*, vector<pair<pool<SigBit>, pool<SigBit>>>>> &cands,
                                  dict<Cell*, pool<Cell*>> &mffc_cache)
{
    cands.clear();
    if (!gate || !IsCombinationalGate(gate)) return;

    // 1) gate 的 MFFC（用于 cand 排除）
    pool<Cell*> mffc_gate; ComputeMFFCFast(gate, mffc_gate);

    // 2) cut 的 TFO（cand 必须在其中），限制搜索深度
    pool<Cell*> tfo_cells; ComputeTFO(cut_bits, tfo_cells, MAX_TFO_DEPTH);

    const int K = 6;
    const bool cut_is_full = (int)cut_bits.size() == K;
    
    // 优化：预先获取gate输出，避免重复调用
    SigBit Z_gate = GetCellOutputCached(gate);

    // 3) 遍历 TFO 里的候选 cand
    // bool early_stop = false;
    int cand_count = 0;
    for (Cell *cand : tfo_cells) {
        if (!cand || !IsCombinationalGate(cand)) continue;
        if (g_current_module && g_current_module->cell(cand->name) != cand) continue; // skip stale cand
        if (cand == gate) continue;
        if (mffc_gate.count(cand)) continue; // cand 不能属于 gate 的 MFFC
        if (visited.count(cand)) continue;
        
        // 优化：预先获取cand输出，避免重复调用
        SigBit Z_cand = GetCellOutputCached(cand);
        if (cut_bits.count(Z_cand)) continue; // cand 输出不能在 cut 内

        // Level window 过滤，减少远距离候选
        if (g_cell_level.count(gate) && g_cell_level.count(cand)) {
            int dg = std::abs(g_cell_level[gate] - g_cell_level[cand]);
            if (dg > MAX_LEVEL_DELTA) {
                // log("Level window filtering reached maximum delta (%d > %d), skipping candidate %s\n", dg, MAX_LEVEL_DELTA, cand->name.c_str());
                continue;
            }
        }
        // 4) 收集 cand 的 cuts，并做筛选 + 排序
        if (!cell2cuts.count(cand) || cell2cuts[cand].empty()) continue;

        struct Item { pool<SigBit> cut; pool<SigBit> potential_I5; int uni; int inter; pool<SigBit> i5_gate; pool<SigBit> i5_cand; };
        vector<Item> items; items.reserve(cell2cuts[cand].size());

        int candCut_count = 0;
        for (auto &kv : cell2cuts[cand]) {
            const pool<SigBit> &cand_cut = kv.first;

            // 优化：按成本从低到高排列检查顺序
            // 1. 最快：整数比较
            if ((int)kv.second.size() > MAX_CONE_CELLS) {
                continue;
            }
            // 2. 较快：单个元素查找
            if (cand_cut.count(Z_gate)) {
                continue;
            }
            // 3. 较慢：需要迭代比较
            // INLINED CutHasRemovedBits for performance
            bool has_removed = false;
            if (!removed_sigbit.empty()) {
                if (removed_sigbit.size() < cand_cut.size()) {
                    for (const auto &bit : removed_sigbit) {
                        if (cand_cut.count(bit)) {
                            has_removed = true;
                            break;
                        }
                    }
                } else {
                    for (const auto &bit : cand_cut) {
                        if (removed_sigbit.count(bit)) {
                            has_removed = true;
                            break;
                        }
                    }
                }
            }
            if (has_removed) {
                continue;
            }
            
            // Skip if cand_cut and cut_bits have no intersection
            auto [inter, uni] = InterUnionSize(cut_bits, cand_cut);
            if (inter == 0 || uni > K) continue; // 合并k-可行检查

            // 4.2) 约束：对于 cut==6，要 cand_cut ⊆ cut_bits 且 cand_cut != cut_bits
            if (cut_is_full) {
                // 优化：先检查大小，避免不必要的循环
                if (cand_cut.size() >= cut_bits.size()) continue;
                
                // 检查 cand_cut 是否为 cut_bits 的子集
                bool subset = true;
                for (auto &x : cand_cut) {
                    if (!cut_bits.count(x)) { 
                        subset = false; 
                        break; 
                    }
                }
                if (!subset) continue;
            }

            // 4.3) 预筛：CombinationPreCheck（仅当 uni==K 时才计算）
            if (uni == K) {
                auto potential_I5 = CombinationPreCheck(Z_gate, Z_cand, cut_bits, cand_cut);
                if (potential_I5.empty()) continue; // 若需要 I5 且无可行 I5，则跳过
                items.push_back(Item{cand_cut, potential_I5, uni, inter});
            }
            else {
                // uni < K 时，直接接受，无需 I5
                items.push_back(Item{cand_cut, pool<SigBit>(), uni, inter});
            }
            candCut_count++;
            if (candCut_count >= MAX_CAND_CUTS_PER_CAND) {
                // log("  Reached maximum candidate cuts per candidate (%d), stopping further cut processing.\n", MAX_CAND_CUTS_PER_CAND);
                break;
            }
        }

        if (items.empty()) continue;

        // 4.3) cand 内部 cand_cuts 排序（或部分排序）：并集大优先，其次交集大优先，其次 cand_cut 小优先
        auto comp = [](const Item &a, const Item &b){
            if (a.uni != b.uni) return a.uni > b.uni;         // 并集越大越优
            if (a.inter != b.inter) return a.inter > b.inter; // 交集越大越优
            return a.cut.size() < b.cut.size();               // 相同时，cut 越小越优
        };
        if ((int)items.size() > MAX_CAND_CUTS_PER_CAND) {
            std::partial_sort(items.begin(), items.begin()+MAX_CAND_CUTS_PER_CAND, items.end(), comp);
            items.resize(MAX_CAND_CUTS_PER_CAND);
        } else {
            std::sort(items.begin(), items.end(), comp);
        }

        // 4.4) 早停条件：如果存在 merge_cut.size()<6（且通过 precheck），则收下当前 cand 并全局早停
        bool found_sub6 = false;
        for (auto &it : items) {
            if (it.uni < K) {
                found_sub6 = true; 
                break; 
            }
        }

        // 4.5) gate 不能属于 cand 的 MFFC（昂贵）：延后到这里，仅对有候选 cand 的情况做一次
        // 优化：使用引用避免拷贝
        pool<Cell*> *mffc_cand_ptr;
        pool<Cell*> mffc_cand_local;
        if (mffc_cache.count(cand)) {
            mffc_cand_ptr = &mffc_cache[cand];
        } else { 
            ComputeMFFCFast(cand, mffc_cand_local); 
            mffc_cache[cand] = std::move(mffc_cand_local);
            mffc_cand_ptr = &mffc_cache[cand];
        }
        const pool<Cell*> &mffc_cand = *mffc_cand_ptr;
        if (mffc_cand.count(gate)) {
            if (cand_count >= MAX_CANDS_PER_CUT) {
                // log("  Reached maximum candidates per cut (%d), stopping further candidate processing.\n", MAX_CANDS_PER_CUT);
                break;
            }
            continue;
        }

        // 4.6) 生成 cand 的有序 cand_cuts 结果
        // 优化：直接构造，避免中间临时对象
        cands.emplace_back();
        auto &new_entry = cands.back();
        new_entry.first = cand;
        new_entry.second.reserve(items.size());
        for (auto &it : items) {
            new_entry.second.emplace_back(std::move(it.cut), std::move(it.potential_I5));
        }
        cand_count++;

        if (found_sub6) {
            // early_stop = true;
            break;
        }
        if (cand_count >= MAX_CANDS_PER_CUT) {
            // log("  Reached maximum candidates per cut (%d), stopping further candidate processing.\n", MAX_CANDS_PER_CUT);
            break;
        }
    }

    if (cands.empty()) return;

    // 5) 候选 cand 列表排序：按其“最优 cand_cut”的评分排序（同一比较规则）
    auto best_key = [&](const vector<pair<pool<SigBit>, pool<SigBit>>> &cuts){
        int best_uni = -1, best_inter = -1; size_t best_size = 0;
        for (auto &p : cuts) {
            auto [inter, uni] = InterUnionSize(cut_bits, p.first);
            if (uni > best_uni || (uni == best_uni && (inter > best_inter || (inter == best_inter && p.first.size() < best_size)))) {
                best_uni = uni; best_inter = inter; best_size = p.first.size();
            }
        }
        return tuple<int,int,size_t>(best_uni, best_inter, best_size);
    };
    std::sort(cands.begin(), cands.end(), [&](const pair<Cell*, vector<pair<pool<SigBit>, pool<SigBit>>>> &A,
                                              const pair<Cell*, vector<pair<pool<SigBit>, pool<SigBit>>>> &B){
        auto ka = best_key(A.second);
        auto kb = best_key(B.second);
        if (std::get<0>(ka) != std::get<0>(kb)) return std::get<0>(ka) > std::get<0>(kb);
        if (std::get<1>(ka) != std::get<1>(kb)) return std::get<1>(ka) > std::get<1>(kb);
        return std::get<2>(ka) < std::get<2>(kb);
    });
}

// 检查目标bit是否在源bit的TFO中 - 高性能优化版本（带缓存）
// 返回true表示target在source的TFO中
static bool IsInTFO(SigBit source, SigBit target, int max_depth = -1)
{
    if (source == target) return true; // 自身算在TFO中
    
    // Cache查询
    auto cache_key = std::make_pair(source, target);
    if (g_tfo_cache.count(cache_key)) {
        return g_tfo_cache[cache_key];
    }
    
    // 优化1：快速排除 - 如果source没有reader，不可能到达target
    if (!bit2reader.count(source)) {
        g_tfo_cache[cache_key] = false;
        return false;
    }
    
    // 优化2：使用队列而非deque，减少内存分配开销
    deque<pair<SigBit, int>> work;
    pool<SigBit> seen;
    work.push_back({source, 0});
    seen.insert(source);
    
    // 优化3：限制遍历的节点数量，避免在大图中无限扩展
    int nodes_visited = 0;
    const int MAX_NODES = 500;  // reduce from 1000 -> 500 for faster early exit
    
    while (!work.empty()) {
        auto [b, d] = work.front();
        work.pop_front();
        
        // 早停1：找到目标
        if (b == target) {
            g_tfo_cache[cache_key] = true;
            return true;
        }
        
        // 早停2：深度限制
        if (max_depth >= 0 && d >= max_depth) continue;
        
        // 早停3：节点数限制
        if (++nodes_visited > MAX_NODES) {
            g_tfo_cache[cache_key] = false;
            return false;
        }
        
        if (!bit2reader.count(b)) continue;
        
        const auto &readers = bit2reader[b];
        // 优化4：预先分配空间，减少动态扩展
        for (Cell *c : readers) {
            if (!c || !IsCombinationalGate(c)) continue;
            if (cell2bits.count(c) == 0) continue;
            
            SigBit ob = GetCellOutputCached(c);
            
            // 早停4：立即检查是否是目标，避免加入队列
            if (ob == target) {
                g_tfo_cache[cache_key] = true;
                return true;
            }
            
            if (!seen.count(ob)) {
                seen.insert(ob);
                work.push_back({ob, d+1});
            }
        }
    }
    
    g_tfo_cache[cache_key] = false;
    return false;
}

// 预筛查是否能够合并 - 高性能优化版本
pool<SigBit> CombinationPreCheck(const SigBit &Z1, const SigBit &Z2, const pool<SigBit> &cut1, const pool<SigBit> &cut2) {
    // potential I5 初始为 cut1 与 cut2 的对称差
    pool<SigBit> diff_bits = cut1;
    for (SigBit bit : cut2) {
        if (diff_bits.count(bit)) diff_bits.erase(bit);
        else diff_bits.insert(bit);
    }

    pool<SigBit> potential_I5 = diff_bits;
    
    // 早期返回：如果没有候选I5，直接返回
    if (potential_I5.empty()) return potential_I5;
    
    // 早期TFO过滤：用一次反向TFI代替多次正向TFO
    // 遍历所有I5候选，确定z5_output，然后检查I5是否在z5的TFI中
    pool<SigBit> valid_I5;
    for (const SigBit &i5_cand : potential_I5) {
        bool i5_in_cut1 = cut1.count(i5_cand);
        bool i5_in_cut2 = cut2.count(i5_cand);
        
        SigBit z5_output;
        if (i5_in_cut1 && !i5_in_cut2) {
            z5_output = Z2;
        } else if (!i5_in_cut1 && i5_in_cut2) {
            z5_output = Z1;
        } else {
            continue; // Should not happen for symmetric difference
        }

        // 检查 z5_output 是否依赖于 i5_cand
        // IsInTFO(i5_cand, z5_output) 等价于 i5_cand 在 z5_output 的扇入锥中
        if (!IsInTFO(i5_cand, z5_output, MAX_TFO_DEPTH)) {
            valid_I5.insert(i5_cand);
        }
    }

    if (valid_I5.empty()) return valid_I5;
    potential_I5 = valid_I5;
    
    // Safety check: ensure Z1 and Z2 are valid
    if (!g_lut6d_bit2States.count(Z1) || !g_lut6d_bit2States.count(Z2)) {
        return potential_I5;
    }
    
    const auto &vz1 = g_lut6d_bit2States[Z1];
    const auto &vz2 = g_lut6d_bit2States[Z2];
    
    // Safety check for vector sizes
    if (vz1.empty() || vz2.empty()) {
        return potential_I5;
    }
    
    size_t rows = std::min<size_t>({(size_t)NUM_PATTERNS, vz1.size(), vz2.size()});

    // 优化3：预先收集所有potential_I5的状态向量，避免重复查找
    dict<SigBit, const vector<State>*> i5_states;
    for (const SigBit &bit : potential_I5) {
        if (g_lut6d_bit2States.count(bit)) {
            i5_states[bit] = &g_lut6d_bit2States[bit];
        }
    }

    for (size_t i = 0; i < rows; i++) {
        // Bounds check before accessing
        if (i >= vz1.size() || i >= vz2.size()) break;
        
        State state_z1 = vz1[i];
        State state_z2 = vz2[i];
        if (state_z1 == State::Sx || state_z2 == State::Sx) continue;
        
        if (state_z1 != state_z2) {
            // 优化4：使用迭代器避免重复查找，直接在遍历中删除
            for (auto it = potential_I5.begin(); it != potential_I5.end(); ) {
                const SigBit &bit = *it;
                
                // 使用预先收集的状态向量
                auto state_it = i5_states.find(bit);
                if (state_it == i5_states.end()) {
                    ++it;
                    continue;
                }
                
                const vector<State> *vb = state_it->second;
                if (i >= vb->size()) {
                    ++it;
                    continue;
                }
                
                State state_bit = (*vb)[i];
                if (state_bit != State::S1) {
                    i5_states.erase(state_it);  // 同时从缓存中移除
                    it = potential_I5.erase(it);
                } else {
                    ++it;
                }
            }
        }
        
        // 早期退出：如果没有候选了
        if (potential_I5.empty()) break;
    }
    
    return potential_I5;
}

void GetPotentialCands(Cell *gate, NestedCandMap* &candidates, const pool<Cell*> &visited) {
    if (!gate) return;
    if (g_current_module && g_current_module->cell(gate->name) != gate) return; // stale gate
    if (!IsCombinationalGate(gate)) return;
    if (gate->get_bool_attribute(ID(_lut6d_zombie))) return;
    if (visited.count(gate)) return;

    candidates = new NestedCandMap();

    // 计算当前gate的MFFC - 使用快速版本
    // log("Computing MFFC for gate %s\n", gate->name.c_str());
    pool<Cell*> mffc;
	ComputeMFFCFast(gate, mffc);
    // log("  MFFC has %zu cells\n", mffc.size());

    // 将cut按照大小排序，优先尝试较大的cut
    // cell2cuts: dict<Cell *, dict<pool<SigBit>, pool<Cell *>>>
    // cell -> { cut_bits -> cone_cells }
    // 内层值 pool<Cell*>（cone_cells）：与该 cut 对应的“可行锥体”中的所有 Cell（从 cut_bits 这些叶子往上一直到 root 之间的内部节点）。
    // log("Sorting cuts for gate %s\n", gate->name.c_str());
    using CutKV = pair<const pool<SigBit>*, const pool<Cell*>*>;
    vector<CutKV> sorted_cuts;
    if (!cell2cuts.count(gate)) return; // gate might be removed or cuts not available
    sorted_cuts.reserve(cell2cuts[gate].size());
    for (auto &kv : cell2cuts[gate]) sorted_cuts.push_back({&kv.first, &kv.second});
    std::sort(sorted_cuts.begin(), sorted_cuts.end(), [](const CutKV &a, const CutKV &b){ return a.first->size() > b.first->size(); });
	// log("  Sorted cuts count: %zu\n", sorted_cuts.size());

    // 遍历这个gate的所有cut
    int cut_count = 0;
    for (auto cutKV : sorted_cuts) {
        // log(" Processing cut %d of size %zu\n", cut_count, cutKV.first->size());
        const pool<SigBit> &cut_bits = *cutKV.first;
        
        // Check if any bit in the cut is in removed_sigbit
        if (CutHasRemovedBits(cut_bits)) {
            // log("  Skipping cut with removed signal bits\n");
            continue; // Skip this cut
        }
        
        // 寻找该cut对应的所有candidate（有序）
        vector<pair<Cell*, vector<pair<pool<SigBit>, pool<SigBit>>>>> cand_list;
        static dict<Cell*, pool<Cell*>> mffc_cache; // cache across cuts of the same gate pass
        getPotentialCandidatesForCut(gate, cut_bits, cand_list, mffc_cache);
        // log("  Found %zu potential candidates for the cut\n", cand_list.size());

        // 填入返回结构：OuterEntry = pair<cut, vector<CandEntry>>
        // log("  Organizing candidates into NestedCandMap structure...\n");
        vector<NestedCandMap::CandEntry> cvec;
        cvec.reserve(cand_list.size());
        for (auto &cand_entry : cand_list) {
            Cell *cand = cand_entry.first;
            NestedCandMap::CandCuts ordered_cuts;
            ordered_cuts.reserve(cand_entry.second.size());
            for (auto &pc : cand_entry.second) {
                // Use move to avoid copying pools
                ordered_cuts.emplace_back(std::move(const_cast<pool<SigBit>&>(pc.first)), 
                                         std::move(const_cast<pool<SigBit>&>(pc.second)));
            }
            cvec.push_back({cand, std::move(ordered_cuts)});
        }
        if (!cvec.empty()) {
            candidates->data.push_back({cut_bits, std::move(cvec)});
        }
        // log("  Organized %zu candidates for the cut\n", cvec.size());
        cut_count++;
        if (cut_count >= MAX_CUTS_PER_GATE) {
            // log("  Reached maximum cuts per gate (%d), stopping further cut processing.\n", MAX_CUTS_PER_GATE);
            break;
        }
    }
    // log("  Gate %s has %zu cuts with potential candidates\n", gate->name.c_str(), candidates->size());
    // Removed detailed debug output for performance
    // 'visited' is passed by value; do not modify it here.
    // visited.insert(gate);

    return;
}

bool SATVerification(Cell *z5, Cell *z, const SigBit &I5, const pool<Cell*> &gate_fanin_cone, const pool<SigBit> &gate_leaves) {
    
    // Check SAT cache first
    SatCacheKey cache_key = {z5, z, I5};
    if (g_sat_cache.count(cache_key)) {
        return g_sat_cache[cache_key];
    }
    
    // 使用通用输出获取，避免硬编码 Y 端口导致空端口访问
    SigBit out1 = GetCellOutputCached(z5);
    SigBit out2 = GetCellOutputCached(z);
    if (out1 == SigBit() || out2 == SigBit()) {
        // log("  [SATVerification] Skip: invalid output bit for %s or %s.\n", z5->name.c_str(), z->name.c_str());
        g_sat_cache[cache_key] = false;
        return false;
    }
    
    // 快速检查：z5不能在I5的TFO中（等价于I5不能在z5的TFI中）
    // 这是最快的检查方式，因为：
    // 1. 只需要一次BFS遍历（从I5出发）
    // 2. 可以早停（一旦找到out1就返回）
    // 3. 避免了昂贵的TFI计算
    if (IsInTFO(I5, out1, MAX_TFO_DEPTH)) {
        // log("  [SATVerification] Failed: z5 output (%s) is in TFO of I5 (%s)\n", 
            // log_signal(out1), log_signal(I5));
        g_sat_cache[cache_key] = false;
        return false;
    }
    
    ezSatPtr ez;
	SatGen satgen(ez.get(), &sigmap, stringf("%s_vs_%s", z5->name.c_str(), z->name.c_str()));
	satgen.model_undef = false;

    SigBit Z, Z5;

    pool<Cell*> cand_fanin_cone;
    pool<SigBit> cand_leaves;
    // 使用有限深度和大小的GetFullFanInCone，防止在大型CPU设计中SAT爆炸
    // 参数: max_depth=20, max_cells=200 - 这会显著减少SAT变量数
    GetFullFanInCone(GetCellOutputCached(z), cand_fanin_cone, cand_leaves, 20, 200);
    
    // 统计cone大小，用于性能分析
    size_t total_cells = gate_fanin_cone.size() + cand_fanin_cone.size();
    if (total_cells > 100) {
        // log("  [SATVerification] Large cone detected: gate=%zu cells, cand=%zu cells, total=%zu\n",
        //     gate_fanin_cone.size(), cand_fanin_cone.size(), total_cells);
    }
    
    for (auto cell : gate_fanin_cone)
        satgen.importCell(cell);
    for (auto cell : cand_fanin_cone)
        satgen.importCell(cell);
    for (auto leaf : gate_leaves)
        satgen.importSigBit(leaf);
    for (auto leaf : cand_leaves)
        satgen.importSigBit(leaf);

    int constrained_var = satgen.importSigBit(I5);
    int out1_var = satgen.importSigBit(out1);
    int out2_var = satgen.importSigBit(out2);

    ez->assume(ez->NOT(constrained_var));
    ez->assume(ez->XOR(out1_var, out2_var));

    //  std::vector<int> solution;
    bool sat_result = ez->solve();
    if (sat_result) {
        // SAT: 找到了使 out1 != out2 的赋值
        // 意味着两个门的输出可以不相等
        // log("SAT: Gates are NOT equivalent under constraint\n");
        return false;  // 验证失败
    } 
    // else {
        // UNSAT: 在约束下，out1 和 out2 必然相等
        // log("UNSAT: Gates are equivalent under constraint\n");
        // return true;   // 验证通过
    // }

    ezSatPtr ez2;
    SatGen satgen2(ez2.get(), &sigmap, stringf("%s_vs_%s_dep", z5->name.c_str(), z->name.c_str()));
    satgen2.model_undef = false;

    // pool<Cell*> cand_fanin_cone;
    // pool<SigBit> cand_leaves;
    // GetFullFanInCone(GetCellOutputCached(z), cand_fanin_cone, cand_leaves);

    for (auto cell : gate_fanin_cone)
        satgen2.importCell(cell);
    for (auto cell : cand_fanin_cone)
        satgen2.importCell(cell);
    for (auto leaf : gate_leaves)
        satgen2.importSigBit(leaf);
    for (auto leaf : cand_leaves)
        satgen2.importSigBit(leaf);

    // To check the dependency, we need to model the logic of 'z5' twice
    // with different values for I5, while keeping all other inputs the same.
    dict<SigBit, int> subst_map2;
    pool<SigBit> all_leaves = gate_leaves;
    all_leaves.insert(cand_leaves.begin(), cand_leaves.end());

    for (auto leaf : all_leaves) {
        if (leaf != I5) {
            // Register the leaf signal with the SatGen instance and use its returned literal/int.
            // ez2->newLit() does not exist on the ezSAT object, so import via satgen2.
            subst_map2[leaf] = satgen2.importSigBit(leaf);
        }
    }

    // Model 1: z5_output when I5 is 0
    // Use the SatGen instance to import a constant sigbit for FALSE (S0)
    subst_map2[I5] = satgen2.importSigBit(State::S0);
    int out1_clone_i5_0 = satgen2.importSigBit(out1);
    
    // Model 2: z5_output when I5 is 1
    // Use the SatGen instance to import a constant sigbit for TRUE (S1)
    subst_map2[I5] = satgen2.importSigBit(State::S1);
    int out1_clone_i5_1 = satgen2.importSigBit(out1);

    // Assert that the output of z5 differs when I5 changes
    ez2->assume(ez2->XOR(out1_clone_i5_0, out1_clone_i5_1));

    bool z5_depends_on_i5 = ez2->solve();

    if (z5_depends_on_i5) {
        // Second check fails: Z5's logic depends on I5.
        // log("  SAT Check 2 failed: z5 (%s) output depends on I5 (%s)\n", z5->name.c_str(), log_signal(I5));
        g_sat_cache[cache_key] = false;
        return false;
    }

    // Both checks passed.
    // log("  SAT Check passed for z5=%s and z=%s\n", z5->name.c_str(), z->name.c_str());
    g_sat_cache[cache_key] = true;
    return true;
}

bool GetCandToMerge(Cell* gate, NestedCandMap* candidates, SigBit &I5_merge, Cell* &z_merge, Cell* &z5_merge, pool<SigBit> &candCut_merge, pool<SigBit> &gateCut_merge) {
    pool<Cell*> gate_fanin_cone;
    pool<SigBit> gate_leaves;
    // 使用有限深度fanin cone，加速SAT验证
    GetFullFanInCone(GetCellOutputCached(gate), gate_fanin_cone, gate_leaves, 20, 200);

    // 遍历候选单元（gate的不同cut）
    for (auto &cand : candidates->data) {
    const pool<SigBit> &cut = cand.first;

        // Check if gate cut has any removed bits (should already be filtered, but double-check)
        if (CutHasRemovedBits(cut)) continue; // Skip this gate cut

        // 遍历该 cut 下的所有候选 candidate cell
        for (auto &cand_entry : cand.second) {
            Cell *cand_cell = cand_entry.first;

            // 遍历该 candidate cell 的所有 cand_cut
            for (auto &cand_cut_pair : cand_entry.second) {
                const pool<SigBit> &cand_cut = cand_cut_pair.first;
                const pool<SigBit> &potential_I5 = cand_cut_pair.second;
                
                // Check if cand_cut has any removed bits
                if (CutHasRemovedBits(cand_cut)) continue; // Skip this cand_cut
                
                if (potential_I5.empty()) {
                    I5_merge = State::S1;
                    z_merge = gate;
                    z5_merge = cand_cell;
                    candCut_merge = cand_cut;
                    gateCut_merge = cut;
                    // log("  Successful merge of gate %s and candidate %s without I5\n",
                    //     gate->name.c_str(), cand_cell->name.c_str());
                    // log("  gate output: %s, candidate output: %s\n",
                    //     log_signal(GetCellOutputCached(gate)), log_signal(GetCellOutputCached(cand_cell)));
                    // log("  gate cut: ");
                    // for (auto &b : gateCut_merge) log("%s ", log_signal(b));
                    // log("\n");
                    // log("  candidate cut: ");
                    // for (auto &b : candCut_merge) log("%s ", log_signal(b));
                    // log("\n");
                    return true;
                } 

                // 遍历该 candidate cell 的所有 potential I5
                for (const SigBit &I5_cand : potential_I5) {
                        Cell* z5;
                        Cell* z;
                        bool i5_in_gate = gateCut_merge.count(I5_merge);
                        bool i5_in_cand = candCut_merge.count(I5_merge);
                        if (i5_in_gate && !i5_in_cand) {
                            z = gate;
                            z5 = cand_cell;
                        } else if (!i5_in_gate && i5_in_cand) {
                            z = cand_cell;
                            z5 = gate;
                        } else {
                            // 两边都不含或都含（极少见），统一按 gate 为 Z 处理，保证确定性
                            z = gate;
                            z5 = cand_cell;
                        }
                    bool sat_result = SATVerification(z5, z, I5_cand, gate_fanin_cone, gate_leaves);
                    if (sat_result) {
                        I5_merge = I5_cand;
                        z_merge = z;
                        z5_merge = z5;
                        gateCut_merge = cut;
                        candCut_merge = cand_cut;
                        // log("  Successful merge of gate %s and candidate %s with I5 %s\n",
                        //     gate->name.c_str(), cand_cell->name.c_str(), log_signal(I5_merge));
                        // log("  z: %s and z5: %s\n",
                        //     z_merge->name.c_str(), z5_merge->name.c_str());
                        // log("  gate output: %s, candidate output: %s\n",
                        //     log_signal(GetCellOutputCached(gate)), log_signal(GetCellOutputCached(cand_cell)));
                        // log("  gate cut: ");
                        // for (auto &b : gateCut_merge) log("%s ", log_signal(b));
                        // log("\n");
                        // log("  candidate cut: ");
                        // for (auto &b : candCut_merge) log("%s ", log_signal(b));
                        // log("\n");
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

// Helper: compute whether Z and Z5 have at least one common dependent input among non-I5 inputs
static bool HasCommonDependentInput(const vector<SigBit> &LUT6Dinput, const vector<bool> &LUT6Dinit)
{
    if (LUT6Dinput.size() != 6 || LUT6Dinit.size() != 64) return false;
    // For indices j=0..4 (non-I5), test dependency for Z (i5=1) and Z5 (i5=0)
    pool<int> depZ, depZ5;
    for (int j = 0; j < 5; j++) {
        // Skip constants as dependency candidates
        if (!LUT6Dinput[j].is_wire()) continue;
        bool z_dep = false, z5_dep = false;
        for (int idx = 0; idx < 64 && !(z_dep && z5_dep); idx++) {
            int i5 = (idx >> 5) & 1;
            int idx_flip = idx ^ (1 << j);
            if (i5 == 1) {
                if (LUT6Dinit[idx] != LUT6Dinit[idx_flip]) z_dep = true;
            } else {
                if (LUT6Dinit[idx] != LUT6Dinit[idx_flip]) {
                    z5_dep = true;
                    z_dep = true;
                }
            }
        }
        if (z_dep) depZ.insert(j);
        if (z5_dep) depZ5.insert(j);
    }
    // Check overlapping indices (representing overlapping non-I5 input wires)
    for (int j = 0; j < 5; j++) if (depZ.count(j) && depZ5.count(j)) return true;
    return false;
}

// Find all cells that will be removed after merging gate and cand_merge
// These cells are:
// 1. In the fan-in cone of gate and cand_merge
// 2. In their MFFC (Maximum Fanout-Free Cone)
// 3. Their TFO only leads to gate and/or cand_merge
static void FindCellsToBeRemoved(Cell *gate, Cell *cand_merge, const pool<SigBit> &lut6d_cut)
{
    pool<Cell*> cells_to_remove;
    // cells_to_remove.clear();
    if (!gate || !cand_merge) return;

    // 1. Get fan-in cones for both cells
    pool<Cell*> gate_fanin_cone, cand_fanin_cone;
    pool<SigBit> gate_leaves, cand_leaves;
    // 使用有限深度fanin cone - FindCellsToBeRemoved不需要完整cone
    GetFullFanInCone(GetCellOutputCached(gate), gate_fanin_cone, gate_leaves, 20, 200);
    GetFullFanInCone(GetCellOutputCached(cand_merge), cand_fanin_cone, cand_leaves, 20, 200);

    // 2. Union of both fan-in cones
    pool<Cell*> combined_fanin_cone = gate_fanin_cone;
    for (auto *c : cand_fanin_cone) combined_fanin_cone.insert(c);

    // 3. Build fanout map for cells in the combined cone
    dict<Cell*, pool<Cell*>> local_fanout;
    for (auto *c : combined_fanin_cone) {
        if (!c || !IsCombinationalGate(c)) continue;
        SigBit out = GetCellOutputCached(c);
        if (!bit2reader.count(out)) continue;
        
        for (Cell *reader : bit2reader[out]) {
            if (!reader || !IsCombinationalGate(reader)) continue;
            local_fanout[c].insert(reader);
        }
    }

    // 4. Find cells whose TFO only leads to gate and/or cand_merge
    // These are cells that:
    // - Are in the combined fan-in cone
    // - Are not in the lut6d_cut (cut inputs should be preserved)
    // - All their fanout cells (recursively) only lead to gate or cand_merge
    
    pool<Cell*> targets;
    targets.insert(gate);
    targets.insert(cand_merge);

    // Helper lambda to check if a cell's TFO only reaches target cells
    std::function<bool(Cell*, pool<Cell*>&)> tfo_only_to_targets = 
        [&](Cell *c, pool<Cell*> &visited_tfo) -> bool {
        if (!c || !IsCombinationalGate(c)) return true;
        if (targets.count(c)) return true;
        if (visited_tfo.count(c)) return true;
        visited_tfo.insert(c);

        // Check if this cell's output is in the cut (if so, it's a boundary)
        SigBit out = GetCellOutputCached(c);
        if (lut6d_cut.count(out)) return false;

        // Check all fanout cells
        if (!local_fanout.count(c)) return false; // No fanout means dead code (should be removed)
        
        for (Cell *fo : local_fanout[c]) {
            if (!tfo_only_to_targets(fo, visited_tfo)) return false;
        }
        return true;
    };

    // 5. Collect all cells that should be removed
    for (auto *c : combined_fanin_cone) {
        if (!c || !IsCombinationalGate(c)) continue;
        if (c == gate || c == cand_merge) continue;

        // Check if this cell's output is in the LUT6D cut
        SigBit out = GetCellOutputCached(c);
        if (lut6d_cut.count(out)) continue; // Cut inputs are boundaries, keep them

        // Check if TFO only leads to targets
        pool<Cell*> visited_tfo;
        if (tfo_only_to_targets(c, visited_tfo)) {
            cells_to_remove.insert(c);
            visited.insert(c);
            removed_sigbit.insert(GetCellOutputCached(c));
        }
    }

    // log("  Found %zu cells to be removed after merging %s and %s\n", 
    //     cells_to_remove.size(), gate->name.c_str(), cand_merge->name.c_str());
    // for (auto *c : cells_to_remove) {
    //     log("    - %s (type: %s)\n", c->name.c_str(), c->type.c_str());
    // }
}

// Cache for GetLUT6DINIT results to avoid recomputation
struct LUT6DCacheKey {
    Cell* z_merge;
    Cell* z5_merge;
    SigBit I5_merge;
    pool<SigBit> candCut;
    pool<SigBit> gateCut;
    
    bool operator==(const LUT6DCacheKey& other) const {
        return z_merge == other.z_merge && 
               z5_merge == other.z5_merge && 
               I5_merge == other.I5_merge &&
               candCut == other.candCut &&
               gateCut == other.gateCut;
    }
    
    template<typename Hasher>
    Hasher hash_into(Hasher h) const {
        h.hash32((uintptr_t)z_merge);
        h.hash32((uintptr_t)z5_merge);
        if (I5_merge.wire) {
            h.hash32((uintptr_t)I5_merge.wire);
            h.hash32(I5_merge.offset);
        }
        // Simple hash for pools - just XOR all elements
        for (const auto &bit : candCut) {
            if (bit.wire) {
                h.hash32((uintptr_t)bit.wire);
                h.hash32(bit.offset);
            }
        }
        for (const auto &bit : gateCut) {
            if (bit.wire) {
                h.hash32((uintptr_t)bit.wire);
                h.hash32(bit.offset);
            }
        }
        return h;
    }
};

static dict<LUT6DCacheKey, pair<bool, tuple<vector<SigBit>, vector<bool>, SigBit, SigBit>>> g_lut6d_init_cache;

bool GetLUT6DINIT(Cell *z_merge, Cell *z5_merge, SigBit I5_merge, const pool<SigBit> &candCut_merge, const pool<SigBit> &gateCut_merge, vector<SigBit> &LUT6Dinput, vector<bool> &LUT6Dinit, SigBit &Z, SigBit &Z5) {
    // Check cache first
    LUT6DCacheKey cache_key = {z_merge, z5_merge, I5_merge, candCut_merge, gateCut_merge};
    if (g_lut6d_init_cache.count(cache_key)) {
        auto &cached = g_lut6d_init_cache[cache_key];
        if (cached.first) {
            LUT6Dinput = std::get<0>(cached.second);
            LUT6Dinit = std::get<1>(cached.second);
            Z = std::get<2>(cached.second);
            Z5 = std::get<3>(cached.second);
        }
        return cached.first;
    }
    
    // log("Calculating INIT...\n");
    // 1) 使用选中的 cut 集合作为 LUT 的输入候选，并去除 I5；随后确定 Z/Z5 的归属
    pool<SigBit> input_set = gateCut_merge;
    for (auto &b : candCut_merge) input_set.insert(b);
    input_set.erase(I5_merge); // 非 I5 列表中不包含 I5

    // 选择 Z/Z5：优先依据 I5 是否出现在对应的 cut 中；若都不包含，默认 Z=gate, Z5=cand
    Z = GetCellOutputCached(z_merge);
    Z5 = GetCellOutputCached(z5_merge);

    // 2) 生成 I0..I5 顺序：先 5 个非 I5，再 I5 作为 I5（索引 5）
    // 优先选择 gate/cand cut 的交集，确保至少包含一个共同输入（占线层面）
    vector<SigBit> inter_vec, rest_vec;
    inter_vec.reserve(input_set.size());
    rest_vec.reserve(input_set.size());
    for (auto &b : input_set) {
        if (gateCut_merge.count(b) && candCut_merge.count(b)) inter_vec.push_back(b);
        else rest_vec.push_back(b);
    }
    auto by_name = [](const SigBit &a, const SigBit &b){ return strcmp(log_signal(a), log_signal(b)) < 0; };
    std::sort(inter_vec.begin(), inter_vec.end(), by_name);
    std::sort(rest_vec.begin(), rest_vec.end(), by_name);

    LUT6Dinput.clear();
    for (auto &b : inter_vec) {
        if (LUT6Dinput.size() >= 5) break;
        LUT6Dinput.push_back(b);
    }
    for (auto &b : rest_vec) {
        if (LUT6Dinput.size() >= 5) break;
        LUT6Dinput.push_back(b);
    }
    // 若非 I5 不足 5 个，用 0 常量补齐（极端情况）
    while (LUT6Dinput.size() < 5) LUT6Dinput.push_back(State::S0);
    // 追加 I5
    LUT6Dinput.push_back(I5_merge);

    // 3) 遍历 64 种输入组合，计算 INIT：
    // 低 32 位（I5=0）输出对应 Z5，高 32 位（I5=1）输出对应 Z。
    LUT6Dinit.assign(64, false);
    for (int idx = 0; idx < 64; idx++) {
        dict<SigBit, State> bit_map;
        for (int j = 0; j < 6; j++) {
            SigBit inb = LUT6Dinput[j];
            State s = ((idx >> j) & 1) ? State::S1 : State::S0;
            bit_map[inb] = s;
        }
        // log("  Bitmap contains: ");
        // for (auto bit : bit_map) {
        //     log("%s, ", log_signal(bit.first));
        // }
        // log("\n");
        // log("  When the bitmap belike: ");
        // for (int j = 0; j < 6; j++) {
        //     SigBit inb = LUT6Dinput[j];
        //     State s = bit_map[inb];
        //     char sc = (s == State::S0) ? '0' : (s == State::S1) ? '1' : 'x';
        //     log("%c ", sc);
        // }
        // log("\n");
        // 计算两路输出
        // log("\tGetting z: %s\n", log_signal(Z));
        State z_main = StateEvalSoft(bit_map, Z);
        if (z_main == State::Sx) {
            // log("Oh no z5 is Sx!\n");
            return false;
        }
        // log("   z = %d\n", z_main == State::S1);
        // log("\tGetting z5: %s\n", log_signal(Z5));
        State z5_aux = StateEvalSoft(bit_map, Z5);
        // log("   z5 = %d\n", z5_aux == State::S1);
        bool i5_is_one = ((idx >> 5) & 1) != 0;
        bool out_bit = false;
        if (i5_is_one) out_bit = (z_main == State::S1);
        else out_bit = (z5_aux == State::S1);
        LUT6Dinit[idx] = out_bit;
    }

    // log("  Generated LUT6D INIT: \n  ");
    // for (int i = 63; i >= 0; i--) {
    //     log("%d", LUT6Dinit[i] ? 1 : 0);
    //     if (i % 8 == 0) log("\n  ");
    // }
    // log("\n");

    // 合并前的自检：Z/Z5 是否存在共同依赖的非 I5 输入
    bool ok = HasCommonDependentInput(LUT6Dinput, LUT6Dinit);
    if (!ok) {
        // log("  [skip] LUT6D candidate rejected: Z and Z5 have no common dependent input.\n");
        g_lut6d_init_cache[cache_key] = {false, {}};
        return false;
    }
    
    // Store successful result in cache
    g_lut6d_init_cache[cache_key] = {true, std::make_tuple(LUT6Dinput, LUT6Dinit, Z, Z5)};
    return true;
}

// Track removed cells to avoid accessing deleted pointers
static pool<Cell*> g_removed_cells;

// Utility: run a pass and report which cells in the current module were removed
static void RunPassWithRemovedCellLog(Module *module, const std::string &cmd)
{
    if (module == nullptr || module->design == nullptr) {
        Pass::call(nullptr, cmd);
        return;
    }

    pool<IdString> before;
    for (auto cell : module->cells()) before.insert(cell->name);

    Pass::call(module->design, cmd);

    pool<IdString> after;
    for (auto cell : module->cells()) after.insert(cell->name);

    std::vector<IdString> removed;
    removed.reserve(before.size());
    for (auto &n : before) if (!after.count(n)) removed.push_back(n);

    // log("  [cleanup] %s removed %zu cells in module %s\n", cmd.c_str(), removed.size(), log_id(module));
    // for (auto &n : removed) {
    //     log("    - %s\n", log_id(n));
    //     g_removed_cell_names.insert(n);
    // }
}

Cell* addLUT6D(Module *module, const vector<SigBit> &inputs, const vector<bool> &init, const SigBit &Z, const SigBit &Z5) {
    log_assert(inputs.size() == 6);
    log_assert(init.size() == 64);

    // Create cell
    Cell *lut6d = module->addCell(NEW_ID, ID(GTP_LUT6D));

    // Connect inputs I0..I5
    lut6d->setPort(ID(I0), SigSpec(inputs[0]));
    lut6d->setPort(ID(I1), SigSpec(inputs[1]));
    lut6d->setPort(ID(I2), SigSpec(inputs[2]));
    lut6d->setPort(ID(I3), SigSpec(inputs[3]));
    lut6d->setPort(ID(I4), SigSpec(inputs[4]));
    lut6d->setPort(ID(I5), SigSpec(inputs[5]));

    // Build INIT as a 64-bit Const (bit i corresponds to address formed by I0..I5 with I5 as MSB)
    std::vector<RTLIL::State> init_states;
    init_states.reserve(64);
    for (int i = 0; i < 64; i++)
        init_states.push_back(init[i] ? RTLIL::State::S1 : RTLIL::State::S0);
    lut6d->setParam(ID(INIT), RTLIL::Const(init_states));

    // Connect outputs
    lut6d->setPort(ID(Z), SigSpec(Z));
    lut6d->setPort(ID(Z5), SigSpec(Z5));

    // Optional: textual INIT for log
    string init_str;
    init_str.reserve(64);
    for (int i = 63; i >= 0; i--) init_str += init[i] ? '1' : '0';
    // log("  Created LUT6D cell %s with INIT %s\n", lut6d->name.c_str(), init_str.c_str());
    return lut6d;
}

void PerformMerge(Module *module, const MergeInfo &merge_info) {
    log_assert(merge_info.gate != nullptr && merge_info.cand_merge != nullptr);

    // 1) Create the LUT6D driving original nets Z and Z5
    Cell *lut6d = addLUT6D(module, merge_info.LUT6Dinput, merge_info.LUT6Dinit, merge_info.Z, merge_info.Z5);

    // 2) Update driver mapping for Z/Z5 to the new LUT6D
    bit2driver[merge_info.Z] = lut6d;
    bit2driver[merge_info.Z5] = lut6d;
    
    // 3) Update cell2bits mapping for the new LUT6D (dual output)
    // Format: {Z, Z5, I0, I1, I2, I3, I4, I5}
    vector<SigBit> lut6d_bits;
    lut6d_bits.push_back(merge_info.Z);
    lut6d_bits.push_back(merge_info.Z5);
    for (auto &inp : merge_info.LUT6Dinput) {
        lut6d_bits.push_back(inp);
    }
    cell2bits[lut6d] = lut6d_bits;

    // 4) Clean up bit2reader: remove gate/cand_merge from ALL their input signals' reader lists
    // Get all input signals from gate and cand_merge before modifying them
    pool<SigBit> gate_inputs, cand_inputs;
    if (cell2bits.count(merge_info.gate)) {
        auto &gate_bits = cell2bits[merge_info.gate];
        for (size_t i = 1; i < gate_bits.size(); i++) gate_inputs.insert(gate_bits[i]);
    }
    if (cell2bits.count(merge_info.cand_merge)) {
        auto &cand_bits = cell2bits[merge_info.cand_merge];
        for (size_t i = 1; i < cand_bits.size(); i++) cand_inputs.insert(cand_bits[i]);
    }
    
    // Remove gate and cand_merge from bit2reader for all their inputs
    for (auto &inp : gate_inputs) {
        if (bit2reader.count(inp)) {
            auto &readers = bit2reader[inp];
            readers.erase(std::remove(readers.begin(), readers.end(), merge_info.gate), readers.end());
            if (readers.empty()) bit2reader.erase(inp);
        }
    }
    for (auto &inp : cand_inputs) {
        if (bit2reader.count(inp)) {
            auto &readers = bit2reader[inp];
            readers.erase(std::remove(readers.begin(), readers.end(), merge_info.cand_merge), readers.end());
            if (readers.empty()) bit2reader.erase(inp);
        }
    }
    
    // Add lut6d as reader for its inputs
    for (auto &inp : merge_info.LUT6Dinput) {
        if (!bit2reader.count(inp)) bit2reader[inp] = vector<Cell*>();
        auto &readers = bit2reader[inp];
        if (std::find(readers.begin(), readers.end(), lut6d) == readers.end()) {
            readers.push_back(lut6d);
        }
    }

    // 5) Remove old cells' outputs from bit2driver (before modifying/deleting cells)
    SigBit gate_out = GetCellOutputCached(merge_info.gate);
    SigBit cand_out = GetCellOutputCached(merge_info.cand_merge);
    
    if (bit2driver.count(gate_out) && bit2driver[gate_out] == merge_info.gate) {
        bit2driver.erase(gate_out);
    }
    if (bit2driver.count(cand_out) && bit2driver[cand_out] == merge_info.cand_merge) {
        bit2driver.erase(cand_out);
    }

    // 6) Remove old cells from bit2reader for their outputs
    if (bit2reader.count(gate_out)) {
        auto &readers = bit2reader[gate_out];
        readers.erase(std::remove(readers.begin(), readers.end(), merge_info.gate), readers.end());
        if (readers.empty()) bit2reader.erase(gate_out);
    }
    if (bit2reader.count(cand_out)) {
        auto &readers = bit2reader[cand_out];
        readers.erase(std::remove(readers.begin(), readers.end(), merge_info.cand_merge), readers.end());
        if (readers.empty()) bit2reader.erase(cand_out);
    }

    // 7) Remove old cells from cell2bits BEFORE modifying ports
    cell2bits.erase(merge_info.gate);
    cell2bits.erase(merge_info.cand_merge);

    // 8) Remove old cells from g_orig_fanout and g_cell_level
    if (g_orig_fanout.count(merge_info.gate)) g_orig_fanout.erase(merge_info.gate);
    if (g_orig_fanout.count(merge_info.cand_merge)) g_orig_fanout.erase(merge_info.cand_merge);
    if (g_cell_level.count(merge_info.gate)) g_cell_level.erase(merge_info.gate);
    if (g_cell_level.count(merge_info.cand_merge)) g_cell_level.erase(merge_info.cand_merge);

    // 9) Remove old cells from cell2cuts if present
    if (cell2cuts.count(merge_info.gate)) cell2cuts.erase(merge_info.gate);
    if (cell2cuts.count(merge_info.cand_merge)) cell2cuts.erase(merge_info.cand_merge);

    // 10) Mark cells as zombies to avoid accessing them before cleanup
    merge_info.gate->set_bool_attribute(ID(_lut6d_zombie), true);
    merge_info.cand_merge->set_bool_attribute(ID(_lut6d_zombie), true);
    
    // 11) Disconnect and remove old cells from module
    // Disconnect all ports first to make cells unreachable
    merge_info.gate->unsetPort(ID::A);
    merge_info.gate->unsetPort(ID::B);
    merge_info.gate->unsetPort(ID::S);
    merge_info.gate->unsetPort(ID::Y);
    merge_info.cand_merge->unsetPort(ID::A);
    merge_info.cand_merge->unsetPort(ID::B);
    merge_info.cand_merge->unsetPort(ID::S);
    merge_info.cand_merge->unsetPort(ID::Y);
    
    // Remove cells from module (will be cleaned up by opt_clean)
    IdString gatename = merge_info.gate->name;
    IdString candname = merge_info.cand_merge->name;
    module->remove(merge_info.gate);
    // log("   Gate %s removed!\n", log_id(gatename));
    module->remove(merge_info.cand_merge);
    // log("   Cand %s removed!\n", log_id(candname));
    // log("  Merged %s and %s into %s (old cells removed).\n",
        // merge_info.gate->name.c_str(), merge_info.cand_merge->name.c_str(), lut6d->name.c_str());
}

void LUT6DMapper(Module *module) {

    log("Perform lut6d mapping\n");
    CheckCellWidth(module);
    g_current_module = module;

    // Clear global state and caches
    removed_sigbit.clear();
    visited.clear();
    g_tfo_cache.clear();
    g_fanin_cone_cache.clear();
    g_sat_cache.clear();
    g_cell_output_cache.clear();

    // Pre-populate cell output cache for all cells - CRITICAL OPTIMIZATION
    // This eliminates 2.5 billion lookups in cell2bits
    log("  Building cell output cache...\n");
    for (auto cell : module->cells()) {
        if (!IsCombinationalGate(cell)) continue;
        SigBit out = GetCellOutput(cell);  // Use original GetCellOutput here
        if (out != SigBit()) {
            g_cell_output_cache[cell] = out;
        }
    }
    log("  Cached %zu cell outputs\n", g_cell_output_cache.size());

    // Count combinational gates for dynamic parameter adjustment
    int num_comb_gates = 0;
    for (auto cell : module->cells()) {
        if (IsCombinationalGate(cell)) num_comb_gates++;
    }
    
    // Adjust parameters based on netlist size
    AdjustParametersBasedOnNetlistSize(num_comb_gates);

    // 初始化分析结构
    BuildOrigFanout();
    BuildCellLevels(module);
    GenerateCuts(module);
    pool<SigBit> prime_inputs, prime_outputs;
    GetPrimeInputOuput(module, prime_inputs, prime_outputs);
    log("  Found %zu prime inputs and %zu prime outputs\n", prime_inputs.size(), prime_outputs.size());
    // log("    Primary inputs: ");
    // for (auto &b : prime_inputs) {
    //     log("%s ", log_signal(b));
    // }
    // log("\n");
    // log("    Primary outputs: ");
    // for (auto &b : prime_outputs) {
    //     log("%s ", log_signal(b));
    // }
    // log("\n");
    log(" --- Generate Bit States --- \n");
    GenerateBitState(module, prime_inputs, false);
    log("  Generated bit states for %zu signals\n", g_lut6d_bit2States.size());

    // 一次 topo 排序，按名称遍历；批量清理间隔
    vector<Cell*> topo_cells;
    vector<Cell*> re_topo_cells;
    log(" topo sorting gates...\n");
    GetTopoSortedGates(module, topo_cells);

    // keep original topo_cells and create a reversed copy for reverse-topo order
    re_topo_cells = topo_cells;
    std::reverse(re_topo_cells.begin(), re_topo_cells.end());
    vector<IdString> re_gate_names; 
    re_gate_names.reserve(re_topo_cells.size());
    for (auto *c : re_topo_cells) if (c && IsCombinationalGate(c)) re_gate_names.push_back(c->name);
    log(" topo sorting done. %zu gates to process.\n", re_gate_names.size());

    int i = 1;

    static dict<Cell*, MergeInfo> moduleMergeInfo;

    for (auto gate : re_topo_cells) {
        if (!gate) continue; // 可能被清理
        if (!IsCombinationalGate(gate)) continue;
        if (gate->get_bool_attribute(ID(_lut6d_zombie))) continue;
        if (visited.count(gate)) continue;

        log("%d. Processing gate %s\n", i, gate->name.c_str());
        NestedCandMap* candidates = nullptr;
        GetPotentialCands(gate, candidates, visited);
        if (!candidates || candidates->empty()) {
            if (candidates) delete candidates;
            i++;
            continue;
        }

        SigBit I5_merge = SigBit();
        Cell *z_merge;
        Cell *z5_merge;
        pool<SigBit> candCut_merge;
        pool<SigBit> gateCut_merge;
        if (!GetCandToMerge(gate, candidates, I5_merge, z_merge, z5_merge, candCut_merge, gateCut_merge)) {
            delete candidates;
            i++;
            continue;
        }
        // log("I5: %s\n", log_signal(I5_merge));
        // log("Z = %s\n", log_id(z_merge->name));
        // log("Z5 = %s\n", log_id(z5_merge->name));

        vector<SigBit> LUT6Dinput;
        vector<bool> LUT6Dinit;
        SigBit Z = SigBit();
        SigBit Z5 = SigBit();
        // log("I5: %s\n", log_signal(I5_merge));
        if (!GetLUT6DINIT(z_merge, z5_merge, I5_merge, candCut_merge, gateCut_merge, LUT6Dinput, LUT6Dinit, Z, Z5)) {
            // 无共同依赖，放弃当前 cand
            delete candidates;
            i++;
            continue;
        }
        
        // Find cells that will be removed after merge
        pool<SigBit> lut6d_cut;
        for (auto &b : gateCut_merge) lut6d_cut.insert(b);
        for (auto &b : candCut_merge) lut6d_cut.insert(b);
        
        // pool<Cell*> cells_to_remove;
        FindCellsToBeRemoved(z_merge, z5_merge, lut6d_cut);
        
        // Create MergeInfo structure with cells to remove
        MergeInfo merge_info(z_merge, z5_merge, LUT6Dinput, LUT6Dinit, Z, Z5);
        
        // Store merge_info for both gate and cand_merge (they share the same merge info)
        moduleMergeInfo[gate] = merge_info;
        // moduleMergeInfo[cand_merge] = merge_info;
        
        delete candidates;
        visited.insert(z_merge);
        visited.insert(z5_merge);

        i++;
        
        // 内存管理：每处理50个gate，清理不再需要的数据（提高频率）
        if (i % 50 == 0) {
            // 清理已处理gate的bit states来节省内存
            for (auto v : visited) {
                SigBit out = GetCellOutputCached(v);
                if (out != SigBit() && g_lut6d_bit2States.count(out)) {
                    g_lut6d_bit2States.erase(out);
                }
            }
            
            // 定期进度报告（每1000个gate）
            if (i % 1000 == 0) {
                log("  Progress: %d/%zu gates processed (%.1f%%)\n", 
                    i, re_gate_names.size(), 
                    100.0 * i / re_gate_names.size());
            }
            
            // 定期清理缓存以防止内存膨胀 - 增大缓存限制以提高性能
            if (g_tfo_cache.size() > 100000) {
                // 清理一半旧的缓存项而不是全部清空
                int to_remove = g_tfo_cache.size() / 2;
                auto it = g_tfo_cache.begin();
                while (to_remove-- > 0 && it != g_tfo_cache.end()) {
                    it = g_tfo_cache.erase(it);
                }
            }
            if (g_fanin_cone_cache.size() > 50000) {
                // 清理一半旧的缓存项
                int to_remove = g_fanin_cone_cache.size() / 2;
                auto it = g_fanin_cone_cache.begin();
                while (to_remove-- > 0 && it != g_fanin_cone_cache.end()) {
                    it = g_fanin_cone_cache.erase(it);
                }
            }
            if (g_sat_cache.size() > 20000) {
                // 清理一半旧的缓存项
                int to_remove = g_sat_cache.size() / 2;
                auto it = g_sat_cache.begin();
                while (to_remove-- > 0 && it != g_sat_cache.end()) {
                    it = g_sat_cache.erase(it);
                }
            }
        }
    }

    // pool<Cell *> merged;
    for (auto name : re_gate_names) {
        Cell *cell = module->cell(name);
        // log(" Trying to process: %s\n", log_id(name));
        if (!cell) continue;
        if (!IsCombinationalGate(cell)) continue;
        // if (merged.count(cell)) continue;
        if (!moduleMergeInfo.count(cell)) continue; // 如果cell不在moduleMergeInfo中，跳过

        const MergeInfo &mi = moduleMergeInfo[cell];
        // log("   Merging gate %s and cand %s\n", log_id(cell->name), log_id(mi.cand_merge->name));
        IdString cellname = mi.gate->name;
        PerformMerge(module, mi);
        // merged.insert(mi.gate);
        // merged.insert(mi.cand_merge);
        log("   Cell %s merged\n", log_id(cellname));
    }
    g_current_module = nullptr;

    // 清理全局数据结构
    g_lut6d_bit2States.clear();
    moduleMergeInfo.clear();

    // 最终清理
    log("  Running final cleanup pass...\n");
    RunPassWithRemovedCellLog(module, "opt -purge");

    return;
}

YOSYS_NAMESPACE_END
