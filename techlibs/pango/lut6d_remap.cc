/*
 *  yosys -- Yosys Open SYnthesis Suite
 *
 *  Copyright (C) 2025  Shenzhen Pango Microsystems Co., Ltd. <marketing@pangomicro.com>
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

/*
Area flow based lut-mapper.
Implement of "Heuristics for Area Minimization in LUT-Based FPGA Technology Mapping".
*/

#include "kernel/celltypes.h"
#include "kernel/consteval.h"
#include "kernel/modtools.h"
#include "kernel/sigtools.h"
#include "kernel/yosys.h"
// #include <algorithm>
#include <queue>
#include <ranges>
#include <string.h>
#include "synth_pango.h"

USING_YOSYS_NAMESPACE
using namespace std;
YOSYS_NAMESPACE_BEGIN

// -----------------------
// global variables delare here

// ---- add code here ----
dict<Cell *, size_t> lut2levels;
dict<size_t, pool<Cell *>> level2luts;
// -----------------------



// ---- add code here ----
// 需要用Resubstitution的算法！！！！！
// 已有：
// 1. Primary Input和 Output
// 2. 所有的门（LUT）按照topological order进行排序
// TODO：
// 1. 计算重汇聚割 √
// 2. 收集除数节点

void GetTopoSortedLUTs(Module *module, vector<Cell *> &luts)
{
    luts.clear();
    lut2levels.clear();
	level2luts.clear();
    
    pool<SigBit> prime_inputs;
    pool<SigBit> prime_outputs;
    GetPrimeInputOuput(module, prime_inputs, prime_outputs);
    
    dict<Cell *, size_t> indegree;
    dict<Cell *, size_t> max_input_level; // 记录每个LUT输入的最大层级
    pool<Cell *> visited;
    queue<Cell *> zero_indegree_nodes;
    
    // 计算每个LUT的入度，并初始化层级
    for (auto &cell_iter : module->cells_) {
        Cell *cell = cell_iter.second;
        if (!IsCombinationalCell(cell)) {
            continue;
        }
        
        log_assert(cell2bits.count(cell));
        auto &bits = cell2bits[cell];
        
        // 根据LUT类型设置初始入度
        if (IsGTP_LUT6D(cell)) {
            indegree[cell] = bits.size() - 2; // 两个输出位：Z和Z5
        } else {
            indegree[cell] = bits.size() - 1; // 一个输出位：Z
        }
        
        // 初始化最大输入层级为0
        max_input_level[cell] = 0;
        
        // 获取输入信号
        vector<SigBit> inputs;
        GetCellInputsVector(cell, inputs);
        
        // 检查每个输入的驱动源
        for (auto b : inputs) {
            Cell *drv = bit2driver[b];
            if (!drv || !IsCombinationalCell(drv)) {
                indegree[cell]--; // 如果驱动源不是组合逻辑单元，入度减1
            }
        }
        
        // 入度为0的节点加入队列（这些是第0层的LUT）
        if (indegree[cell] == 0) {
            zero_indegree_nodes.push(cell);
            lut2levels[cell] = 0; // 第0层
			level2luts.operator[](0).insert(cell);
        }
    }
    
    
    // BFS拓扑排序，同时计算层级
    while (!zero_indegree_nodes.empty()) {
        Cell *current = zero_indegree_nodes.front();
        zero_indegree_nodes.pop();
        
        if (!IsCombinationalCell(current) || visited.count(current)) {
            continue;
        }
        
        luts.push_back(current);
        visited.insert(current);
        
        // 记录当前LUT到对应层级
        size_t current_level = lut2levels[current];
        // 获取当前LUT的所有读者
        pool<Cell *> all_readers;
        
        // 处理Z输出
        SigBit outbit = GetCellOutput(current);
        if (bit2reader.count(outbit)) {
            auto readers = bit2reader[outbit];
            for (Cell *reader : readers) {
                all_readers.insert(reader);
            }
        }
        
        // 如果是LUT6D，还需要处理Z5输出
        if (IsGTP_LUT6D(current)) {
            auto &bits = cell2bits[current];
            if (bits.size() >= 2) {
                SigBit outbit_z5 = bits[1]; // Z5输出
                if (bit2reader.count(outbit_z5)) {
                    auto readers_z5 = bit2reader[outbit_z5];
                    for (Cell *reader : readers_z5) {
                        all_readers.insert(reader);
                    }
                }
            }
        }
        
        for (Cell *neighbor : all_readers) {
            if (!IsCombinationalCell(neighbor)) {
                continue; // 只处理组合逻辑单元
            }
            
            if (visited.count(neighbor)) {
                log_error("toposort found loop, LUT %s reader %s\n", 
                    current->name.c_str(), neighbor->name.c_str());
                continue;
            }
            
            // 更新邻居的最大输入层级
            if (lut2levels.count(current)) {
                max_input_level[neighbor] = std::max(max_input_level[neighbor], 
                                                    lut2levels[current] + 1);
            }
            
            // 更新邻居的入度
            if (--indegree[neighbor] == 0) {
                // 当所有输入都准备好时，设置该LUT的层级
                lut2levels[neighbor] = max_input_level[neighbor];
                zero_indegree_nodes.push(neighbor);
                
				// 将LUT添加到对应层级的集合中
				size_t level = max_input_level[neighbor];
				level2luts[level].insert(neighbor);
            }
        }
	}
    // 输出层级信息（调试用）
	log("LUT Level Information:\n");
	for (auto &level_pair : level2luts) {
		size_t level = level_pair.first;
		auto &level_luts = level_pair.second;
		log("  Level %lu: %lu LUTs\n", level, level_luts.size());
		for (auto lut : level_luts) {
			log("    - %s\n", lut->name.c_str());
		}
	}

    log("Total LUTs processed: %lu\n", luts.size());
}

RTLIL::Cell *addLut6d(Module *module, const string &name, const vector<SigBit> &inputs, const vector<SigBit> &outputs, Const INIT=Const(0, 64))
{
    // 验证参数
    log_assert(outputs.size() == 2); // LUT6D必须有两个输出：Z和Z5
    log_assert(inputs.size() <= 6);  // LUT6D最多6个输入
    
    // 创建新的LUT6D单元
    string lutName = name.empty() ? NEW_ID.str() : name;
    Cell *lut = module->addCell(IdString(lutName), ID(GTP_LUT6D));
	log("Successfully added LUT6D: %s\n", lutName.c_str());
    // 设置输出端口
    lut->setPort(ID(Z), outputs[0]);   // 主输出Z
    lut->setPort(ID(Z5), outputs[1]);  // 第二输出Z5
    log("Adding LUT6D: %s with %lu inputs\n", lutName.c_str(), inputs.size());
    // 设置输入端口（I0到I5）
    IdString input_ports[] = {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)};
    for (size_t i = 0; i < 6; i++) {
        if (i < inputs.size()) {
            lut->setPort(IdString(input_ports[i]), inputs[i]);
        } else {
            lut->setPort(IdString(input_ports[i]), State::S0);
        }
    }
	// log("LUT6D %s inputs: ", lutName.c_str());
	// for (size_t i = 0; i < inputs.size(); i++) {
	// 	log("%s ", log_signal(inputs[i]));
	// }
    // 初始化LUT表为全0（需要后续通过真值表计算设置）
    // INIT参数是64位，对应6输入LUT的所有可能组合
    lut->setParam(ID(INIT), INIT);
    log("LUT6D %s INIT: %s\n", lutName.c_str(), INIT.as_string().c_str());
    // 将新创建的LUT添加到相关数据结构中
    vector<SigBit> all_bits;
    all_bits.push_back(outputs[0]); // Z输出
    all_bits.push_back(outputs[1]); // Z5输出
    for (const auto& input : inputs) {
        all_bits.push_back(input);
    }
    cell2bits[lut] = all_bits;
    log("LUT6D %s cell2bits: ", lutName.c_str());
	for (const auto& bit : all_bits) {
		log("%s ", log_signal(bit));
	}
	log("\n");
    // 更新bit2driver映射
    bit2driver[outputs[0]] = lut;
    bit2driver[outputs[1]] = lut;
    
    return lut;
}

void GetCellOutputsVector(Cell *cell, vector<SigBit> &outputs)
{
	log_assert(cell && cell2bits.count(cell));
	auto &bits = cell2bits[cell];
	outputs.clear();
	int offset = 0;
	if (IsGTP_LUT6D(cell)) {
		offset = 2; // LUT6D有两个输出
	} else {
		offset = 1; // 其他LUT只有一个输出
	}
	for (int i = 0; i < offset; i++) {
		outputs.push_back(bits[i]);
	}
}

pair<Const, vector<SigBit>> checkCombinable(Cell *lut1, Cell *lut2, SigBit I5) {
	// 获取两个LUT的输入和INIT值
	vector<SigBit> lut1_inputs, lut2_inputs;
	GetCellInputsVector(lut1, lut1_inputs);
	GetCellInputsVector(lut2, lut2_inputs);
	// 按照I0~I5的顺序对输入进行排序
	std::sort(lut1_inputs.begin(), lut1_inputs.end(), [&](const SigBit &a, const SigBit &b) {
		IdString ports[] = {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)};
		int idx_a = -1, idx_b = -1;
		for (int i = 0; i < 6; i++) {
			if (lut1->hasPort(ports[i]) && lut1->getPort(ports[i]) == a) {
				idx_a = i;
			}
			if (lut1->hasPort(ports[i]) && lut1->getPort(ports[i]) == b) {
				idx_b = i;
			}
		}
		return idx_a < idx_b;
	});

	std::sort(lut2_inputs.begin(), lut2_inputs.end(), [&](const SigBit &a, const SigBit &b) {
		IdString ports[] = {ID(I0), ID(I1), ID(I2), ID(I3), ID(I4), ID(I5)};
		int idx_a = -1, idx_b = -1;
		for (int i = 0; i < 6; i++) {
			if (lut2->hasPort(ports[i]) && lut2->getPort(ports[i]) == a) {
				idx_a = i;
			}
			if (lut2->hasPort(ports[i]) && lut2->getPort(ports[i]) == b) {
				idx_b = i;
			}
		}
		return idx_a < idx_b;
	});
	vector<SigBit> merged_inputs;

	if (!lut1->hasParam(ID::INIT) || !lut2->hasParam(ID::INIT)) {
		return make_pair(Const(0, 64), vector<SigBit>());
	}

	Const lut1_init = lut1->getParam(ID::INIT);
	Const lut2_init = lut2->getParam(ID::INIT);

	// 构建合并后的输入列表
	pool<SigBit> all_inputs_set;
	for (auto input : lut1_inputs) {
		all_inputs_set.insert(input);
	}
	for (auto input : lut2_inputs) {
		all_inputs_set.insert(input);
	}

	// 确保I5在输入列表中
	all_inputs_set.insert(I5);

	if (all_inputs_set.size() > 6) {
		return {Const(0, 64), vector<SigBit>()};
	}

	// 将输入集合转换为有序向量，I5放在最后（作为第6个输入）
	for (auto input : all_inputs_set) {
		if (input != I5) {
			merged_inputs.push_back(input);
		}
	}
	merged_inputs.push_back(I5); // I5作为最高位输入

	// 补齐到6个输入
	while (merged_inputs.size() < 6) {
		merged_inputs.insert(merged_inputs.end() - 1, SigBit()); // 在I5前插入无效输入
	}

	// 优化：预先构建merged_inputs的反向索引，避免重复的std::find调用
	dict<SigBit, int> merged_input_index;
	for (int i = 0; i < merged_inputs.size(); i++) {
		if (merged_inputs[i].wire != nullptr || merged_inputs[i].data != State::S0) {
			merged_input_index[merged_inputs[i]] = i;
		}
	}

	// 创建输入映射：从合并后的输入到原LUT输入的映射
	dict<SigBit, int> lut1_input_map, lut2_input_map;
	for (int i = 0; i < lut1_inputs.size(); i++) {
		if (merged_input_index.count(lut1_inputs[i])) {
			lut1_input_map[lut1_inputs[i]] = merged_input_index[lut1_inputs[i]];
		}
	}
	for (int i = 0; i < lut2_inputs.size(); i++) {
		if (merged_input_index.count(lut2_inputs[i])) {
			lut2_input_map[lut2_inputs[i]] = merged_input_index[lut2_inputs[i]];
		}
	}

	// 计算LUT6D的INIT值
	vector<bool> lut6d_init_bits(64, false);
	bool compatible = true;

	// 遍历所有64种输入组合
	for (int addr = 0; addr < 64; addr++) {
		// 提取I5的值（最高位）
		bool i5_val = (addr >> 5) & 1;

		// 计算lut1和lut2在当前输入组合下的输出
		bool lut1_output = false, lut2_output = false;

		// 计算lut1的输出
		int lut1_addr = 0;
		for (int i = 0; i < lut1_inputs.size(); i++) {
			if (lut1_input_map.count(lut1_inputs[i])) {
				int merged_idx = lut1_input_map[lut1_inputs[i]];
				bool bit_val = (addr >> merged_idx) & 1;
				lut1_addr |= (bit_val << i);
			}
		}
		if (lut1_addr < lut1_init.size()) {
			lut1_output = lut1_init.bits()[lut1_addr] == State::S1;
		}

		// 计算lut2的输出
		int lut2_addr = 0;
		for (int i = 0; i < lut2_inputs.size(); i++) {
			if (lut2_input_map.count(lut2_inputs[i])) {
				int merged_idx = lut2_input_map[lut2_inputs[i]];
				bool bit_val = (addr >> merged_idx) & 1;
				lut2_addr |= (bit_val << i);
			}
		}
		if (lut2_addr < lut2_init.size()) {
			lut2_output = lut2_init.bits()[lut2_addr] == State::S1;
		}

		// 检查兼容性规则
		// 特殊处理：如果I5是常量1，则跳过某些检查
		bool is_i5_const_one = (I5.wire == nullptr && I5.data == State::S1);
		
		if (i5_val == 0 && !is_i5_const_one) {
			// I5=0时，两个LUT的输出必须相同
			if (lut1_output != lut2_output) {
				compatible = false;
				break;
			}
			// Z输出 = lut1的输出（也等于lut2的输出）
			lut6d_init_bits[addr] = lut1_output;
		} else {
			// I5=1时，输出可以不同
			// Z5输出 = lut2的输出，Z输出 = lut1的输出
			// 但这里我们只计算一个INIT值，所以选择其中一个
			if (addr < 32) {
				lut6d_init_bits[addr] = lut1_output;
			} else {
				lut6d_init_bits[addr] = lut2_output;
			}
		}
	}

	if (!compatible) {
		return {Const(0, 64), vector<SigBit>()};
	}

	return {RTLIL::Const(lut6d_init_bits), merged_inputs};
}

void updateDataStructuresAfterMerge(Module *module, Cell *old_lut1, Cell *old_lut2, Cell *new_lut6d) {
    log("更新合并后的数据结构\n");
    
    // 获取原LUT的层级信息
    size_t lut1_level = lut2levels.count(old_lut1) ? lut2levels[old_lut1] : 0;
    size_t lut2_level = lut2levels.count(old_lut2) ? lut2levels[old_lut2] : 0;
    
    size_t new_level = std::max(lut1_level, lut2_level);
    lut2levels[new_lut6d] = new_level;
    level2luts[new_level].insert(new_lut6d);
    // 从数据结构中移除原LUT
    if (lut2levels.count(old_lut1)) {
        size_t old_level1 = lut2levels[old_lut1];
        level2luts[old_level1].erase(old_lut1);
        lut2levels.erase(old_lut1);
    }
    
    if (lut2levels.count(old_lut2)) {
        size_t old_level2 = lut2levels[old_lut2];
        level2luts[old_level2].erase(old_lut2);
        lut2levels.erase(old_lut2);
    }
    // 更新bit2driver映射
    vector<SigBit> old1_outputs, old2_outputs, new_outputs;
    GetCellOutputsVector(old_lut1, old1_outputs);
    GetCellOutputsVector(old_lut2, old2_outputs);
    GetCellOutputsVector(new_lut6d, new_outputs);
    
    // 更新驱动映射
    for (auto bit : old1_outputs) {
        if (bit2driver.count(bit)) {
            bit2driver.erase(bit);
        }
    }
    for (auto bit : old2_outputs) {
        if (bit2driver.count(bit)) {
            bit2driver.erase(bit);
        }
    }
    for (auto bit : new_outputs) {
        bit2driver[bit] = new_lut6d;
    }
	// 清理cell2bits中的引用
	if (cell2bits.count(old_lut1)) {
		cell2bits.erase(old_lut1);
	}
	if (cell2bits.count(old_lut2)) {
		cell2bits.erase(old_lut2);
	}
	
	// 从module中删除旧的LUT cells
	module->remove(old_lut1);
	module->remove(old_lut2);
}

Cell *performLUT6DMerge(Module *module, Cell *lut1, Cell *lut2, const pool<SigBit> &common_inputs) {
    log("执行LUT6D合并: %s + %s\n", lut1->name.c_str(), lut2->name.c_str());
    
	vector<SigBit> lut1_inputs, lut2_inputs;
    GetCellInputsVector(lut1, lut1_inputs);
    GetCellInputsVector(lut2, lut2_inputs);
	SigBit lut1_output, lut2_output;
	lut1_output = GetCellOutput(lut1);
	lut2_output = GetCellOutput(lut2);

	// 计算合并后的输入集合
	pool<SigBit> lut1_unique, lut2_unique;
    for (auto input : lut1_inputs) {
        if (!common_inputs.count(input)) {
            lut1_unique.insert(input);
        }
    }
    for (auto input : lut2_inputs) {
        if (!common_inputs.count(input)) {
            lut2_unique.insert(input);
        }
    }
	
	vector<SigBit> lut6dInput;
	vector<SigBit> lut6dOutput;

	Const init = Const(0, 64);

	if (common_inputs.size() + lut1_unique.size() + lut2_unique.size() < 6) {
		lut6dOutput.push_back(lut2_output); // Z输出 (lut2作为主输出)
		lut6dOutput.push_back(lut1_output); // Z5输出 (lut1作为第二输出)
		auto result = checkCombinable(lut1, lut2, SigBit(State::S1));
		init = result.first;
		lut6dInput = result.second;
	}
	if (lut2_inputs.size() < 6 && init.bits() == Const(0, 64).bits()) {
		for (auto I5 : lut1_unique) {
			auto result = checkCombinable(lut2, lut1, I5);
			if (init != Const(0, 64)) {
				lut6dOutput.push_back(lut1_output);
				lut6dOutput.push_back(lut2_output);
				init = result.first;
				lut6dInput = result.second;
				break;
			}
		}
	}

	if (lut1_inputs.size() < 6 && init.bits() == Const(0, 64).bits()) {
		for (auto I5 : lut2_unique) {
			auto result = checkCombinable(lut1, lut2, I5); // Disabled for debugging
			
			if (init != Const(0, 64)) {
				lut6dOutput.push_back(lut2_output); // Z输出 (lut2作为主输出)
				lut6dOutput.push_back(lut1_output); // Z5输出 (lut1作为第二输出)
				init = result.first;
				lut6dInput = result.second;
				break;
			}
		}
	}
	if (init == Const(0, 64)) {
		log("无法合并LUT %s 和 %s\n", lut1->name.c_str(), lut2->name.c_str());
		return nullptr; // 返回空指针表示无法合并
	}
	
	log("成功合并LUT %s 和 %s，共有 %lu个输入\n", lut1->name.c_str(), lut2->name.c_str(), lut6dInput.size());
	string name = string(lut1->name.c_str()) + "_plus_" + string(lut2->name.c_str()) + "_LUT6D";
	Cell* newlut = addLut6d(module, name, lut6dInput, lut6dOutput, init);
	updateDataStructuresAfterMerge(module, lut1, lut2, newlut);

	return newlut; // 临时返回nullptr，实际应该返回新创建的LUT6D
}

void LUT6DRemapMain(Module *module, Cell *root) {
	// 获取当前节点的所有输入
	vector<SigBit> root_inputs;
	if (IsGTP_LUT6D(root)) {
		log("Root cell %s is a GTP_LUT6D, skipping remap.\n", root->name.c_str());
		return;
	}
	GetCellInputsVector(root, root_inputs);

	// 定义搜索的层数
	const size_t SEARCH_DEPTH = 2;

	log("Root LUT %s inputs: ", root->name.c_str());
	for (auto input : root_inputs) {
		log("%s ", log_signal(input));
	}
	log("\n");

	// 获取当前节点的输出
	vector<SigBit> root_outputs;
	GetCellOutputsVector(root, root_outputs);

	log("Root LUT %s outputs: ", root->name.c_str());
	for (auto output : root_outputs) {
		log("%s ", log_signal(output));
	}
	log("\n");

	// 获取root在拓扑结构中的层级
	size_t root_level = lut2levels[root];
	log("Root LUT %s is at level %lu\n", root->name.c_str(), root_level);

	dict<Cell *, pair<pool<SigBit>, pool<SigBit>>> compatible_cells; // Cell -> (common_inputs, all_inputs)

	// 收集往下n层拓扑结构中的所有节点
	std::queue<Cell *> worklist;
	pool<SigBit> root_input_set(root_inputs.begin(), root_inputs.end());
	worklist.push(root);
	
	if (SEARCH_DEPTH) {
		while(!worklist.empty() && compatible_cells.size() < 1000) {// 限制最大节点数，防止爆炸
			Cell *current = worklist.front();
			worklist.pop();
			if (lut2levels[current] > root_level + SEARCH_DEPTH) {
				continue; // 超过搜索深度，停止扩展
			}

			vector<SigBit> outputs;
			GetCellOutputsVector(current, outputs);

			for (auto &outbit : outputs) {
				if (bit2reader.count(outbit)) {
					for (auto reader : bit2reader[outbit]) {
						if (IsCombinationalCell(reader) && 
							lut2levels.count(reader) && 
							lut2levels[reader] > root_level) {
							if (lut2levels[reader] <= root_level + SEARCH_DEPTH) {
								vector<SigBit> cell_inputs;
								GetCellInputsVector(reader, cell_inputs);
								
								// 计算共同输入和总输入数量
								size_t common_count = 0;
								size_t total_count = root_input_set.size();
								
								for (const auto& input : cell_inputs) {
									if (root_input_set.count(input)) {
										common_count++;
									} else {
										total_count++;
									}
								}

								if (common_count > 0 && total_count <= 6) {
									// 只有在需要时才构建完整的pool
									pool<SigBit> cell_input_set(cell_inputs.begin(), cell_inputs.end());
									pool<SigBit> common_inputs;
									for (const auto& input : root_input_set) {
										if (cell_input_set.count(input)) {
											common_inputs.insert(input);
										}
									}
									
									pool<SigBit> all_inputs = root_input_set;
									all_inputs.insert(cell_input_set.begin(), cell_input_set.end());
									
									compatible_cells[reader].first = std::move(common_inputs);
									compatible_cells[reader].second = std::move(all_inputs);
									log("Found compatible cell %s at level %lu with %lu common inputs, total inputs would be %lu\n",
										reader->name.c_str(), lut2levels[reader], common_count, total_count);
								}
							}
							if (lut2levels[reader] == root_level + SEARCH_DEPTH) {
								worklist.push(reader);
							}
						}
					}
				}
			}
		}
	}

	// 在同一层中找到与root有相同输入的节点
	for (auto cell : level2luts[root_level]) {
		if (cell == root) continue;
		
		vector<SigBit> cell_inputs;
		GetCellInputsVector(cell, cell_inputs);
		
		// 先快速计算数量，避免不必要的pool构建
		size_t common_count = 0;
		size_t total_count = root_input_set.size();
		
		for (const auto& input : cell_inputs) {
			if (root_input_set.count(input)) {
				common_count++;
			} else {
				total_count++;
			}
		}
		
		// 只有在满足条件时才构建完整的pool
		if (total_count < 6 || (common_count > 0 && total_count <= 6 && !(root_inputs.size() == 6 && cell_inputs.size() == 6))) {
			pool<SigBit> cell_input_set(cell_inputs.begin(), cell_inputs.end());
			pool<SigBit> common_inputs;
			for (const auto& input : root_input_set) {
				if (cell_input_set.count(input)) {
					common_inputs.insert(input);
				}
			}
			
			pool<SigBit> all_inputs = root_input_set;
			all_inputs.insert(cell_input_set.begin(), cell_input_set.end());
			
			compatible_cells[cell].first = std::move(common_inputs);
			compatible_cells[cell].second = std::move(all_inputs);
			log("Found compatible same-level cell %s with %lu common inputs, total inputs would be %lu\n",
				cell->name.c_str(), common_count, total_count);
		}
	}

	log("Total compatible cells found: %lu\n", compatible_cells.size());

	if (compatible_cells.empty()) {
		return;
	}
		
	// 将compatible_cells转换为vector以便排序
	vector<pair<Cell *, pair<pool<SigBit>, pool<SigBit>>>> sorted_compatible_cells;
	sorted_compatible_cells.reserve(compatible_cells.size());
	for (auto &entry : compatible_cells) {
		sorted_compatible_cells.emplace_back(entry.first, std::move(entry.second));
	}
	
	// 按照共同输入数量由大到小、总输入数量由大到小排序
	std::sort(sorted_compatible_cells.begin(), sorted_compatible_cells.end(),
		[](const auto &a, const auto &b) {
			size_t common_a = a.second.first.size();    // 共同输入数量
			size_t common_b = b.second.first.size();
			size_t total_a = a.second.second.size();    // 总输入数量
			size_t total_b = b.second.second.size();
			
			// 首先按共同输入数量降序
			if (common_a != common_b) {
				return common_a > common_b;
			}
			// 共同输入数量相同时，按总输入数量降序
			return total_a > total_b;
		});
	
	while(!sorted_compatible_cells.empty()) {
		Cell *candidate = sorted_compatible_cells[0].first;
		
		// 检查candidate是否仍然存在（可能已被之前的合并删除）
		if (!module->cells_.count(candidate->name)) {
			log("Candidate LUT %s has been removed, skipping\n", candidate->name.c_str());
			sorted_compatible_cells.erase(sorted_compatible_cells.begin());
			continue;
		}
		
		// 检查root是否仍然存在（虽然不太可能，但为了安全）
		if (!module->cells_.count(root->name)) {
			log("Root LUT %s has been removed, stopping merge attempts\n", root->name.c_str());
			return;
		}
		
		pool<SigBit> common_inputs = sorted_compatible_cells[0].second.first;
		pool<SigBit> all_inputs = sorted_compatible_cells[0].second.second;
		if (performLUT6DMerge(module, root, candidate, common_inputs) != nullptr) {
			return;
		}
		sorted_compatible_cells.erase(sorted_compatible_cells.begin());
	}
}

// Internal implementation (static, not exported)
static void LUT6DRemap_impl(Module *module){
	log("Perform remapping after mapping\n");

	MapperInit(module);
	CheckCellWidth(module);
	vector<Cell *> initial_luts;
	GetTopoSortedLUTs(module, initial_luts);
	
	// 按层级处理，从低层级到高层级
	for (size_t level = 0; level < level2luts.size(); level++) {
		if (!level2luts.count(level)) continue;
		
		log("Processing level %lu with %lu LUTs\n", level, level2luts[level].size());
		
		// 收集当前层级的LUT名称（避免指针问题）
		vector<string> level_lut_names;
		for (auto lut : level2luts[level]) {
			if (lut && module->cells_.count(lut->name)) {
				level_lut_names.push_back(lut->name.c_str());
			}
		}
		
		// 处理当前层级的每个LUT
		for (const auto& lut_name : level_lut_names) {
			// 每次都重新从module获取Cell指针，确保有效性
			if (!module->cells_.count(IdString(lut_name))) {
				// log("LUT %s has been removed, skipping\n", lut_name.c_str());
				continue;
			}
			
			Cell *lut = module->cell(IdString(lut_name));
			if (!lut) {
				// log("LUT %s pointer is null, skipping\n", lut_name.c_str());
				continue;
			}
			
			if (IsGTP_LUT6D(lut)) {
				// log("LUT %s is a GTP_LUT6D, skipping\n", lut->name.c_str());
				continue;
			}
			
			if (!lut2levels.count(lut)) {
				// log("LUT %s not found in lut2levels, skipping\n", lut->name.c_str());
				continue;
			}
			
			// log("Processing LUT %s at level %lu\n", lut->name.c_str(), lut2levels[lut]);
			LUT6DRemapMain(module, lut);
			// log("Finished processing LUT %s\n", lut->name.c_str());
		}
	}
	


	return;
}

// Public wrapper function declared in lut6d_remap.h
bool LUT6DRemap(Module *module) {
	LUT6DRemap_impl(module);
	return true;
}

YOSYS_NAMESPACE_END
