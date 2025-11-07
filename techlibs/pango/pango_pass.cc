#include "kernel/yosys.h"
#include "kernel/sigtools.h"
#include "kernel/celltypes.h"
#include "kernel/modtools.h"
#include "kernel/satgen.h"
// Use headers instead of including implementation .cc files
// #include "techlibs/pango/lut6d_remap.h"
#include "techlibs/pango/synth_pango.h"
#include "techlibs/pango/lut6d_map.h"
#include "techlibs/pango/lut6d_remap.h"
#include <queue>
#include <algorithm>

USING_YOSYS_NAMESPACE
using namespace std;
YOSYS_NAMESPACE_BEGIN

struct MapperPass : public Pass {
	MapperPass() : Pass("mapper", "synthesis for Pango FPGA. Mapper main function.") {}
	void help() override
	{
		//   |---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|
		log("\n");
		log("    mapper [options] [selection]\n");
	}
	bool write_out_black_list;
	void clear_flags() override { write_out_black_list = false; }
	void execute(std::vector<std::string> args, RTLIL::Design *design) override
	{
		log_header(design, "Start MapperPass\n");

		size_t argidx;
		for (argidx = 1; argidx < args.size(); argidx++) {
			if (args[argidx] == "-ilut") {
				using_internel_lut_type = true;
				continue;
			}
			if (args[argidx] == "-interation" && argidx + 1 < args.size()) {
				MAX_INTERATIONS = max(atoi(args[++argidx].c_str()), 3);
				continue;
			}
			break;
		}
		extra_args(args, argidx, design);

		Module *module = design->top_module();
		if (module == nullptr)
			log_cmd_error("No top module found.\n");

		log_header(design, "Continuing MapperPass pass.\n");
		MapperInit(module);
		MapperMain(module);
		log_pop();
	}
} MapperPass;

struct SynthPangoPass : public ScriptPass {
	SynthPangoPass() : ScriptPass("synth_pango", "synthesis script pass for Pango FPGA. Map gate to GTP_LUT.") {}
	void help() override
	{
		//   |---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|
		log("\n");
		log("    synth_pango [options] [selection]\n");
	}

	string input_verilog_file;
	string output_verilog_file;
	string top_module_name;
	void clear_flags() override
	{
		using_internel_lut_type = false;
		output_verilog_file = "";
		top_module_name = "";
	}
	void execute(std::vector<std::string> args, RTLIL::Design *design) override
	{
		string run_from, run_to;
		clear_flags();

		size_t argidx;
		for (argidx = 1; argidx < args.size(); argidx++) {
			if (args[argidx] == "-top" && argidx + 1 < args.size()) {
				top_module_name = args[++argidx];
				continue;
			}
			if (args[argidx] == "-input" && argidx + 1 < args.size()) {
				input_verilog_file = args[++argidx];
				continue;
			}
			if (args[argidx] == "-interation" && argidx + 1 < args.size()) {
				MAX_INTERATIONS = max(atoi(args[++argidx].c_str()), 3);
				continue;
			}
			if (args[argidx] == "-run" && argidx + 1 < args.size()) {
				size_t pos = args[argidx + 1].find(':');
				if (pos == std::string::npos)
					break;
				run_from = args[++argidx].substr(0, pos);
				run_to = args[argidx].substr(pos + 1);
				continue;
			}
			break;
		}

		extra_args(args, argidx, design);
		if (!design->full_selection())
			log_cmd_error("This command only operates on fully selected designs!\n");

		log_header(design, "Start synth_pango\n");
		log_push();

		run_script(design, run_from, run_to);

		log_pop();
	}
	void script() override
	{
		if (check_label("begin")) {
#if defined(_WIN32)
			run("read_verilog -lib ./techlibs/pango/pango_lib.v");
#else
			run("read_verilog -lib +/pango/pango_lib.v");
			// run("read_verilog -lib ./pango_lib.v");
#endif
			run(stringf("read_verilog -icells %s", input_verilog_file.c_str()));
			if (top_module_name.size() > 0) {
				run(stringf("hierarchy -check -top %s", top_module_name.c_str()));
			}
		}

		Module *module = active_design->top_module();
		if (module == nullptr)
			log_cmd_error("No top module found.\n");
		if (top_module_name.size() == 0) {
			top_module_name = module->name.c_str();
		}

		if (check_label("pango")) {
			run(stringf("hierarchy -check -top %s;;", top_module_name.c_str()));
			MapperInit(module);
			LUT6DMapper(module);
			MapperInit(module);
			MapperMain(module);
			// run("lut6d_map");
			// run("lut6d_remap");
		}
		if (check_label("check")) {
			run("check -mapped");
		}
		if (check_label("verilog")) {
			run(stringf("write_verilog -noexpr -noattr +/../techlibs/pango/outputs/%s_syn.v", top_module_name.c_str() + 1));
		}
		if (check_label("score")) {
			run(stringf("score -before %s -after %s_syn.v -out +/../techlibs/pango/outputs/%s_score.txt", input_verilog_file.c_str(), top_module_name.c_str() + 1, top_module_name.c_str() + 1));
		}
	}
} SynthPangoPass;

struct LUT6DMapperPass : public Pass {
    LUT6DMapperPass() : Pass("lut6d_map", "Map logic to dual-output LUT6D cells") { }
    
    void help() override {
        log("\n");
        log("    lut6d_map [options] [selection]\n");
        log("\n");  
        log("This pass maps single-output logic gates to dual-output LUT6D cells\n");
        log("to reduce area by merging compatible LUTs.\n");
        log("\n");
        log("The algorithm:\n");
        log("  1. Generates cuts for each combinational gate\n");
        log("  2. Finds candidate gates that can be merged into LUT6D\n");
        log("  3. Verifies functional compatibility using SAT\n");
        log("  4. Creates LUT6D cells and removes merged gates\n");
        log("\n");
        log("Options:\n");
        log("  (none currently supported)\n");
        log("\n");
    }
    
    void execute(std::vector<std::string> args, RTLIL::Design *design) override {
        log_header(design, "Executing LUT6D Mapper\n");
        
        size_t argidx;
        for (argidx = 1; argidx < args.size(); argidx++) {
            // Future: add options here
            break;
        }
        extra_args(args, argidx, design);
        
        for (auto module : design->selected_modules()) {
            if (module->has_processes_warn())
                continue;
                
            log("Processing module %s for LUT6D mapping\n", log_id(module));
            
            // Initialize mapper infrastructure
            MapperInit(module);
            
            // Run LUT6D mapping algorithm
            LUT6DMapper(module);
            
            log("Module %s LUT6D mapping complete\n", log_id(module));
        }
        
        log_header(design, "LUT6D Mapper completed\n");
    }
} LUT6DMapperPass;

struct LUT6DRemapPass : public Pass {
    LUT6DRemapPass() : Pass("lut6d_remap", "简化的双输出LUT重组优化") { }
    
    void help() override {
        log("\n");
        log("    lut6d_remap [selection]\n");
        log("\n");
        log("将兼容的单输出LUT合并为双输出LUT6D以减少面积\n");
        log("\n");
    }
    
    void execute(std::vector<std::string> args, RTLIL::Design *design) override {
        log_header(design, "执行LUT6D重组优化\n");
        
        size_t argidx;
        for (argidx = 1; argidx < args.size(); argidx++) {
            break;
        }
        extra_args(args, argidx, design);
        
        for (auto module : design->selected_modules()) {
            if (module->has_processes_warn())
                continue;
                
            // log("处理模块 %s\n", log_id(module));
            
            // 简化的重组逻辑
            LUT6DRemap(module);
            
            // log("模块 %s 重组完成\n", log_id(module));
        }
    }
} LUT6DRemapPass;

YOSYS_NAMESPACE_END