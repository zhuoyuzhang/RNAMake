//
//  thermo_simulation.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "thermo_simulation.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_topology.h"
#include "motif_data_structures/motif_graph.h"

void
ThermoSimulationApp::setup_options() {
    add_option("mg", String(""), OptionType::STRING, true);

    /*add_option("f", String(""), OptionType::STRING, true);
    add_option("break_point", String(""), OptionType::STRING, true);
    add_option("connections", String(""), OptionType::STRING, false);
    add_option("o", String("sequences.dat"), OptionType::STRING, false);
    add_cl_options(optimizer_.options(), "optimizer");*/
}

void
ThermoSimulationApp::parse_command_line(
    int argc,
    const char ** argv) {
    
    Application::parse_command_line(argc, argv);
    //cl_parser_.assign_options(cl_options_, optimizer_.options(), "optimizer");
}

void
ThermoSimulationApp::run() {
    
    
    
}

int main(int argc, const char * argv[]) {
    
    //load tectos
    auto tecto_dir = String(base_dir()+"/rnamake/lib/RNAMake/apps/simulate_tectos");
    RM::instance().add_motif(tecto_dir+"/resources/GAAA_tetraloop");
    RM::instance().add_motif(tecto_dir+"/resources/GGAA_tetraloop");
    
    auto app = ThermoSimulationApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();
    
    return 0;

}