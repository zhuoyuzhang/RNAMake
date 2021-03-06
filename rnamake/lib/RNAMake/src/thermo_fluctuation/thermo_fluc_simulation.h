//
//  thermo_fluc_simulation.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_simulation__
#define __RNAMake__thermo_fluc_simulation__

#include <stdio.h>

#include "base/types.h"
#include "base/option.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_scorer.h"


class ThermoFlucSimulationException : public std::runtime_error {
public:
    ThermoFlucSimulationException(
        String const & message) :
    std::runtime_error(message)
    {}
};

class ThermoFlucSimulation : public OptionClass {
public: // constructors
    ThermoFlucSimulation();

    ~ThermoFlucSimulation() {}
    
public:
    
    void
    setup(
        MotifStateEnsembleTreeOP const &,
        int,
        int,
        int,
        int);
    
public: //run
    
    inline
    int
    _check_sterics() {
        clash_ = 0;
        for(auto const & i : check_nodes_1_) {
            for(auto const & j : check_nodes_2_) {
                for(auto const & b2 : sampler_.mst()->get_node(i)->data()->cur_state->beads()) {
                    for(auto const & b1 : sampler_.mst()->get_node(j)->data()->cur_state->beads()) {
                        if(b1.distance(b2) < 2.2) { clash_ = 1; }
                    }
                    
                    if(clash_) { break; }
                }
                if(clash_) { break; }
            }
            if(clash_) { break; }
        }
        return clash_;
    }
    
    int
    run();
    
public: //option wrappers
    
    inline
    Options &
    options() { return options_; }
    
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
    void
    update_var_options();
    
protected:
    void
    setup_options();
    
private:
    ThermoFlucScorerOP scorer_;
    ThermoFlucSampler sampler_;
    BasepairStateOP end_state_1_, end_state_2_;
    Options options_;
    int ni1_, ni2_, ei1_, ei2_;
    int clash_;
    float score_;
    Ints check_nodes_1_, check_nodes_2_;
    int setup_;
    //option vars
    bool record_;
    float temperature_, cutoff_;
    int steps_;
};


#endif /* defined(__RNAMake__thermo_fluc_simulation__) */
