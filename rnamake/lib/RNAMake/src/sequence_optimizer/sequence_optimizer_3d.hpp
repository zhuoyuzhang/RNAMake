//
//  sequence_optimizer_3d.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/30/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef sequence_optimizer_3d_hpp
#define sequence_optimizer_3d_hpp

#include <stdio.h>

#include "base/types.h"
#include "base/option.h"
#include "util/random_number_generator.h"
#include "eternabot/scorer.h"
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_data_structures/motif_state_tree.h"
#include "motif_data_structures/motif_state_graph.hpp"


class SequenceOptimizerScorer {
public:
    virtual
    float
    score(MotifStateGraphOP const &) = 0;
};

typedef std::shared_ptr<SequenceOptimizerScorer> SequenceOptimizerScorerOP;


class ExternalTargetScorer : public SequenceOptimizerScorer {
public:
    ExternalTargetScorer(
        BasepairStateOP const & target,
        int ni,
        int ei):
        target_(target),
        ni_(ni),
        ei_(ei) {}
    
public:
    float
    score(MotifStateGraphOP const & msg) {
        state_ = msg->get_node(ni_)->data()->get_end_state(ei_);
        return target_->diff(state_);
    }
    
private:
    BasepairStateOP target_, state_;
    int ni_, ei_;
    
};


class InternalTargetScorer : public SequenceOptimizerScorer {
public:
    InternalTargetScorer(
        int ni1,
        int ei1,
        int ni2,
        int ei2):
    ni1_(ni1),
    ei1_(ei1),
    ni2_(ni2),
    ei2_(ei2) {}
    
public:
    float
    score(MotifStateGraphOP const & msg) {
        state1_ = msg->get_node(ni1_)->data()->get_end_state(ei1_);
        state2_ = msg->get_node(ni2_)->data()->get_end_state(ei2_);
        return state1_->diff(state2_);
    }
    
private:
    int ni1_, ni2_;
    int ei1_, ei2_;
    BasepairStateOP state1_, state2_;
};


class SequenceOptimizer3D  {
public:
    
    SequenceOptimizer3D();
    
    ~SequenceOptimizer3D() {}
    
public: //setup
    
    void
    set_scorer(SequenceOptimizerScorerOP const & scorer) {
        scorer_ = scorer;
    }
    
private:
    
    struct OptimizedSequence {
        String sequence;
        float dist_score, eterna_score;
    };
    
    typedef std::shared_ptr<OptimizedSequence> OptimizedSequenceOP;
    typedef std::vector<OptimizedSequenceOP> OptimizedSequenceOPs;
    
    struct DesignableBP {
        inline
        DesignableBP(
            sstruct::BasepairOP const & nbp):
            bp(nbp),
            last_state(Strings{"", ""}),
            m_id_bot(nullptr),
            m_id_top(nullptr)
            {}
        
        void
        update_state(
            Strings const & bp_name) {
            last_state[0] = bp->res1()->name();
            last_state[1] = bp->res2()->name();
            bp->res1()->name(bp_name[0]);
            bp->res2()->name(bp_name[1]);
        }
        
        void
        revert_state() {
            bp->res1()->name(last_state[0]);
            bp->res2()->name(last_state[1]);
        }
        
        
        sstruct::BasepairOP bp;
        Strings last_state;
        std::shared_ptr<Uuid> m_id_bot, m_id_top;
        
    };
    
    typedef std::shared_ptr<DesignableBP> DesignableBPOP;
    typedef std::vector<DesignableBPOP> DesignableBPOPs;
    
public:
    
    inline
    OptimizedSequenceOPs
    get_optimized_sequences(
        MotifGraphOP const & mg,
        SequenceOptimizerScorerOP const & scorer) {
        set_scorer(scorer);
        return get_optimized_sequences(mg);
    }
    
    OptimizedSequenceOPs
    get_optimized_sequences(
        MotifGraphOP const &);
    
    inline
    MotifGraphOP
    get_optimized_mg(
        MotifGraphOP const & mg,
        SequenceOptimizerScorerOP const & scorer) {
        set_scorer(scorer);
        return get_optimized_mg(mg);
    }
    
    MotifGraphOP
    get_optimized_mg(
        MotifGraphOP const &);
    
private:
    void
    _update_designable_bp(
        DesignableBPOP const &,
        MotifStateGraphOP &,
        sstruct::PoseOP &);
    
    String
    _validate_sequence(
        MotifStateGraphOP const &,
        sstruct::PoseOP const &);
    
    DesignableBPOPs
    _get_designable_bps(
        sstruct::PoseOP &);
    
    void
    _initiate_sequence_in_msg(
        MotifStateGraphOP &,
        sstruct::PoseOP const &);
    
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
    Options options_;
    eternabot::Scorer eterna_scorer_;
    RandomNumberGenerator rng_;
    SequenceOptimizerScorerOP scorer_;
    std::vector<Strings> possible_bps_;
    // option vars
    int solutions_;
    int steps_;
    float cutoff_, eterna_cutoff_;
    bool verbose_, return_lowest_;
    
    
};

#endif /* sequence_optimizer_3d_hpp */
