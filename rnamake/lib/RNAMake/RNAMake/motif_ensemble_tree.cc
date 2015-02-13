//
//  motif_ensemble_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <queue>
#include "motif_ensemble_tree.h"
#include "motif_tree_state.h"
#include "motif_tree_state_tree.h"

MotifEnsembleTree::MotifEnsembleTree() {
    MotifEnsemble me;
    MotifTreeStateOP mts (new MotifTreeState(ref_mts()));
    MotifState ms ( mts, 1.0);
    me.add_motif_state(ms);
    MotifEnsembleTreeNodeOP head ( new MotifEnsembleTreeNode(me, NULL, 0));
    nodes_ = MotifEnsembleTreeNodeOPs();
    nodes_.push_back(head);
    last_node_ = head;
}

MotifEnsembleTree::MotifEnsembleTree(MotifEnsemble const & me) {
    MotifEnsembleTreeNodeOP head ( new MotifEnsembleTreeNode(me, NULL, 0));
    nodes_ = MotifEnsembleTreeNodeOPs();
    nodes_.push_back(head);
    last_node_ = head;
}

MotifEnsembleTreeNodeOP
MotifEnsembleTree::add_ensemble(
    MotifEnsemble const & ensemble,
    MotifEnsembleTreeNodeOP parent,
    int parent_end_index) {
    
    if (parent == NULL) { parent = last_node_; }
    if (parent_end_index == -1) { parent_end_index = parent->available_ends()[0]; }
    
    MotifEnsembleTreeNodeOP new_node ( new MotifEnsembleTreeNode(ensemble, parent, (int)nodes_.size()));
    parent->add_child(new_node, parent_end_index);
    nodes_.push_back(new_node);
    last_node_ = new_node;
    return new_node;
}

MotifTreeStateTree
MotifEnsembleTree::get_mtst() {
    MotifTreeStateTree mtst;
    String mts_name = nodes_[0]->motif_ensemble().motif_states()[0].mts->name();
    if( mts_name.compare("start") != 0) {
        MotifTreeStateOP mts = nodes_[0]->motif_ensemble().motif_states()[0].mts;
        mtst = MotifTreeStateTree(mts);
    }
    mtst.sterics(0);
    std::queue<MotifEnsembleTreeNodeOP> open_nodes;
    open_nodes.push(nodes_[0]);
    int i = 0;
    while (! open_nodes.empty()) {
        MotifEnsembleTreeNodeOP current = open_nodes.front();
        open_nodes.pop();
        i = -1;
        for( auto const & n : current->children() ) {
            i++;
            if (n == NULL) { continue; }
            MotifTreeStateOP mts = n->motif_ensemble().motif_states()[0].mts;
            MotifTreeStateNodeOP parent = mtst.nodes() [ current->index() ];
            MotifTreeStateNodeOP mstn = mtst.add_state(mts, parent, i);
            open_nodes.push(n);
        }
    }
    
    return mtst;
}


Ints
MotifEnsembleTreeNode::available_ends() {
    Ints indices;
    int i = -1;
    for( auto const & c : children_) {
        i++;
        if(i == start_index_) { continue; }
        if(c == NULL) { indices.push_back(i); }
    }
    return indices;
}


void
MotifEnsembleTreeNode::add_child(
    MotifEnsembleTreeNodeOP const & node,
    int end_direction) {
    children_[end_direction] = node;
}





