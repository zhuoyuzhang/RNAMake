#ifndef __tools_build_motif_tree
#define __tools_build_motif_tree

#include <map>

#include "builder_graph.hpp"

#include "motif_data_structures/motif_tree.h"
#include "resources/motif_sqlite_library.h"


class MotifTreeBuilder {
public:
    MotifTreeBuilder() {
        g_ = helix_and_two_way();
        mlibs_ = std::map<String, MotifSqliteLibraryOP>();
    }
    
    MotifTreeBuilder(
        BuilderGraphOP const & g) {
        
        g_ = g;
        mlibs_ = std::map<String, MotifSqliteLibraryOP>();
    }
    
    
public:
    MotifTreeOP
    build(int repeat=1) {
        auto mt = std::make_shared<MotifTree>();
        
        for(int i = 0; i < repeat; i++) {
            auto index_map = std::map<int, int>();
            int j = 0, pos = 0, parent_index = -1, parent_end_index = -1;
            for(auto const & n : *g_) {
                auto lib_name = _get_lib_name(n->data());
                auto mlib = _get_lib(lib_name);
                
                if(j != 0) {
                    parent_index = n->parent()->index();
                    parent_end_index = n->parent_end_index();
                    
                    parent_index = index_map[parent_index];
                }
                
                pos = _add_motif_to_tree(mt, mlib, parent_index, parent_end_index);
                
                if(pos == -1) { return mt; }
                index_map[n->index()] = pos;
                j++;
                
            }
        }
        
        return mt;
    }
    
private:
    String
    _get_lib_name(MotifType const & type) {
        if     (type == HELIX)   { return "ideal_helices"; }
        else if(type == TWOWAY)  { return "twoway"; }
        else if(type == NWAY)    { return "nway"; }
        else if(type == HAIRPIN) { return "hairpin"; }
        else if(type == TCONTACT){ return "tcontact"; }
        else { throw std::runtime_error("cannot convert type"); }
        
    }
    
    MotifSqliteLibraryOP const &
    _get_lib(String const & name) {
        if(mlibs_.find(name) == mlibs_.end()) {
            mlibs_[name] = std::make_shared<MotifSqliteLibrary>(name);
        }
        
        return mlibs_[name];
    }
    
    inline
    int
    _add_motif_to_tree(
        MotifTreeOP & mt,
        MotifSqliteLibraryOP const & mlib,
        int parent_index,
        int parent_end_index) {
        
        int pos = -1, count = 0;
        while(pos == -1) {
            auto m = mlib->get_random();
            
            pos = mt->add_motif(m, parent_index, parent_end_index);
         
            count++;
            if(count > 100) {
                return 0;
            }
        }

        return pos;
    }
    
    
private:
    BuilderGraphOP g_;
    std::map<String, MotifSqliteLibraryOP> mlibs_;
    
};


#endif





















