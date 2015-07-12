import unittest
import build
import rnamake.motif_tree
import rnamake.resource_manager as rm
import rnamake.secondary_structure_factory as ssfactory
import util

class MotifTreeUnittest(unittest.TestCase):

    def test_creation(self):
        mt = rnamake.motif_tree.MotifTree()

    def test_add_motif(self):
        mt = rnamake.motif_tree.MotifTree()
        m = rm.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)
        if len(mt) != 2:
            self.fail("did not add motifs properly")

    def test_remove_node(self):
        mt = rnamake.motif_tree.MotifTree()
        m = rm.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)

        mt.remove_node()
        if len(mt) != 1:
            self.fail("did not remove node correctly")

    def test_merge(self):
        mt = rnamake.motif_tree.MotifTree()
        m = rm.manager.get_motif("HELIX.IDEAL.2")
        mt.add_motif(m)
        mt.add_motif(m)
        p = mt.to_pose(chain_closure=0)


        """return
        f = open("test.out")
        s = f.readline()
        f.close()

        mt = rnamake.motif_tree.str_to_motif_tree(s)
        pose = mt.get_pose()
        if len(pose.ends) != 2:
            self.fail("did not merge properly")"""

    def test_readd(self):
        mt = rnamake.motif_tree.MotifTree()
        rm = rnamake.resource_manager.ResourceManager()
        m = rm.get_motif("HELIX.IDEAL")
        mt.add_motif(m)
        m = mt.nodes[1].motif
        mt.remove_node(mt.last_node)
        node = mt.add_motif(m)
        print node
        mt.write_pdbs()

    def test_replace_ideal_helices(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        mt.write_pdbs()
        p = mt.to_pose()
        seq = p.optimized_sequence()
        db = p.dot_bracket()

        ss = ssfactory.factory.get_structure(seq, db)
        connectivity = ss.motif_topology_from_end(ss.ends[0])
        mt1 = rnamake.motif_tree.motif_tree_from_topology(connectivity)

    def test_secondary_structure(self):
        builder = build.BuildMotifTree()
        mt = builder.build()
        p = mt.to_pose()
        print p.secondary_structure.copy().ends





def main():
    unittest.main()

if __name__ == '__main__':
    main()
