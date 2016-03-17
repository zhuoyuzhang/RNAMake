import unittest
import time
from rnamake import motif_state, motif, util
from rnamake import resource_manager as rm


class MotifStateUnittest(unittest.TestCase):
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print "%s: %.3f" % (self.id(), t)

    def test_creation(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms = motif_state.get_motif_state(m)

    def test_speed_test(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms = motif_state.get_motif_state(m)

        for i in range(100):
            ms_copy = ms.copy()

    def test_speed_test_2(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms = m.get_state()

        for i in range(100):
            ms_copy = ms.copy()

    def test_align(self):
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        m2 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        motif.align_motif(m1.ends[1].state(), m2.ends[0], m2)

        ms1 = motif_state.get_motif_state(m1)
        ms2 = motif_state.get_motif_state(m2)
        motif_state.align_motif_state(ms1.ends[1], ms2)

        m2_new = ms2.to_motif()
        atoms1 = m2.structure.atoms()
        atoms2 = m2_new.structure.atoms()
        for i, a in enumerate(atoms1):
            dist = util.distance(a.coords, atoms2[i].coords)
            if dist > 0.1:
                self.fail("did not get the right coordinate transformation")

    def test_align_speed(self):
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms1 = motif_state.get_motif_state(m1)
        ms2 = motif_state.get_motif_state(m1)

        for i in range(100):
            motif_state.align_motif_state(ms1.ends[1], ms2)

    def test_align_speed_2(self):
        m1 = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms1 = m1.get_state()
        ms2 = m1.get_state()

        for i in range(100):
            motif.align_motif_state(ms1.end_states[1], ms2)

    def test_to_str_residue(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms = motif_state.get_motif_state(m)

        r = ms.residues()[1]
        r_str = r.to_str()
        r_new = motif_state.str_to_residue(r_str)

        for i in range(len(r.beads)):
            dist = util.distance(r.beads[i].center, r_new.beads[i].center)
            if dist > 0.1:
                self.fail()

    def test_to_str_chain(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms = motif_state.get_motif_state(m)

        c = ms.structure.chains[0]
        c_str = c.to_str()
        c_new = motif_state.str_to_chain(c_str)

    def test_to_str(self):
        m = rm.manager.get_motif(name="HELIX.IDEAL.2")
        ms = motif_state.get_motif_state(m)
        s = ms.to_str()

        ms_copy = motif_state.str_to_motif(s)
        motif_state.align_motif_state(ms.ends[1], ms_copy)




if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(MotifStateUnittest)
    unittest.TextTestRunner(verbosity=0).run(suite)