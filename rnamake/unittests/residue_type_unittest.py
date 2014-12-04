import unittest
import rnamake.residue_type


class ResidueTypeUnittest(unittest.TestCase):
    def test_creation(self):
        try:
            rtypes = rnamake.residue_type.rtypes
        except:
            self.fail("test_creation test failed")

    def test_slots(self):
        rtypes = rnamake.residue_type.rtypes
        try:
            rtypes.v1 = 1
        except AttributeError:
            pass
        except:
            self.fail("ddi not expect this error in test_slots")

    def test_get_rtype(self):
        gua_rtype = rnamake.residue_type.get_rtype("GUA")
        print gua_rtype
        if gua_rtype is None:
            self.fail("cannot get GUA restype in test_get_rtype")

        none_rtype = rnamake.residue_type.get_rtype("GUA1")
        if none_rtype is not None:
            self.fail("did not return None properly in test_get_rtype")

        altname_rtype = rnamake.residue_type.get_rtype("I")
        if altname_rtype is None:
            self.fail("cannot get altname restype in test_get_rtype")


def main():
    unittest.main()

if __name__ == '__main__':
    main()
