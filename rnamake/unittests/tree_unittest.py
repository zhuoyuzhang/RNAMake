import unittest
import rnamake.tree as tree

class TreeUnittest(unittest.TestCase):

    def test_creation(self):
        t = tree.TreeDynamic()
        t.add_data(0)
        t.add_data(1)

        if len(t) != 2:
            self.fail("did not add nodes properly")

        t = tree.TreeStatic()
        t.add_data(0, 2)
        t.add_data(1, 2, 0, 0)
        t.add_data(2, 2, 0, 1)

        t = tree.TreeDynamic()
        t.add_data(0)
        t.add_data(1)
        t.add_data(2, parent_index=0)
        t.add_data(3, parent_index=0)

        self.failIf(len(t.get_node(0).children) != 3,
                    "node does not have the correct number of children")

    def test_remove(self):
        t = tree.TreeDynamic()
        t.add_data(0)
        t.add_data(1)
        t.add_data(2)
        t.remove_node(index=2)
        if len(t) != 2:
            self.fail("did not properly remove node")

    def test_remove_node_level(self):
        t = tree.TreeDynamic()
        t.add_data(0)
        t.increase_level()
        t.add_data(1)
        t.add_data(2)
        t.remove_node_level()
        if len(t) != 1:
            self.fail("did not properly remove node level")



def main():
    unittest.main()

if __name__ == '__main__':
    main()
