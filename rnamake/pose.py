import motif

class Pose(motif.Motif):
    #TODO finish implementing
    def __init__(self, mdir=None, pdb=None):
        self.structure = None
        self.mdir, self.name, self.ends = "", "", []
        self.beads, self.score, self.basepairs = [], 0, []
        self.designable = {}
        self.motif_bounds = []
        self._setup(mdir, pdb)
