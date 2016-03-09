import primitives
import basic_io
import util

import uuid
import numpy as np
import motif
import resource_manager as rm

def get_motif_state(m):
    #get structure state
    c_states = []
    for c in m.structure.chains:
        r_states = []
        for r in c.residues:
            r_state = Residue(r.name, r.num, r.chain_id, r.i_code,
                              r.uuid, r.get_beads())
            r_states.append(r_state)
        c_state = Chain(r_states)
        c_states.append(c_state)
    struc = Structure(c_states)

    ends = []
    for end in m.ends:
        res1 = struc.get_residue(uuid=end.res1.uuid)
        res2 = struc.get_residue(uuid=end.res2.uuid)
        bp_state = Basepair(end.r(), end.d(), end.state().sugars, res1, res2,
                            end.name())
        ends.append(bp_state)

    return Motif(struc, ends, end_ids=m.end_ids, name=m.name)


class Residue(primitives.Residue):

    __slots__ = [
        "name",
        "num",
        "chain_id",
        "i_code",
        "uuid",
        "beads"
    ]

    def __init__(self, name, num, chain_id, i_code="",
                 r_uuid=None, beads=None):

        super(self.__class__, self).__init__(name, num, chain_id, i_code)
        if r_uuid is None:
           self.uuid = uuid.uuid1()
        else:
            self.uuid = r_uuid
        if beads is None:
            self.beads = []
        else:
            self.beads = beads

    def __str__(self):
        return self.name + "," + "," + str(self.num) + "," + \
               str(self.chain_id) + "," + str(self.i_code) + \
               "," + basic_io.points_to_str(self.beads)

    def to_str(self):
        return str(self)

    def copy(self):
        return Residue(self.name, self.num, self.chain_id, self.i_code,
                       self.uuid, [b.copy() for b in self.beads])


class Chain(primitives.Chain):

    __slots__ = [
        "residues"
    ]

    def __init__(self, residues=None):
        super(self.__class__, self).__init__(residues)

    def __str__(self):
        s = ""
        for r in self.residues:
            s += r.to_str() + ";"
        return s

    def to_str(self):
        return str(self)

    def copy(self):
        residues = []
        for r in self.residues:
            residues.append(r.copy())
        return Chain(residues)


class Basepair(object):

    """
    A small container class to hold the "State" of a basepair for finding
    matches in the database. The critical features are:

    :params r: Reference Frame of basepair
    :type r: np.array
    :params d: Center of mass of basepair
    :type d: np.array
    :params sugars: C1` atom coords for both residues in basepair
    :type sugars: List of Np.Arrays

    Attributes
    ----------
    `r` : Np.Matrix
        Reference Frame of basepair
    `d` : Np.Array
        Center of mass of basepair
    `sugars` : List of Np.Arrays
        C1` atom coords for both residues in basepair

    """

    __slots__ = ["r", "d", "sugars", "res1", "res2", "name"]

    def __init__(self, r, d, sugars, res1, res2, name):
        self.r, self.d, self.sugars = r, d, sugars
        self.res1, self.res2 = res1, res2
        self.name = name

    def flip(self):
        self.r[1] = -self.r[1]
        self.r[2] = -self.r[2]

    def get_transforming_r_and_t(self, r, t, sugars):
        r1 = self.r
        r2 = r
        r_trans = util.unitarize(r1.T.dot(r2))
        t_trans = -t

        new_sugars_2 = np.dot(sugars, r_trans.T) + t_trans + self.d

        if sugars is not None:
            diff = (((self.sugars[0] - new_sugars_2[0]) +
                     (self.sugars[1] - new_sugars_2[1]))/2)
        else:
            diff = 0
        return r_trans, t_trans+diff

    def get_transforming_r_and_t_w_state(self, state):
		return self.get_transforming_r_and_t(state.r, state.d, state.sugars)

    def get_transformed_state(self, r, t):

        r_T = r.T

        new_r = util.unitarize(np.dot(self.r, r_T))
        new_sugars = np.dot(self.sugars, r_T) + t
        new_origin = np.dot(self.d, r_T) + t

        return new_r, new_origin, new_sugars

    def set(self, r, d, sug):
        self.r = r
        self.d = d
        self.sugars = sug

    def copy(self):
        """
        returns a deep copy of this BasepairState object
        """
        return Basepair(np.array(self.r, copy=True),
                        np.array(self.d, copy=True),
                       [coord[:] for coord in self.sugars],
                        self.res1,
                        self.res2,
                        self.name)

    def diff(self, state):
        diff  = util.distance(self.d, state.d)
        diff += self._rot_diff(state) * 2
        #diff += self._sugar_diff(state)
        return diff

    def _rot_diff(self, state):
        r_diff = util.matrix_distance(self.r, state.r)
        state.flip()
        r_diff_2 = util.matrix_distance(self.r, state.r)
        state.flip()
        if r_diff > r_diff_2:
            r_diff = r_diff_2
        return r_diff

    def _sugar_diff(self, state):
        diff_1 = util.distance(self.sugars[0], state.sugars[0]) + \
                 util.distance(self.sugars[1], state.sugars[1])
        diff_2 = util.distance(self.sugars[1], state.sugars[0]) + \
                 util.distance(self.sugars[0], state.sugars[1])
        if diff_1 > diff_2:
            diff_1 = diff_2
        return diff_1

    def to_str(self):
        """
        converts basepairstate into a string
        """

        s = basic_io.point_to_str(self.d) + ";" + \
            basic_io.matrix_to_str(self.r) + ";" +\
            basic_io.matrix_to_str(self.sugars)
        return s


class Structure(primitives.Structure):
    __slots__ = [
        "chains"
    ]

    def __init__(self, chains=None):
        super(self.__class__, self).__init__(chains)

    def copy(self):
        new_chains = [ c.copy() for c in self.chains ]
        return Structure(new_chains)


class Motif(primitives.RNAStructure):
    def __init__(self, struct=None, ends=None, name=None, end_names=None,
                 end_ids=None, score=0, block_end_add=0):
        self.structure = struct
        self.ends = ends
        self.name = name
        self.end_names = end_names
        self.end_ids = end_ids
        self.score = score
        self.block_end_add = 0

    def copy(self):
        struc = self.structure.copy()
        ends = []
        for end in self.ends:
            res1 = struc.get_residue(uuid=end.res1.uuid)
            res2 = struc.get_residue(uuid=end.res2.uuid)
            new_end = end.copy()
            new_end.res1 = res1
            new_end.res2 = res2
            ends.append(new_end)
        return Motif(struc, ends, self.name, self.end_names, self.score,
                     self.block_end_add)

    def beads_to_pdb(self):
        pass

    def to_motif(self):
        m = rm.manager.get_motif(name=self.name,
                                 end_name=self.ends[0].name,
                                 end_id=self.end_ids[0])

        motif.align_motif(self.ends[0], m.ends[0], m)
        return m

    def to_pdb(self, filename):
        return self.to_motif().to_pdb(filename)



def align_motif_state(ref_bp_state, org_state):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(org_state.ends[0])
    t += ref_bp_state.d
    r_T = r.T

    for i, s in enumerate(org_state.ends):
        new_r, new_d, new_sug = s.get_transformed_state(r, t)
        org_state.ends[i].set(new_r,new_d,new_sug)
