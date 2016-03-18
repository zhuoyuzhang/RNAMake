import primitives
import basic_io
import util
import basic_io

import uuid
import numpy as np
import residue
import motif
import motif_type
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
        bp_state = Basepair(end.r(), end.d(), end.state().sugars, res1, res2)
        bp_state.uuid = end.uuid
        ends.append(bp_state)

    ends[0].res1.beads = []
    ends[0].res2.beads = []

    ms = Motif(struc, ends, end_ids=m.end_ids, name=m.name)
    ms.uuid = m.id
    return ms


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

    def to_str(self):
        return str(self.name) + "," + str(self.num) + "," + \
               str(self.chain_id) + "," + str(self.i_code) + \
               "," + basic_io.beads_to_str(self.beads)

    def copy(self):
        return Residue(self.name, self.num, self.chain_id, self.i_code,
                       self.uuid, [b.copy() for b in self.beads])


class Chain(primitives.Chain):

    __slots__ = [
        "residues"
    ]

    def __init__(self, residues=None):
        super(self.__class__, self).__init__(residues)

    def to_str(self):
        s = ""
        for r in self.residues:
            s += r.to_str() + ";"
        return s

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

    __slots__ = ["r", "d", "sugars", "res1", "res2", "name", "uuid"]

    def __init__(self, r, d, sugars, res1, res2):
        self.r, self.d, self.sugars = r, d, sugars
        self.res1, self.res2 = res1, res2
        self.name = self._get_name()

    def _get_name(self):
        str1 = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        str2 = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        if str1 < str2:
            return str1+"-"+str2
        else:
            return str2+"-"+str1

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
                        self.res2)

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

        s = self.name+ ";" + basic_io.point_to_str(self.d) + ";" + \
            basic_io.matrix_to_str(self.r) + ";" +\
            basic_io.matrix_to_str(self.sugars) + ";"

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

    def to_str(self):
        s = ""
        for c in self.chains:
            s += c.to_str() + "|"
        return s


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
        self.mtype = motif_type.UNKNOWN
        self.path = ""
        self.uuid = uuid.uuid1()

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
        return Motif(struc, ends, self.name, self.end_names, self.end_ids,
                     self.score, self.block_end_add)

    def beads(self):
        beads = []
        for r in self.residues():
            beads.extend(r.beads)
        return beads

    def beads_to_pdb(self, filename):
        beads = self.beads()
        basic_io.beads_to_pdb(filename, beads)

    def to_motif(self):
        m = rm.manager.get_motif(name=self.name,
                                 end_name=self.ends[0].name,
                                 end_id=self.end_ids[0])

        motif.align_motif(self.ends[0], m.ends[0], m)
        return m

    def to_pdb(self, filename):
        return self.to_motif().to_pdb(filename)

    def to_str(self):
        s = self.path + "&" + self.name + "&" + str(self.score) + "&" +  \
            str(self.block_end_add) + "&" + \
            str(self.mtype) + "&" + self.structure.to_str() + "&"
        for end in self.ends:
            s += end.to_str() + "@"
        s += "&"
        for end_id in self.end_ids:
            s += end_id + " "
        s += "&"
        return s

    def copy_uuids_from(self, m):
        self.uuid = m.uuid
        for i, end in enumerate(m.ends):
            self.ends[i].uuid = end.uuid
        res = self.residues()
        for i, r in enumerate(m.residues()):
            res[i] = r.uuid

    def get_end_state(self, name):
        for i, end in enumerate(self.ends):
            if i == self.block_end_add:
                continue
            if end.name == name:
                return end
        raise ValueError("cannot get end state of name: " + name)


def str_to_residue(s):
    spl = s.split(",")
    beads = []
    for i in range(4, len(spl)-1):
        b = residue.str_to_bead(spl[i])
        beads.append(b)
    return Residue(spl[0], int(spl[1]), spl[2], spl[3], None, beads)


def str_to_chain(s):
    spl = s.split(";")
    residues = []
    for r_str in spl[:-1]:
        r = str_to_residue(r_str)
        residues.append(r)
    return Chain(residues)


def str_to_structure(s):
    spl = s.split("|")
    chains = []
    for c_str in spl[:-1]:
        c = str_to_chain(c_str)
        chains.append(c)
    return Structure(chains)


def str_to_motif(s):
    spl = s.split("&")
    m = Motif()
    m.path = spl[0]
    m.name = spl[1]
    m.score = float(spl[2])
    m.block_end_add = int(spl[3])
    m.mtype = int(spl[4])
    m.structure = str_to_structure(spl[5])
    m.ends = []
    ends_str = spl[6].split("@")
    for end_str in ends_str[:-1]:
        bp_spl = end_str.split(";")
        res_spl = bp_spl[0].split("-")
        res1_id, res1_num = res_spl[0][0], int(res_spl[0][1:])
        res2_id, res2_num = res_spl[1][0], int(res_spl[1][1:])
        res1 = m.get_residue(num=res1_num, chain_id=res1_id)
        res2 = m.get_residue(num=res2_num, chain_id=res2_id)
        d = basic_io.str_to_point(bp_spl[1])
        r = basic_io.str_to_matrix(bp_spl[2])
        sugars = basic_io.str_to_matrix(bp_spl[3])
        end = Basepair(r, d, sugars, res1, res2)
        m.ends.append(end)
    m.end_ids = spl[7].split()

    return m


def align_motif_state(ref_bp_state, org_state):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(org_state.ends[0])
    t += ref_bp_state.d
    r_T = r.T

    for i, s in enumerate(org_state.ends):
        new_r, new_d, new_sug = s.get_transformed_state(r, t)
        org_state.ends[i].set(new_r, new_d, new_sug)

    for r in org_state.residues():
        for b in r.beads:
            b.center = np.dot(b.center, r_T) + t


def get_aligned_motif_state(ref_bp_state, cur_state, org_state):
    r, t = ref_bp_state.get_transforming_r_and_t_w_state(org_state.ends[0])
    t += ref_bp_state.d
    r_T = r.T

    for i, s in enumerate(org_state.ends):
        new_r, new_d, new_sug = s.get_transformed_state(r, t)
        cur_state.ends[i].set(new_r, new_d, new_sug)

    for r in cur_state.residues():
        for b in r.beads:
            b.center = np.dot(b.center, r_T) + t



