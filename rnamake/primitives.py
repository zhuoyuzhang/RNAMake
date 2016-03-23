

class Residue(object):
    def __init__(self, name, num, chain_id, i_code=""):
        self.name = name
        self.num = num
        self.chain_id = chain_id
        self.i_code = i_code

    def copy(self):
        pass

    def to_str(self):
        pass


class Basepair(object):
    def __init__(self, res1, res2, bp_type="c..."):
        self.res1 = res1
        self.res2 = res2
        self.bp_type = bp_type

    def residues(self):
        return [self.res1, self.res2]

    def name(self):
        str1 = self.res1.chain_id+str(self.res1.num)+str(self.res1.i_code)
        str2 = self.res2.chain_id+str(self.res2.num)+str(self.res2.i_code)

        if str1 < str2:
            return str1+"-"+str2
        else:
            return str2+"-"+str1

    def partner(self, res):
        if res == self.res1:
            return self.res2
        elif res == self.res2:
            return self.res1
        else:
            raise ValueError("residue not in basepair")

class Chain(object):
    def __init__(self, residues=None):
        self.residues = []
        if residues is not None:
            self.residues = residues

    def __len__(self):
        return len(self.residues)

    def subchain(self, start, end=None):
        if start < 0:
            raise ValueError("start cannot be less then 0")

        if end is None:
            end = len(self.residues)

        return Chain(self.residues[start:end])

    def copy(self):
        pass

    def first(self):
        return self.residues[0]

    def last(self):
        return self.residues[-1]

    def to_str(self):
        pass

    def sequence(self):
        seq = ""
        for r in self.residues:
            seq += r.name
        return seq

class Structure(object):
    def __init__(self, chains=None):
        self.chains = chains
        if self.chains is None:
            self.chains = []
        self.name = "N/A"

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        """
        find a residue based on residue num, chain_id, insert_code and uuid
        will return an error if more then one residue matches search to avoid
        confusion

        :param num: residue number
        :param chain id: what chain the residue belongs to
        :param i_code: the insertation code of the residue
        :param uuid: the unique indentifier that each residue is given

        :type num: int
        :type chain_id: str
        :type i_code: str
        :type uuid: uuid
        """

        found = []
        for c in self.chains:
            for r in c.residues:
                if uuid is not None and uuid != r.uuid:
                    continue
                if num is not None and num != r.num:
                    continue
                if i_code is not None and i_code != r.i_code:
                    continue
                if chain_id is not None and chain_id != r.chain_id:
                    continue
                found.append(r)

        if len(found) == 0:
            return None

        return found[0]

    def residues(self):
        residues = []
        for c in self.chains:
            residues.extend(c.residues)
        return residues

    def to_str(self):
        pass

    def copy(self):
        pass

    def sequence(self):
        sequences = [x.sequence() for x in self.chains]
        return "&".join(sequences)


class RNAStructure(object):
    def __init__(self):
        self.structure = Structure()

    def get_residue(self, num=None, chain_id=None, i_code=None, uuid=None):
        return self.structure.get_residue(num, chain_id, i_code, uuid)

    def get_basepair(self, bp_uuid=None, res1=None, res2=None, uuid1=None,
                     uuid2=None, name=None):
        """
        locates a Basepair object based on residue objects or uuids if nothing
        is supplied you will get back all the basepairs in the motif. The way
        to make sure you will only get one basepair back is to supply BOTH
        res1 and res2 OR uuid1 and uuid2, I have left it open like this
        because it is sometimes useful to get all basepairs that a single
        residue is involved

        :param res1: First residue
        :param res2: Second residue
        :param uuid1: First residue uuid
        :param uuid2: Second residue uuid

        :type res1: Residue object
        :type res2: Residue object
        :type uuid1: uuid object
        :type uuid2: uuid object
        """
        alt_name = None
        if name:
            name_spl = name.split("-")
            alt_name = name_spl[1] + "-" + name_spl[0]

        found = []
        for bp in self.basepairs:
            if bp_uuid is not None and bp_uuid != bp.uuid:
                continue
            if res1 is not None and (res1 != bp.res1 and res1 != bp.res2):
                continue
            if res2 is not None and (res2 != bp.res1 and res2 != bp.res2):
                continue
            if uuid1 is not None and \
               (uuid1 != bp.res1.uuid and uuid1 != bp.res2.uuid):
                continue
            if uuid2 is not None and \
               (uuid2 != bp.res1.uuid and uuid2 != bp.res2.uuid):
                continue
            if name is not None and \
               (name != bp.name() and alt_name != bp.name()):
                continue
            found.append(bp)
        return found

    def residues(self):
        return self.structure.residues()

    def chains(self):
        return self.structure.chains

    def sequence(self):
        return self.structure.sequence()

class Motif(RNAStructure):
    def __init__(self, struct=None, basepairs=None, ends=None):
        self.structure = struct
        self.basepairs = basepairs
        self.ends = ends
