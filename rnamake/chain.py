import exceptions


class Chain(object):

    """
    Stored chain information from pdb file. Stores all residues in chain.
    Implementation is designed to be extremely lightweight. To connect residues
    into chains it highly adviced that you use
    :func:`connect_residues_into_chains`

    :param residues: the residues that are to be included in this chain
    :type residues: list of residue.Residue objects

    :attributes:

    `residues` : List of Residue objects
        The list of residues that belong to this chain will always be in
        5' to 3' order

    :examples:

    ..  code-block:: python

        >>> import rnamake.unittests.instances
        >>> c = rnamake.unittests.instances.chain()
        >>> c.first()
        <Residue('G103 chain A')>

        >>> c.last()
        <Residue('C260 chain A')>

        >>> len(c)
        157

        >>> cs = c.subchain(1, 10)
        >>> len(cs)
        9

        >>> cs.first()
        <Residue('A104 chain A')>

        >>> cs2 = c.subchain(start_res=c.residues[10], end_res=c.residues[15])
        >>> len(cs2)
        6

        >>> c.to_pdb_str()
        ATOM      1 O5'  G   A 103     -26.469 -47.756  84.669  1.00  0.00
        ATOM      2 C5'  G   A 103     -25.050 -47.579  84.775  1.00  0.00
        ATOM      3 C4'  G   A 103     -24.521 -48.156  86.068  1.00  0.00
        ATOM      4 O4'  G   A 103     -24.861 -49.568  86.118  1.00  0.00
        ATOM      5 C3'  G   A 103     -23.009 -48.119  86.281  1.00  0.00
        ATOM      6 O3'  G   A 103     -22.548 -46.872  86.808  1.00  0.00
        ATOM      7 C1'  G   A 103     -23.806 -50.289  86.732  1.00  0.00
        ATOM      8 C2'  G   A 103     -22.812 -49.259  87.269  1.00  0.00
        .
        .
        .

    """
    __slots__ = ["residues"]

    def __init__(self, residues=None):
        self.residues = residues
        if residues is None:
            self.residues = []

    def __repr__(self):
        if len(self.residues) == 0:
            return "<Chain( First: None\n\t  Last:  None\n\t  Size: 0)>"

        return "<Chain( First: %s\n\t  Last:  %s\n\t  Size: %s)>" %\
            (self.first(), self.last(), len(self.residues))

    def __len__(self):
        return len(self.residues)

    def first(self):
        """
        returns residue at 5' end of chain

        :returns: first residue in chain
        :rtype: residue.Residue

        :examples:

        ..  code-block:: python

            >>> import rnamake.unittests.instances
            >>> c = rnamake.unittests.instances.chain()
            >>> c.first()
            <Residue('G103 chain A')>

        """
        if len(self.residues) == 0:
            raise exceptions.ChainException("cannot call first there are no "
                                            "residues in chain")
        return self.residues[0]

    def last(self):
        """
        returns 3' end of chain
        """
        if len(self.residues) == 0:
            raise exceptions.ChainException("cannot call first there are no "
                                            "residues in chain")
        return self.residues[-1]

    def subchain(self, start=None, end=None, start_res=None, end_res=None):
        """
        Creates a new chain from a subsection of the residues in the current
        chain.

        :param start: start position in residues object list, default:None
        :param end: end position in residues object list, default:None
        :param start_res: The 5' residue of sub chain, default:None
        :param end_res:  The 3' resiude of sub chain, default:None

        :type start: int
        :type end: int

        :return: Chain object

        :examples:

        ..  code-block:: python

            >>> cs = c.subchain(1, 10)
            >>> len(cs)
            9

            >>> cs.first()
            <Residue('A104 chain A')>

            >>> cs2 = c.subchain(start_res=c.residues[10], end_res=c.residues[15])
            >>> len(cs2)
            6

        """

        if start_res is not None and end_res is not None:
            try:
                start = self.residues.index(start_res)
                end = self.residues.index(end_res)
            except:
                raise exceptions.ChainException("supplied start_res and end_res "
                                                "but they are not members of "
                                                "chain")

            if start > end:
                start, end = end, start
            end += 1

        elif start_res is not None and end_res is None:
            raise exceptions.ChainException("supplied start_res but not end_res")

        elif start_res is not None and start is not None:
            raise exceptions.ChainException("cannot supply start and start_res")

        if start < 0:
            raise exceptions.ChainException("start pos cannot be less then 0")

        if end is None:
            end = len(self.residues)

        return Chain(self.residues[start:end])

    def contain_res(self, r):
        if r in self.residues:
            return 1
        else:
            return 0

    def copy(self):
        """
        Creates a deepcopy of the this chain object.

        :return: Chain object
        """
        residues = [r.copy() for r in self.residues]
        return Chain(residues)

    def to_str(self):
        """
        Stringifes Chain object

        :return: str
        """

        s = ""
        for r in self.residues:
            s += r.to_str() + ";"
        return s

    def to_pdb_str(self, acount=1, return_acount=0, rnum=-1, chain_id=""):
        """
        creates a PDB string formatted verision of this Chain object.

        :param acount: current atom index, default: 1
        :param return_acount: final atom index after current atoms, default: 0
        :param rnum: starting residue number, default: -1
        :param chain_id: the chain id of the chain, i.e. "A", "B" etc

        :type  acount: int
        :type  return_acount: int
        :type  rnum: int
        :type  chain_id: str

        :return: str
        """

        s = ""
        for r in self.residues:
            r_str, acount = r.to_pdb_str(acount, 1, rnum, chain_id)
            if rnum != -1:
                rnum += 1
            s += r_str

        # TODO fix returns should only be one possibility
        if return_acount:
            return s, acount
        else:
            return s

    def to_pdb(self, fname="chain.pdb", rnum=-1, chain_id=""):
        """
        Writes a PDB string formmated verision of this Chain object to file

        :param fname: filename of output PDB file, default="chain.pdb"
        :param rnum: starting residue number, default: -1
        :param chain_id: the chain id of the chain, i.e. "A", "B" etc

        :type  fname: str
        :type  rnum: int
        :type  chain_id: str

        :return: None
        """
        f = open(fname, "w")
        f.write(self.to_pdb_str(rnum=rnum, chain_id=chain_id))
        f.close()


def connect_residues_into_chains(residues):
        """
        takes all residues and puts into the correct order in chains checking
        for physical connection between O5' and P atoms between residues

        :param residues: residue objects that belong in this structure
        :type residues: List of Residue objects

        :return: List of Chain objects
        """

        chains = []
        # sort residues so check residues for connection quicker as the next on
        # in the array will be closest to it by number
        residues.sort(key=lambda x: x.num)

        while True:
            current = None
            # find next 5' end, all chains go from 5' to 3'
            for i, r in enumerate(residues):
                five_prime_end = 1
                for j, r2 in enumerate(residues):
                    if r.connected_to(r2) == -1:
                        five_prime_end = 0
                        break
                if five_prime_end:
                    current = r
                    break
            if not current:
                break
            residues.remove(current)
            current_chain_res = []
            # extend chain until 3' end
            while current is not None:
                current_chain_res.append(current)
                found = 0
                for r in residues:
                    if current.connected_to(r) == 1:
                        current = r
                        found = 1
                        break
                if found:
                    residues.remove(current)
                else:
                    # no more residues to add, make chain object
                    chains.append(Chain(current_chain_res))
                    current = None

        return chains

