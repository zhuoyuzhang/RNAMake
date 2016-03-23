
import graph
import motif_merger
import motif_type
import motif_state
import secondary_structure

class MotiftoSecondaryStructure(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.chains = []
        self.open_chains = []
        self.seen_res = {}
        self.seen_bp = {}

    def _get_next_chain(self, motif):
        best_score = -1

        for c in self.chains:
            score = 0
            for r in c.residues:
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in self.seen_bp:
                        score += 1
            if score > best_score:
                best_score = score

        best_chains = []
        for c in self.chains:
            score = 0
            for r in c.residues:
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in self.seen_bp:
                        score += 1
            if score == best_score:
                best_chains.append(c)

        best_chain = None
        best_score = 10000
        for c in best_chains:
            pos = 1000
            for i, r in enumerate(c.residues):
                bps = motif.get_basepair(res1=r)
                for bp in bps:
                    if bp in self.seen_bp:
                        pos = i
                        break
            if pos < best_score:
                best_score = pos
                best_chain = c

        return best_chain

    def _setup_basepairs_and_ends(self, struct, motif):
        ss_bps = []
        for bp in self.seen_bp.keys():
            res1 = struct.get_residue(uuid=bp.res1.uuid)
            res2 = struct.get_residue(uuid=bp.res2.uuid)
            ss_bps.append(secondary_structure.Basepair(res1, res2, bp.uuid))
        ss_ends = []
        for end in motif.ends:
            res1 = struct.get_residue(uuid=end.res1.uuid)
            res2 = struct.get_residue(uuid=end.res2.uuid)
            end_bp = None
            for bp in ss_bps:
                if bp.res1 == res1 and bp.res2 == res2:
                    end_bp = bp
                    break
            if end_bp is None:
                motif.to_pdb("test.pdb")
                raise ValueError("did not properly find end in generating ss")
            ss_ends.append(end_bp)

        return secondary_structure.Motif(struct, ss_bps, ss_ends)

    def to_secondary_structure(self, motif):
        saved_bp = None
        ss_chains = []

        self.chains = motif.chains()[::]
        self.open_chains = [self.chains.pop(0)]

        while len(self.open_chains) > 0:
            c = self.open_chains.pop(0)
            ss_res = []

            for r in c.residues:
                ss = "."
                bps = motif.get_basepair(res1=r)
                is_bp = 0
                for bp in bps:
                    partner_res = bp.partner(r)
                    passes = 1
                    saved_bp = None
                    if passes:
                        saved_bp = bp
                        if   bp not in self.seen_bp and \
                              r not in self.seen_res and \
                              partner_res not in self.seen_res:
                            self.seen_res[r] = 1
                            ss = "("
                        elif partner_res in self.seen_res:
                            if self.seen_res[partner_res] > 1:
                                ss = "."
                            else:
                                ss = ")"
                                self.seen_res[r] = 1
                                self.seen_res[partner_res] += 1
                                break
                    elif r not in self.seen_res:
                        ss = "."

                if saved_bp is not None:
                    self.seen_bp[saved_bp] = 1

                ss_res.append(secondary_structure.Residue(r.name, ss, r.num,
                                                          r.chain_id, r.uuid, r.i_code))
            ss_chains.append(secondary_structure.Chain(ss_res))
            best_chain = self._get_next_chain(motif)

            if best_chain is None:
                break
            self.chains.remove(best_chain)
            self.open_chains.append(best_chain)

        struct = secondary_structure.Structure(chains=ss_chains)
        m = self._setup_basepairs_and_ends(struct, motif)

        return m


class MotifStateMerger(object):
    def __init__(self):
        self.all_bps = {}
        self.motifs = {}
        self.res_overrides = {}
        self.bp_overrides = {}

        self.rebuild_structure = 1
        self.chain_graph = graph.GraphStatic()

    def add_state(self, ms, ms_end=None, parent=None, parent_end=None):
        new_chains = [c.subchain(0) for c in ms.chains()]

        for c in new_chains:
            data = motif_merger.ChainNodeData(c, ms.uuid)
            self.chain_graph.add_data(data, n_children=2, orphan=1)

        for bp in ms.basepairs:
            self.all_bps[bp.uuid] = bp

        if parent is not None:
            self._link_motifs(parent, ms, parent_end, ms_end)

        self.motifs[ms.uuid] = ms

    def _link_motifs(self, m1, m2, m1_end, m2_end):
        m1_end_nodes = self._get_end_nodes(self.chain_graph.nodes, m1_end)
        m2_end_nodes = self._get_end_nodes(self.chain_graph.nodes, m2_end)

        if m2.mtype == motif_type.HELIX and m1.mtype != motif_type.HELIX:
            self._link_chains(m1_end_nodes, m2_end_nodes)
            self.bp_overrides[m2_end.uuid] = m1_end.uuid
        else:
            self._link_chains(m2_end_nodes, m1_end_nodes)
            self.bp_overrides[m1_end.uuid] = m2_end.uuid

    def _link_chains(self, dominant_nodes, auxiliary_nodes):
        if dominant_nodes[0] == dominant_nodes[1]:
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[0], 1, 0)
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[1], 0, 1)

        elif auxiliary_nodes[0] == auxiliary_nodes[1]:
            print "did this happen check to make sure it worked!!!"
            self._connect_chains(dominant_nodes[1], auxiliary_nodes[0], 1, 0)
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[0], 0, 1)

        else:
            self._connect_chains(dominant_nodes[1], auxiliary_nodes[0], 1, 0)
            self._connect_chains(dominant_nodes[0], auxiliary_nodes[1], 0, 1)

    def _connect_chains(self, d_node, a_node, d_i, a_i):
        if a_i == 0:
            a_node.data.prime5_override = 1
            self.res_overrides[a_node.data.c.first().uuid] = \
                d_node.data.c.last().uuid
        else:
            a_node.data.prime3_override = 1
            self.res_overrides[a_node.data.c.last().uuid] = \
                d_node.data.c.first().uuid
        self.chain_graph.connect(d_node.index, a_node.index, d_i, a_i)

    def _get_end_nodes(self, nodes, end):
        end_nodes = [None, None]

        for n in nodes:
            for r in end.residues():
                if n.data.c.first().uuid == r.uuid and end_nodes[0] == None:
                    end_nodes[0] = n
                elif n.data.c.first().uuid == r.uuid:
                    raise ValueError("cannot build chain map two residues are assigned"
                                     "to 5' chain")

                if n.data.c.last().uuid == r.uuid and end_nodes[1] == None:
                    end_nodes[1] = n
                elif n.data.c.last().uuid == r.uuid:
                    raise ValueError("cannot build chain map two residues are assigned"
                                     "to 3' chain")

        if end_nodes[0] == None or end_nodes[1] == None:
            raise ValueError("did not build map properly, both chains are not found")

        return end_nodes

    def merge(self, mst):
        for i, n in enumerate(mst.tree.nodes):
            if i == 0:
                self.add_state(n.data.cur_state)
            else:
                self.add_state(n.data.cur_state,
                               n.data.cur_state.ends[0],
                               n.parent.data.cur_state,
                               n.parent.data.cur_state.ends[n.parent_end_index()])

    def _build_structure(self):
        starts = []

        for n in self.chain_graph.nodes:
            if n.available_pos(0):
                starts.append(n)

        chains = []
        for n in starts:
            res = []
            cur = n
            while cur is not None:
                res.extend(cur.data.included_res())
                if cur.available_pos(1):
                    cur = None
                else:
                    con = cur.connections[1]
                    cur = con.partner(cur.index)

            c = motif_state.Chain(res)
            chains.append(c)

        rna_structure = motif_state.RNAStructure()
        rna_structure.structure.chains = chains
        res = rna_structure.residues()
        uuids = [r.uuid for r in res]

        current_bps = []
        for bp in self.all_bps.values():
            if bp.res1.uuid in uuids and bp.res2.uuid in uuids:
                current_bps.append(bp)
        rna_structure.basepairs = current_bps
        rna_structure.ends = self.ends_from_basepairs(rna_structure.structure, current_bps)

        return rna_structure

    def ends_from_basepairs(self, structure, basepairs):
        chain_ends = []
        for c in structure.chains:
            chain_ends.append(c.first())
            if len(c) > 1:
                chain_ends.append(c.last())

        ends = []
        for bp in basepairs:
            if not (self.gu_bp(bp) or self.wc_bp(bp)):
                continue

            if bp.res1 in chain_ends and bp.res2 in chain_ends:
                ends.append(bp)

        return ends

    def sequence(self):
        rna_struc = self._build_structure()
        return rna_struc.sequence()

    def secondary_structure(self):
        converter = MotiftoSecondaryStructure()
        ss = converter.to_secondary_structure(self._build_structure())
        print ss


    def wc_bp(self, bp):
        bp_str = bp.res1.name[0] + bp.res2.name[0]
        wc = "GC,CG,AU,UA".split(",")
        if bp_str in wc:
            return 1
        else:
            return 0

    def gu_bp(self, bp):
        bp_str = bp.res1.name[0] + bp.res2.name[0]
        if bp_str == "GU" or bp_str == "UG":
            return 1
        else:
            return 0

