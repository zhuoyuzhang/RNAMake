
import graph
import motif_merger
import motif_type
import motif_state

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
