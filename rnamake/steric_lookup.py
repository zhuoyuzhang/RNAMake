import util
import settings

class StericLookup(object):
    def __init__(self, m=None):
        self.bhash = {}
        self.grid_size = 0.50
        self._setup_additions(self.grid_size, 5)
        if m is not None:
            if len(m.beads) == 0:
                m.get_beads(m.ends)
            self.add_beads(m.beads)

    def _setup_additions(self, grid_size, search_radius):
        self.additions = []
        add = []
        for i in range(1, search_radius+1):
            add.append(-i*grid_size)
        add.append(0)
        for i in range(1, search_radius+1):
            add.append(i*grid_size)

        for x in add:
            for y in add:
                for z in add:
                    dist = util.distance([x, y, z], [0, 0, 0])
                    if dist > 3.5:
                        continue
                    self.additions.append([x, y, z])

    def _generate_keys(self, center):
        rounded = []
        for d in center:
            pos = round(d / self.grid_size)*self.grid_size
            rounded.append(pos)
        final = []
        for a in self.additions:
            final.append([a[0]+rounded[0], a[1]+rounded[1], a[2]+rounded[2]])
        final_strs = []
        for a in final:
            poss = []
            for pos in a:
                spos = str(pos)
                if spos == "-0.0":
                    spos = "0.0"
                poss.append(spos)
            final_strs.append(" ".join(poss))
        return final_strs

    def add_centers(self, centers):
        for c in centers:
            keys = self._generate_keys(c)
            for key in keys:
                if key not in self.bhash:
                    self.bhash[key] = 1
                else:
                    self.bhash[key] += 1

    def add_beads(self, beads):
        centers = []
        for b in beads:
            if b.btype == 0:
                continue
            centers.append(b.center)
        self.add_centers(centers)

    def clash(self, centers):
        poss = ["", "", ""]
        for j, c in enumerate(centers):
            for i, d in enumerate(c):
                pos = round(d / self.grid_size)*self.grid_size
                spos = str(pos)
                if spos == "-0.0":
                    spos = "0.0"
                poss[i] = spos
            key = " ".join(poss)
            if key in self.bhash:
                return 1
        return 0
