import motif_library
import rnamake.motif_factory as motif_factory
import rnamake.motif_ensemble as motif_ensemble
import rnamake.motif as motif
import rnamake.motif_type as motif_type
import rnamake.util as util
import rnamake.cluster as cluster
import rnamake.settings as settings
import rnamake.sqlite_library as sqlite_library
import math
import subprocess
import os
import numpy as np

def test_align():
    mlib = sqlite_library.MotifLibrarySqlite(libname="ideal_helices")
    m1 = mlib.get_motif("HELIX.IDEAL.12")
    m2 = m1.copy()

    aligned = motif.get_aligned_motif(m1.ends[1], m2.ends[0], m2)

    m1.to_pdb("test_1.pdb")
    aligned.to_pdb("test_2.pdb")

def setup_start_motif():
    start_path = settings.RESOURCES_PATH + "start"
    m = motif_factory.factory.motif_from_file(start_path)
    s = m.to_str()
    m_path = settings.MOTIF_DIRS + "ref.motif"
    print m_path
    f = open(m_path, "w")
    f.write(s)
    f.close()

    m2 = motif.file_to_motif(m_path)

class SSandSeqCluster(object):

    def __init__(self, end_id):
        self.end_id = end_id
        self.motif_and_ends = []

    def motif_matches_end(self, m, ei):
        class MotifandEnd(object):
            def __init__(self, m, ei):
                self.motif, self.end_index = m, ei

        if self.end_id == m.end_ids[ei]:
            self.motif_and_ends.append(MotifandEnd(m, ei))
            return 1
        else:
            return 0

    def motif_matches(self, m):
        for i, end_id in enumerate(m.end_ids):
            r = self.motif_matches_end(m, i)
            if r == 1:
                return 1
        return 0


class BuildSqliteLibraries(object):

    def build_ideal_helices(self):
        mlib = motif_library.ideal_helix_lib()
        aligned_motifs = []
        path = settings.RESOURCES_PATH +"/motif_libraries_new/ideal_helices.db"
        data = []
        rdata = []
        keys = ['data', 'name', 'end_name', 'end_id', 'id']
        count = 0
        for i, m in enumerate(mlib.motifs()):
            #if m.name == "HELIX.IDEAL":
            #    aligned_motif = motif_factory.factory.align_motif_to_common_frame(m, 0)
            #else:
            aligned_motif = motif_factory.factory.align_motif_to_common_frame(m, 1)

            data.append([aligned_motif.to_str(), aligned_motif.name,
                         aligned_motif.ends[0].name(), aligned_motif.end_ids[0], count])

            aligned_motif = motif_factory.factory.can_align_motif_to_end(m, 0)
            aligned_motif = motif_factory.factory.align_motif_to_common_frame(aligned_motif, 0)
            aligned_motif.name = aligned_motif.name
            motif_factory.factory._setup_secondary_structure(aligned_motif)
            rdata.append([aligned_motif.to_str(), aligned_motif.name,
                         aligned_motif.ends[0].name(), aligned_motif.end_ids[0], count])
            count += 1

        sqlite_library.build_sqlite_library(path, data, keys, 'id')
        path = settings.RESOURCES_PATH +"/motif_libraries_new/ideal_helices_reversed.db"
        sqlite_library.build_sqlite_library(path, rdata, keys, 'id')

        mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
        m = mlib.get(name="HELIX.IDEAL.3")
        s = m.to_str()
        path = settings.RESOURCES_PATH + "/motifs/base.motif"
        f = open(path, "w")
        f.write(s)
        f.close()

    def build_basic_libraries(self):

        types = [motif_type.TWOWAY, motif_type.NWAY, motif_type.HAIRPIN,
                 motif_type.TCONTACT]
        #types = [motif_type.HAIRPIN]
        for t in types:
            count = 0
            mlib = motif_library.MotifLibrary(t)
            mlib.load_all()

            succeses = []
            for k, m in enumerate(mlib.motifs()):
                if m is None:
                    continue

                if t == motif_type.HAIRPIN:
                    m.block_end_add = -1

                if t != motif_type.HAIRPIN and len(m.ends) == 1:
                    continue
                elif len(m.ends) == 0:
                    continue

                for ei in range(len(m.ends)):
                    m_added = motif_factory.factory.can_align_motif_to_end(m, ei)
                    if m_added is None:
                        continue
                    else:
                        succeses.append([m_added, ei])

            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']
            for i, s in enumerate(succeses):
                m, ei = s
                m_added = motif_factory.factory.align_motif_to_common_frame(m, ei)
                data.append([m_added.to_str(), m_added.name,
                            m_added.ends[0].name(), m_added.end_ids[0], i])

            path = settings.RESOURCES_PATH +"/motif_libraries_new/"+\
                   motif_type.type_to_str(t).lower()+".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_helix_ensembles(self):
        helix_mlib = motif_library.MotifLibrary(motif_type.HELIX)
        helix_mlib.load_all()

        clusters = []
        steps = []

        kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
        kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

        for m in helix_mlib.motifs():
            spl = m.name.split(".")
            if spl[1] == "IDEAL" or spl[1] == "LE":
                continue
            for i in range(len(m.basepairs)-1):
                if m.basepairs[i].bp_type != "cW-W" or m.basepairs[i+1].bp_type != "cW-W":
                    continue
                bps = [m.basepairs[i], m.basepairs[i+1]]

                res = []
                for bp in bps:
                    for r in bp.residues():
                        if r not in res:
                            res.append(r)
                if len(res) != 4:
                    continue
                try:
                    m_bps = motif_factory.factory.motif_from_bps(bps)
                except:
                    continue
                if len(m_bps.end_ids[0]) != 11:
                    continue
                matched = 0
                for c in clusters:
                    if c.motif_matches(m_bps):
                        matched = 1

                if not matched:
                    clusters.append(SSandSeqCluster(m_bps.end_ids[0]))
                    clusters[-1].motif_matches(m_bps)
                    if m_bps.end_ids[0] != m_bps.end_ids[1]:
                        clusters.append(SSandSeqCluster(m_bps.end_ids[1]))
                        clusters[-1].motif_matches(m_bps)

        mes_keys = ['data', 'name', 'id']
        mes_data = []

        all_mes_keys = ['data', 'name', 'id']
        all_mes_data = []
        all_count = 1

        motif_data = []
        motif_keys = ['data', 'name', 'end_name', 'end_id', 'id']
        count = 0
        for c in clusters:
            all_count += 1
            spl = c.end_id.split("_")
            aligned_motifs = []
            for i, m_and_e in enumerate(c.motif_and_ends):
                m, ei = m_and_e.motif, m_and_e.end_index
                try:
                    m_a = motif_factory.factory.can_align_motif_to_end(m, ei)
                except:
                    continue
                if m_a is None:
                    continue
                m_a = motif_factory.factory.align_motif_to_common_frame(m_a, ei)

                #catch basepairs that have 5' - 5' interactions or 3' - 3'
                fail = 0
                for bp in m_a.basepairs:
                    vec1 = bp.res1.get_atom("C2'").coords - bp.res1.get_atom("O4'").coords
                    vec2 = bp.res2.get_atom("C2'").coords - bp.res2.get_atom("O4'").coords
                    if np.dot(vec1, vec2) > 2:
                        fail = 1
                        break
                if fail:
                    continue

                aligned_motifs.append(m_a)

            motifs_and_energies = []
            """sum = open(c.end_id + ".scores", "w")
            for m in aligned_motifs:
                m.to_pdb("motif.pdb")
                subprocess.call("rna_score -s motif.pdb", shell=True)
                f = open("default.out")
                lines = f.readlines()
                f.close()
                energy_spl = lines[3].split()
                os.remove("default.out")
                motifs_and_energies.append([m, float(energy_spl[1])])
                sum.write(energy_spl[1] + "\n")

            sum.close()"""

            f = open(c.end_id + ".scores")
            lines = f.readlines()
            f.close()

            #best 0.65
            m_clusters = cluster.cluster_motifs(aligned_motifs, 0.65)
            print "start, ", c.end_id
            #m_clusters = cluster.cluster_motifs_3(scored_motifs, 0.65, 150)
            #print "last", len(m_clusters)
            clustered_motifs = []
            energies = []
            for j, c_motifs in enumerate(m_clusters):
                m = c_motifs.motifs[0]
                m.mtype = motif_type.HELIX
                m.name = spl[0][0]+spl[2][1]+"="+spl[0][1]+spl[2][0] + "." + str(j)
                motif_data.append([m.to_str(), m.name, m.ends[0].name(), c.end_id, count])
                count += 1
                clustered_motifs.append(m)

                pop = float(len(c_motifs.motifs)) / float(len(aligned_motifs))
                energy = -kBT*math.log(pop)
                energies.append(energy)

            me = motif_ensemble.MotifEnsemble()
            me.setup(c.end_id, clustered_motifs, energies)
            mes_data.append([me.to_str(), me.id, count])

            motif = me.members[0].motif
            motif.name = spl[0][0]+spl[2][1]+"="+spl[0][1]+spl[2][0]

            motif_data.append([motif.to_str(), motif.name, motif.ends[0].name(),
                               me.id, count])
            count += 1



        path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/bp_steps.db"
        sqlite_library.build_sqlite_library(path, mes_data, mes_keys, 'id')
        path = settings.RESOURCES_PATH +"/motif_libraries_new/bp_steps.db"
        sqlite_library.build_sqlite_library(path, motif_data, motif_keys, 'id')

    def build_motif_state_libraries(self):
        for libname in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']

            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()
            motif_states = []
            names = []
            for i, m in enumerate(mlib.all()):
                ms = m.get_state()
                data.append([ms.to_str(), ms.name,
                             ms.end_names[0], ms.end_ids[0], i])



            path = settings.RESOURCES_PATH + "/motif_state_libraries/" + libname + ".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_unique_twoway_library(self):
        mlib = sqlite_library.MotifSqliteLibrary("twoway")
        mlib.load_all()
        clusters = cluster.cluster_motifs(mlib.all(), 9.0)
        motif_arrays = []
        motif_array_names = []

        data = []
        keys = ['data', 'name', 'end_name', 'end_id', 'id']

        mes_keys = ['data', 'name', 'id']
        mes_data = []
        f = open("sim_list_new", "w")

        count = 0
        for i, c in enumerate(clusters):
            lowest = c.motifs[0]
            for m in c.motifs:
                if lowest.score > m.score:
                    lowest = m

            if lowest.name == "TWOWAY.1GID.6" or \
               lowest.name == "TWOWAY.1GID.2" or \
               lowest.name == "TWOWAY.2GDI.4":
                continue
            count += 1

            f.write(lowest.name + "," + lowest.ends[0].name() + " | ")

            for m in c.motifs:
                f.write(m.name + "," + m.ends[0].name() + " | ")
            f.write("\n")
            #print i, lowest.score, lowest.name
            #lowest.to_pdb("m."+str(i)+".pdb")

            #motif_arrays.append(motif.MotifArray(c.motifs))
            #motif_array_names.append(lowest.name)

            data.append([lowest.to_str(), lowest.name,
                        lowest.ends[0].name(), lowest.end_ids[0], count])


            #remove duplicate sequences
            motifs = []
            end_ids = {}
            for m in c.motifs:
                if m.end_ids[0] not in end_ids:
                    end_ids[m.end_ids[0]] = m
                    motifs.append(m)
                    continue

                org_m = end_ids[m.end_ids[0]]
                if m.score < org_m.score:
                    motifs.remove(org_m)
                    motifs.append(m)
                    end_ids[m.end_ids[0]] = m

            scores = [1 for x in motifs]

            me = motif_ensemble.MotifEnsemble()
            me.setup(lowest.name, motifs, scores)
            mes_data.append([me.to_str(), me.id, count])
        f.close()

        path = settings.RESOURCES_PATH +"/motif_libraries_new/unique_twoway.db"
        sqlite_library.build_sqlite_library(path, data, keys, 'id')

        #path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/twoway_clusters.db"
        #sqlite_library.build_sqlite_library(path, mes_data, mes_keys, 'id')

    def build_ss_and_seq_libraries(self):
        libnames = ["twoway", "tcontact", "hairpin", "nway"]
        #libnames = ["twoway"]

        for libname in libnames:
            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()

            clusters = []

            mes = []
            mes_names = []

            motifs = []

            kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
            kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)

            for m in mlib.all():
                matched = 0
                for c in clusters:
                    if c.motif_matches_end(m, 0):
                        matched = 1

                if not matched:
                    clusters.append(SSandSeqCluster(m.end_ids[0]))
                    clusters[-1].motif_matches_end(m, 0)

            keys = ['data', 'name', 'id']
            data = []

            for i, c in enumerate(clusters):

                all_motifs = []
                for m_e in c.motif_and_ends:
                    all_motifs.append(m_e.motif)
                if libname != "twoway":
                    energies = [1 for x in all_motifs]
                    me = motif_ensemble.MotifEnsemble()
                    me.setup(c.end_id, all_motifs, energies)
                    mes.append(me)
                    motifs.append(me.members[0].motif)
                    motifs[-1].name = me.id
                    mes_names.append(me.id)
                    continue
                m_clusters = cluster.cluster_motifs(all_motifs)
                clustered_motifs = []
                energies = []
                for j, c_motifs in enumerate(m_clusters):
                    clustered_motifs.append(c_motifs.motifs[0])
                    pop = float(len(c_motifs.motifs)) / float(len(all_motifs))
                    energy = -kBT*math.log(pop)
                    energies.append(energy)

                me = motif_ensemble.MotifEnsemble()
                me.setup(c.end_id, clustered_motifs, energies)
                data.append([me.to_str(), me.id, i])

            print libname, len(mlib.all()), len(clusters)

            path = settings.RESOURCES_PATH +"/motif_ensemble_libraries/"+libname+".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_motif_ensemble_state_libraries(self):

        for libname in sqlite_library.MotifEnsembleSqliteLibrary.get_libnames().keys():
            print libname
            me_lib = sqlite_library.MotifEnsembleSqliteLibrary(libname)
            me_lib.load_all()

            keys = ['data', 'name', 'id']
            data = []

            for i, me in enumerate(me_lib.all()):
                mse = me.get_state()
                #print len(me.members)
                data.append([mse.to_str(), mse.id, i])

            path = settings.RESOURCES_PATH +"/motif_state_ensemble_libraries/"+libname+".db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

    def build_trimmed_ideal_helix_library(self):
        for libname in sqlite_library.MotifSqliteLibrary.get_libnames().keys():
            if libname != "ideal_helices":
                continue
            data = []
            keys = ['data', 'name', 'end_name', 'end_id', 'id']

            mlib = sqlite_library.MotifSqliteLibrary(libname)
            mlib.load_all()
            motif_states = []
            names = []
            j = 0
            for i, m in enumerate(mlib.all()):
                name = m.name
                spl = name.split(".")
                if len(spl) < 3:
                    continue
                num = int(spl[2])
                if num < 2:
                    continue
                ms = m.get_state()
                data.append([ms.to_str(), ms.name,
                             ms.end_names[0], ms.end_ids[0], j])
                j += 1



            path = settings.RESOURCES_PATH + "/motif_state_libraries/" + libname + "_min.db"
            sqlite_library.build_sqlite_library(path, data, keys, 'id')

#setup_start_motif()
builder = BuildSqliteLibraries()

builder.build_ideal_helices()
#builder.build_basic_libraries()
#builder.build_helix_ensembles()
#builder.build_ss_and_seq_libraries()
#builder.build_unique_twoway_library()
#builder.build_motif_state_libraries()
#builder.build_motif_ensemble_state_libraries()
#builder.build_trimmed_ideal_helix_library()

#mlib = sqlite_library.MotifSqliteLibrary("ideal_helices")
#m = mlib.get("HELIX.IDEAL")

#me_lib = sqlite_library.MotifStateEnsembleSqliteLibrary("bp_steps")
#me = me_lib.get("GG_LL_CC_RR")
#print len(me.members)

#print ref_motif.name



#builder.build_ideal_helices()
