import rnamake

##############################################################################
# LOADING MOTIFS                                                             #
##############################################################################
# motifs are the main objects of rnamake. They can be loaded in multiple ways
# from a preformated motif directory
m1 = rnamake.motif.Motif("resources/TWOWAY.1GID.4")

#from a pdb
m2 = rnamake.motif.Motif(pdb="resources/motif.pdb")

#from the resource manager
rm  = rnamake.resource_manager.ResourceManager()
m3 = rm.get_motif("HELIX.IDEAL")
m4 = rm.get_motif("TWOWAY.1GID.4")

#the resource manager can load any motif that is the libraries that rnamake
#contains. It can also get other resources that we will dicuss later.

##############################################################################
# BASIC MOTIF FEATURES                                                       #
##############################################################################
# every object in rnamake that contains structural information can printed
# out in pdb format, make sure to include to .pdb at the end
m1.to_pdb("test_output.pdb")

# In addition to pdb output, all objects including motif can be formatted to
# a text string. This allows for storage of motifs in text files and also
# in sqlite databases
s = m1.to_str()
mnew = rnamake.motif.str_to_motif(s)

#take a look at the format if your are curious
f = open("motif.txt", "w")
f.write(s)
f.close()

#possibly the most important feature of motifs is that they track the
#basepair ends of each motif in .ends variable. These ends are how motifs
#are combined together by overlapping these edges

print "m.ends => basepair ends of motifs"
for bp in m1.ends:
    print bp.name()

print "##############################################################################"
print "m.chains() => gets all residues in motif"
for c in m1.chains():
    print c


print "##############################################################################"
print "m.residues() => gets all residues in motif"
for r in m1.residues():
    print r
print "##############################################################################"


print "m.sequence() => ", m1.sequence()
print "m.secondary_structure() => ", m1.secondary_structure()
print "m.mtype => ", rnamake.motif_type.type_to_str(m4.mtype)