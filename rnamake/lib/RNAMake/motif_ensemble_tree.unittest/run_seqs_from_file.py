import rnamake.secondary_structure_tree as secondary_structure_tree
import subprocess

f = open("processed_chip_seq.txt")
lines = f.readlines()
f.close()

def extract_steps_from_ss_tree(ss_tree):
    node = None
    for n in ss_tree.nodes:
        if n.ss_type == "Bulge":
            node = n
            break

    steps = []
    while node != None:
        if len(node.children) == 0:
            break
        c = node.children[0]
        if c.ss_type != "Basepair":
            break
        step = list(c.bp_type)
        if step[0] == "T":
            step[0] = "U"
        if step[1] == "T":
            step[1] = "U"
        steps.append("".join(step))
        node = c
    return steps


def get_step_motifs(steps):
    motif_names = []
    for i in range(1,len(steps)):
        full_step = steps[i-1] + "=" + steps[i]
        motif_names.append(full_step)
    return motif_names



for l in lines:
    spl = l.split()
    chip_ss_tree = secondary_structure_tree.SecondaryStructureTree(spl[2] ,spl[1])
    chip_steps = extract_steps_from_ss_tree(chip_ss_tree)
    motif_names = get_step_motifs(chip_steps[1:])
    chip_str = ",".join(motif_names)
    for i in range (0, 2):
        subprocess.call("./a.out "+chip_str+"> out",shell=True)
        f = open("out")
        result = f.readline()
        f.close()
        print spl[0],result,
