# Quick script to compute cylinder. I have seen that the following knots having linearized
# homology in negative degree don't have any interesting higher augmentations:
# m4_1, m6_1_2


import legendrian_links as ll
front = [1,2,2,0,0,1,1,0,2,1,1] # m6_1_2
pd = ll.PlatDiagram(n_strands=4, front_crossings=front)
pd.set_lch()
cyl = pd.lch_dga.get_cylinder()
