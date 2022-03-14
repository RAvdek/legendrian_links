import legendrian_links as ll

# Twisted two-copy of 10_146_1
front = [1,0,0,5,9,13,7,8,6,7,9,10,8,9,9,10,8,9,5,6,4,5,7,8,6,7,11,12,10,11,13,14,12,13,3,4,2,3,1,2,0,1,9,10,8,9,11,12,10,11,5,6,4,5,3,4,2,3,3,4,2,3,7,8,6,7,11,12,10,11,1,2,0,1,5,6,4,5,9,10,8,9,13,14,12,13,1,2,0,1,3,4,2,3,5,6,4,5,13,14,12,13,11,12,10,11,5,6,4,5,7,8,6,7,1,5,9,13]

pd = ll.PlatDiagram(n_strands=16, front_crossings=front, n_copy=1, lazy_disks=False, lazy_lch=True, lazy_rsft=False)
pd.rsft_dga.pickle('10_46_double.pk')
