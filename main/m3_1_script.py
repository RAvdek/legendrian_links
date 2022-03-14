import legendrian_links as ll

# Twisted two-copy of 10_146_1
front = [1,1,1]

pd = ll.PlatDiagram(n_strands=4, front_crossings=front, n_copy=1, lazy_disks=False, lazy_lch=True, lazy_rsft=False)
pd.rsft_dga.pickle('m3_1.pk')
