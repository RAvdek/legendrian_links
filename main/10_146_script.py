import pickle
import legendrian_links as ll

# Twisted two-copy of 10_46
front = [1,5,9,13,17,7,8,6,7,9,10,8,9,11,12,10,11,13,14,12,13,15,16,14,15,3,4,2,3,5,6,4,5,7,8,6,7,9,10,8,9,11,12,10,11,13,14,12,13,13,14,12,13,11,12,10,11,9,10,8,9,9,10,8,9,13,14,12,13,3,4,2,3,7,8,6,7,11,12,10,11,15,16,14,15,1,5,9,13,18,18,17]

pd = ll.PlatDiagram(n_strands=20, front_crossings=front, n_copy=1, lazy_disks=False, lazy_lch=True, lazy_rsft=False)

with open('10_46_double.pk', 'wb') as f:
    pickle.dump(pd, f)
