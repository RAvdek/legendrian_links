import legendrian_links as ll
import utils

LOG = utils.LOG
PK_NAME = 'm3_1.pk'

# Twisted two-copy of 10_146_1
front = [1,1,1]

pd = ll.PlatDiagram(n_strands=16, front_crossings=front, n_copy=1, lazy_disks=False, lazy_lch=True, lazy_rsft=True)
pd.set_rsft(lazy_augs=True, lazy_bilin=True)
dga = pd.rsft_dga
LOG.info("Pickling after initial setup")
dga.pickle(PK_NAME)
dga.set_aug_data()
LOG.info("Pickling after set_aug_data")
dga.pickle(PK_NAME)
dga.set_augmentations()
LOG.info("Pickling after set_augmentations")
dga.pickle(PK_NAME)
dga.set_all_bilin()
