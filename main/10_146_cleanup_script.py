import algebra
import utils

LOG = utils.LOG
PK_NAME_IN = '10_46_double.pk'
PK_NAME_OUT = '10_146_double_fillna_0.pk'

dga = algebra.DGA.from_pickle(file_name=PK_NAME_IN, decompress=False)
dga.decompress_augmentations(fill_na=0)
dga.set_all_bilin()
dga.pickle(PK_NAME_OUT)
