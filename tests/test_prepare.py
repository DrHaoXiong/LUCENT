import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

from lucent.prepare import prepare_h5ad_section
from lucent.spatial import SpatialData, prepare_hypothesis


@unittest.skipUnless(importlib.util.find_spec('anndata'), 'optional anndata not installed')
class H5adPreparationTests(unittest.TestCase):
    def test_whole_section_and_saved_graph_are_retained(self):
        import anndata as ad
        import pandas as pd
        from scipy.sparse import csr_matrix
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            counts=np.array([[0,1,7],[1,0,8],[2,2,9],[3,1,10]],dtype=np.int32)
            a=ad.AnnData(csr_matrix(counts),
                obs=pd.DataFrame({'type':['a','a','b','b']},index=['c1','c2','c3','c4']),
                var=pd.DataFrame(index=['L','R','unused']),
                obsm={'spatial':np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.]])})
            a.layers['raw_counts']=csr_matrix(counts*2)
            a.write_h5ad(root/'section.h5ad')
            sender=np.array([0,1,2,3]); receiver=np.array([2,3,0,1])
            np.savez(root/'graph.npz',unit_ids=np.array(['c1','c2','c3','c4']),
                     sender_index=sender,receiver_index=receiver)
            prepare_h5ad_section(root/'section.h5ad',['R','L'],root/'graph.npz',
                                 root/'out.npz',label_key='type')
            data=SpatialData.load(root/'out.npz')
            np.testing.assert_array_equal(data.counts,counts[:,[1,0]])
            np.testing.assert_array_equal(data.sender,sender)
            np.testing.assert_array_equal(data.receiver,receiver)
            with np.load(root/'out.npz',allow_pickle=False) as z:
                np.testing.assert_array_equal(z['unit_ids'],a.obs_names.to_numpy())
            h=prepare_hypothesis(data,'L','R','a','b')
            self.assertEqual(h['summary']['n_section'],4)
            self.assertEqual(h['summary']['n_receiver'],2)
            prepare_h5ad_section(root/'section.h5ad',['R','L'],root/'graph.npz',
                                 root/'layer.npz',label_key='type',layer='raw_counts')
            np.testing.assert_array_equal(SpatialData.load(root/'layer.npz').counts,counts[:,[1,0]]*2)
            np.savez(root/'bad_graph.npz',unit_ids=np.array(['c2','c1','c3','c4']),
                     sender=sender,receiver=receiver)
            with self.assertRaisesRegex(ValueError,'observation order'):
                prepare_h5ad_section(root/'section.h5ad',['L'],root/'bad_graph.npz',
                                     root/'bad.npz',label_key='type')
            with self.assertRaisesRegex(ValueError,'absent from H5AD'):
                prepare_h5ad_section(root/'section.h5ad',['missing'],root/'graph.npz',
                                     root/'missing.npz',label_key='type')
            with self.assertRaisesRegex(ValueError,'already exists'):
                prepare_h5ad_section(root/'section.h5ad',['L'],root/'graph.npz',
                                     root/'out.npz',label_key='type')
