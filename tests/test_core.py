import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from lucent import SpatialData, CountCalibrationConfig, prepare_hypothesis, test_hypothesis
from lucent.spatial import directed_knn, positive_ranks
from lucent.calibration.refit import positive_ecdf, calibrate_count_reference, _rng
from lucent.calibration.copula import CountGeometry, GRID_FRACTIONS, GRID_RANGES
from lucent.calibration.geometry import LazyGeometryModel
from lucent.calibration.backend import CPUBackend
from lucent.provenance import core_mapping, nominate_channels


class CoreTests(unittest.TestCase):
    def fixture(self):
        rng=np.random.default_rng(7)
        xy=rng.uniform(size=(36,2))
        raw=rng.poisson(1.5,size=(36,2))
        labels=np.array(['A']*18+['B']*18)
        s,r=directed_knn(xy,4)
        return SpatialData(raw,np.array(['L','R']),xy,labels,s,r)

    def test_centering_and_complete_receiver(self):
        data=self.fixture()
        h=prepare_hypothesis(data,'L','R','A','B')
        self.assertEqual(len(h['raw']),18)
        self.assertEqual(h['outside'].sum(),18)
        self.assertAlmostEqual(h['contrast'].sum(),0,places=14)
        summary=h['summary']
        self.assertAlmostEqual(summary['centered_score'],summary['product_score']-summary['marginal_product'],places=14)
        y=positive_ecdf(h['raw'][:,None],h['values'],h['outside'])[:,0]
        np.testing.assert_array_equal(y,positive_ranks(data.counts[:,1])[18:])

    def test_no_positive_conditioning(self):
        data=self.fixture()
        data.counts[18:22,1]=0
        h=prepare_hypothesis(data,'L','R','A','B')
        self.assertEqual(len(h['raw']),18)
        self.assertTrue(np.all(h['raw'][:4]==0))

    def test_nonpositive_and_missing_are_distinct(self):
        data=self.fixture()
        data.counts[:,1]=0
        out=test_hypothesis(data,'L','R','A','B',seed=1,threads=1)
        self.assertFalse(out['calibration']['evaluable'])
        self.assertIsNone(out['calibration']['p'])

    def test_nonpositive_score_skips_simulation(self):
        class UnusedBackend:
            def fit(self, raw):
                raise AssertionError('nonpositive scores must not fit a reference')
        result=calibrate_count_reference(np.array([0,1,2]), np.array([1.,0.,-1.]),
                                         np.array([0,1,2]), np.zeros(3,dtype=int),
                                         UnusedBackend(),seed=1)
        self.assertTrue(result['evaluable'])
        self.assertEqual(result['p'],1.)
        self.assertEqual(result['outer_done'],0)

    def test_exact_inner_stopping_matches_full_draws(self):
        # Independent finite-budget loop: every inner draw is generated.
        data=self.fixture()
        h=prepare_hypothesis(data,'L','R','A','B')
        c=h['contrast'].copy()
        observed=positive_ecdf(h['raw'][:,None],h['values'],h['outside'])
        if float(np.einsum('i,ij->j',c,observed,optimize=False)[0]) < 0:
            c=-c
        geom=CountGeometry(h['coordinates'],LazyGeometryModel(h['coordinates'],GRID_FRACTIONS,GRID_RANGES))
        backend=CPUBackend(geom,1)
        cfg=CountCalibrationConfig(inner_draws=19,outer_draws=23,inner_batch=4)
        seed=901
        score=lambda y:np.einsum('i,ij->j',c,y,optimize=False)
        transform=lambda raw:positive_ecdf(raw,h['values'],h['outside'])
        try:
            actual=calibrate_count_reference(h['raw'],c,h['values'],h['outside'],backend,seed=seed,config=cfg)
            fit=backend.fit(h['raw'])
            draw=backend.bootstrap(fit,_rng(seed,'observed_inner').normal(size=(19,len(c))).T)
            k0=np.count_nonzero(score(transform(draw))>=score(observed)[0])
            hits=0
            for b in range(23):
                raw=backend.bootstrap(fit,_rng(seed,'outer',b).normal(size=(len(c),1)))
                y=transform(raw);t=score(y)[0]
                if np.ptp(y)<=1e-15 or t<=1e-12:
                    qualifies=k0==19
                else:
                    fitted=backend.fit(raw[:,0])
                    eps=_rng(seed,'inner',b).normal(size=(19,len(c))).T
                    count=np.count_nonzero(score(transform(backend.bootstrap(fitted,eps)))>=t)
                    qualifies=count<=k0
                hits+=int(qualifies)
            self.assertEqual(actual['outer_hits'],hits)
            self.assertEqual(actual['p'],(hits+1)/24)
            self.assertEqual(actual['outer_done'],23)
        finally:
            backend.close()

    def test_core_fallback_is_locus_level(self):
        genes=[dict(gene_id='G1',chrom='1',start=100,end=200),dict(gene_id='G2',chrom='1',start=400,end=500)]
        variants=[dict(variant_id='v1',locus_id='locusA',chrom='1',position=150,p_value=1e-9),
                  dict(variant_id='v2',locus_id='locusA',chrom='1',position=300,p_value=1e-9),
                  dict(variant_id='v3',locus_id='locusB',chrom='1',position=300,p_value=1e-9,pip=.9)]
        result=core_mapping(variants,genes,window_bp=20)
        self.assertEqual([(x['variant_id'],x['gene_id']) for x in result],[('v1','G1'),('v3','G1'),('v3','G2')])
        self.assertEqual([x['mapping_basis'] for x in result],['gene_body','nearest_gene','nearest_gene'])
        channels,ledger=nominate_channels(result,[dict(channel_id='c',ligand='G1',receptor='G2')])
        self.assertEqual(len(channels),1)
        self.assertEqual(len(ledger),3)
        self.assertEqual({x['endpoint_role'] for x in ledger},{'ligand','receptor'})

if __name__=='__main__':
    unittest.main()
