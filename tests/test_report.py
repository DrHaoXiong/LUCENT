from pathlib import Path
import json
import tempfile
import unittest
from lucent.report import export_report
from lucent.report.verify import verify_report


class ReportTests(unittest.TestCase):
    def test_offline_vector_assets_and_manifest(self):
        bundle={'metadata':{'schema_version':'1.0','title':'Example','trait':'Example',
                'dataset':'Example','claim':'Spatial evidence','source_status':'Example',
                'default_version':'main','default_section_id':None,
                'versions':[dict(id='main',label='LUCENT',status='complete',p_label='p',
                                 q_label='Adjusted p',adjustment_method=None,family_size=0)],
                'sources':[],'notes':[]},'units':[],'sections':[],
                'genetic_links':[],'donor_evidence':[]}
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)/'report'
            manifest=export_report(bundle,root)
            html=(root/'index.html').read_text(encoding='utf-8')
            self.assertIn('id="export-pdf"',html)
            self.assertIn('id="export-svg"',html)
            self.assertLess(html.index('assets/vector.js'),html.index('assets/app.js'))
            self.assertTrue((root/'assets/vector.js').is_file())
            self.assertIn('assets/vector.js',[x['path'] for x in manifest['outputs']])
            self.assertFalse(verify_report(root)['statistics_recomputed'])

    def test_underflow_genetic_provenance_survives_export(self):
        bundle={'metadata':{'schema_version':'1.0','title':'Example','trait':'Example',
                'dataset':'Example','claim':'Spatial evidence','source_status':'Example',
                'default_version':'main','default_section_id':None,
                'versions':[dict(id='main',label='LUCENT',status='complete',p_label='p',
                                 q_label='Adjusted p',adjustment_method=None,family_size=1)],
                'sources':[],'notes':[]},
                'units':[dict(unit_id='u',channel_id='c',lr_pair='L → R',context='tissue',
                              sender_state='A',receiver_state='B',ligand_genes=['L'],receptor_genes=['R'],
                              versions={'main':dict(score=.1,p=.01,q=None,status='complete',evaluable=True)})],
                'sections':[], 'donor_evidence':[],
                'genetic_links':[dict(channel_id='c',locus_id='l',variant_id='v',gene_id='L',
                                      pip=None,pip_status='Not provided',gwas_p=None,
                                      p_value_raw='0.0',p_value_status='source_reported_zero',
                                      neglog10p='549.014')]}
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)/'report'
            export_report(bundle,root)
            link=json.loads((root/'downloads'/'bundle.json').read_text(encoding='utf-8'))['genetic_links'][0]
            self.assertIsNone(link['gwas_p'])
            self.assertEqual(link['p_value_raw'],'0.0')
            self.assertEqual(link['neglog10p'],'549.014')
            self.assertIn('source_reported_zero',(root/'downloads'/'genetic_links.tsv').read_text(encoding='utf-8'))
            self.assertFalse(verify_report(root)['statistics_recomputed'])
