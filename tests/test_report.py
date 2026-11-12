from pathlib import Path
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
