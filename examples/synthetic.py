"""Small synthetic input → genetic mapping → calibration → offline report example."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from lucent import SpatialData, CountCalibrationConfig, test_hypothesis
from lucent.spatial import directed_knn, positive_ranks
from lucent.provenance import core_mapping, nominate_channels
from lucent.report import export_report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('outputs/example'))
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    xy = rng.uniform(size=(48,2))
    labels = np.where(xy[:,0]<.5, 'Stroma', 'Epithelium')
    counts = rng.poisson(np.exp(1.2*xy[:,1,None])*np.array([1.,.8]), size=(48,2))
    genes = np.array(['DEMO_L','DEMO_R'])
    sender, receiver = directed_knn(xy)
    data = SpatialData(counts, genes, xy, labels, sender, receiver)
    np.savez_compressed(out/'section.npz', counts=counts, genes=genes, coordinates=xy,
                        labels=labels, sender=sender, receiver=receiver)
    config = dict(ligand='DEMO_L', receptor='DEMO_R', sender_label='Stroma',
                  receiver_label='Epithelium', seed=1234, threads=1, inner_draws=39, outer_draws=49)
    (out/'hypothesis.json').write_text(json.dumps(config,indent=2)+'\n')
    variants = [dict(variant_id='1:150:A:G', locus_id='demo_locus', chrom='1',position=150,p_value=1e-10,pip=.7)]
    intervals = [dict(gene_id='DEMO_L',gene_symbol='DEMO_L',chrom='1',start=100,end=200)]
    catalogue = [dict(channel_id='demo_channel',ligand='DEMO_L',receptor='DEMO_R',evidence='Direct')]
    for name, rows in [('variants',variants),('genes',intervals),('catalogue',catalogue)]:
        with (out/(name+'.tsv')).open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
            writer.writeheader();writer.writerows(rows)
    _, links = nominate_channels(core_mapping(variants, intervals),catalogue)
    result = test_hypothesis(data,'DEMO_L','DEMO_R','Stroma','Epithelium',seed=1234,
                            config=CountCalibrationConfig(inner_draws=39,outer_draws=49),threads=1)
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    calibrated=result['calibration']
    bundle = {
        'metadata': {'schema_version':'1.0','title':'LUCENT synthetic example','trait':'Synthetic trait',
                     'dataset':'Synthetic section','claim':'Explore genetic provenance and spatial coordination.',
                     'source_status':'Synthetic example','default_version':'lucent','default_section_id':'demo',
                     'versions':[{'id':'lucent','label':'LUCENT','status':'complete','p_label':'Calibrated association p',
                                  'q_label':'Adjusted p','adjustment_method':None,'family_size':1}],
                     'sources':[], 'notes':['Synthetic measurements and genetic identifiers; 39 inner and 49 outer draws.']},
        'units':[{'unit_id':'demo_unit','channel_id':'demo_channel','lr_pair':'DEMO_L → DEMO_R','context':'demo',
                  'sender_state':'Stroma','receiver_state':'Epithelium','ligand_genes':['DEMO_L'],
                  'receptor_genes':['DEMO_R'],'versions':{'lucent':{'score':calibrated['score'],'p':calibrated['p'],
                  'q':None,'status':'complete','evaluable':calibrated['evaluable']}}}],
        'genetic_links':[{**x,'gwas_p':x['p_value'],'pip_status':'Synthetic annotation','build':'synthetic'} for x in links],
        'donor_evidence':[],
        'sections':[{'id':'demo','donor_id':'demo_donor','context':'demo','label':'Synthetic section',
                     'coordinate_system':{'unit':'arbitrary','y_direction':'up'},
                     'spots':[{'id':str(i),'state':str(labels[i]),'x':float(x),'y':float(y)} for i,(x,y) in enumerate(xy)],
                     'edges':np.column_stack((sender,receiver)).tolist(),
                     'expression':{g:positive_ranks(counts[:,j]).tolist() for j,g in enumerate(genes)},
                     'expression_label':'Positive expression ranks'}]}
    (out/'report_bundle.json').write_text(json.dumps(bundle,indent=2)+'\n')
    export_report(bundle,out/'report')
    print(json.dumps({'result':str(out/'result.json'),'report':str(out/'report/index.html'),
                      'score':calibrated['score'],'p':calibrated['p']}))

if __name__=='__main__':
    main()
