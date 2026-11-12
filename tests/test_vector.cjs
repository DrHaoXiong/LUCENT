const assert=require('node:assert/strict');
const {SVGContext,escape}=require('../src/lucent/report/assets/vector.js');
const context=new SVGContext();
context.save();context.beginPath();context.rect(1,2,30,40);context.clip();
context.fillStyle='#0072b2';context.beginPath();context.arc(5,6,2,0,2*Math.PI);context.fill();
context.restore();context.fillText('A<&" μm',12,20);
const svg=context.markup();
assert(svg.includes('clip-path="url(#clip0)"'));
assert(svg.includes('A2 2 0 1 0'));
assert(svg.includes('A&lt;&amp;&quot; μm'));
assert(!svg.includes('<image'));
assert.equal(escape('</script>'), '&lt;/script&gt;');
assert.throws(()=>context.lineTo(Infinity,0),/Invalid vector coordinate/);
console.log('SVG recording tests PASS');
// Display simplification must preserve adjustments and the applicable family in downloads.
const {unitsTSV,familySize}=require('../src/lucent/report/assets/app.js');
const version={id:'main',adjustment_method:'BY',family_size:null,family_sizes:{S3:2848,S5:8722}};
const units=['S3','S5'].map(context=>({unit_id:context,channel_id:'channel',lr_pair:'L → R',context,
  sender_state:'A',receiver_state:'B',ligand_genes:['L'],receptor_genes:['R'],
  versions:{main:{score:.02,p:.005,q:.324,status:'complete',evaluable:true}}}));
const before=JSON.stringify(units),lines=unitsTSV(units,version).trimEnd().split('\r\n').map(line=>line.split('\t'));
for(let i=1;i<lines.length;i++){
  assert.equal(lines[i][lines[0].indexOf('q')],'0.324');
  assert.equal(lines[i][lines[0].indexOf('family_size')],i===1?'2848':'8722');
}
assert.equal(JSON.stringify(units),before);
assert.equal(familySize({...units[0],versions:{main:{family_size:123}}},version),123);
assert.equal(familySize(units[0],{id:'main',family_size:42}),42);
console.log('Report adjustment download tests PASS');
const {gwasProbabilityDisplay}=require('../src/lucent/report/assets/app.js');
assert.equal(gwasProbabilityDisplay({gwas_p:null,p_value_status:'source_reported_zero',neglog10p:'549.014'}), 'Source reported 0; −log₁₀(p) = 549.014');
assert.equal(gwasProbabilityDisplay({gwas_p:null,p_value_status:'below_float_range',p_value_raw:'1e-500'}), '1e-500 (source)');
assert.equal(gwasProbabilityDisplay({gwas_p:null,p_value_status:'rounded_at_threshold',p_value_raw:'4.9999999999999999e-8'}), '4.9999999999999999e-8 (source)');
assert.equal(gwasProbabilityDisplay({gwas_p:1e-9}), '1.00e-9');
assert.equal(gwasProbabilityDisplay({gwas_p:null}), 'Not provided');
console.log('Genetic source probability display tests PASS');
