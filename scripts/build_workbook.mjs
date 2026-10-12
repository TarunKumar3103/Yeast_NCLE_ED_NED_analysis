import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const root = process.argv[2];
if (!root) throw new Error('Pass the absolute package directory.');
const p = JSON.parse(await fs.readFile(path.join(root, 'results/workbook_payload.json'), 'utf8'));
const wb = Workbook.create();
const results = wb.worksheets.add('Results');
const models = wb.worksheets.add('Model coefficients');
const flow = wb.worksheets.add('Cohort flow');
const records = wb.worksheets.add('Protein records');
const navy = '#243D59', ink = '#172B40', pale = '#EAF0F6';
for (const s of [results, models, flow, records]) s.showGridLines = false;
results.tabColor = navy;
function style(s, range) {
  const f = s.getRange(range).format;
  f.font = {name:'Arial',size:11,color:ink};
  f.verticalAlignment = 'center';
  f.rowHeight = 23;
}
function title(s, text, last='I') {
  s.getRange('A2').values = [[text]];
  s.getRange('A2').format.font = {name:'Arial',size:16,bold:true,color:ink};
  s.getRange('A2').format.rowHeight = 30;
  s.getRange(`A3:${last}3`).format.borders = {bottom:{style:'thin',color:'#A9B9CA'}};
}
function header(s, range) {
  const f = s.getRange(range).format;
  f.fill = navy;
  f.font = {name:'Arial',size:11,bold:true,color:'#FFFFFF'};
  f.horizontalAlignment = 'center';
  f.wrapText = true;
  f.rowHeight = 38;
  f.borders = {insideVertical:{style:'thin',color:'#FFFFFF'}};
}
function section(s, cell, text) {
  s.getRange(cell).values = [[text]];
  s.getRange(cell).format.font = {name:'Arial',size:12,bold:true,color:ink};
}
function col(i) {
  let text='';
  for(let n=i+1;n>0;n=Math.floor((n-1)/26)) text=String.fromCharCode(65+(n-1)%26)+text;
  return text;
}

const keys=Object.keys(p.audit[0]);
const last=col(keys.length-1);
const matrix=[keys,...p.audit.map(r=>keys.map(k=>r[k]??null))];
const end=matrix.length;
style(records,`A1:${last}${end}`);
records.getRangeByIndexes(0,0,end,keys.length).values=matrix;
records.tables.add(`A1:${last}${end}`,true,'YeastProteinRecords');
header(records,`A1:${last}1`);
records.getRange(`A1:${last}1`).format.rowHeight=55;
records.freezePanes.freezeRows(1);
records.freezePanes.freezeColumns(1);
for(let i=0;i<keys.length;i++) {
  const c=col(i), k=keys[i];
  records.getRange(`${c}1:${c}${end}`).format.columnWidth=
    k==='UniProt_ID'?20:k==='Disposition'?65:k==='Protein name'?72:k==='Gene name'?24:16;
  if(k==='R2') records.getRange(`${c}2:${c}${end}`).setNumberFormat('0.000000');
  else if(k==='Mean_pLDDT') records.getRange(`${c}2:${c}${end}`).setNumberFormat('0.00');
  else if(!['UniProt_ID','Disposition','Protein name','Gene name','Degradation_profile'].includes(k))
    records.getRange(`${c}2:${c}${end}`).setNumberFormat('0');
}
const included=col(keys.indexOf('Included'));
const ncle=col(keys.indexOf('NCLE'));
const r2=col(keys.indexOf('R2'));
const refs={included:`'Protein records'!$${included}$2:$${included}$${end}`,
            ncle:`'Protein records'!$${ncle}$2:$${ncle}$${end}`,
            r2:`'Protein records'!$${r2}$2:$${r2}$${end}`};
records.getRange(`${included}2:${included}${end}`).conditionalFormats.add('cellIs',{
  operator:'equal',formula:0,format:{fill:'#FBE8E5',font:{color:'#A1362B'}}});
const src=[
  ['Source','Fields','Frozen file or definition','Method'],
  ['Christiano 2014 supplement','Accession, R², gene/protein names','Christiano_yeast_mmc2.xlsx, Sheet1','https://doi.org/10.1016/j.celrep.2014.10.065'],
  ['Supplied native annotations','NCLE, mean pLDDT, protein length','AF_Yeast_entanglement.pkl','entanglement_AF Yes=1, No=0; supplied plddt and length_AF'],
  ['Supplied topology exclusions','Knot and covalent-lasso flags','knot_all.txt; AlphaLasso.zip','Exact accession match; normalize terminal -F<number> only'],
  ['Frozen UniProt snapshot','Membrane and Secreted keywords','taxonomy_id:4932; retrieved 2026-10-10','KW-0472 and KW-0964; metadata retained in raw_sources/'],
  ['Quyen instruction','Primary ED/NED classification','R² >=0.90: ED; R² <0.90: NED','Operational single-exponential-fit classification'],
  ['Cutoff sensitivity','ED/NED thresholds','0.75, 0.85 and 0.95','Same cohort; every result reported'],
  ['Length adjustment','z(protein length)','Eligible cohort mean and sample SD','NED ~ NCLE + z(length); length_AF preserved'],
  ['Interpretation','R²-defined NED','Poor fit to a single-exponential model','Does not establish bi-exponential decay or misfolding'],
];
style(records,'Z1:AC9');
records.getRange('Z1:AC9').values=src;
header(records,'Z1:AC1');
for(const [c,w] of [['Z',33],['AA',38],['AB',54],['AC',74]]) records.getRange(`${c}1:${c}9`).format.columnWidth=w;

style(results,'A1:I38');
title(results,'Yeast NCLE–NED association');
const primaryCutoff=p.manifest.Primary_R2_cutoff.toFixed(2);
results.getRange('A4').values=[[`Primary: R² <${primaryCutoff} = NED; R² ≥${primaryCutoff} = ED. Eligible cohort: ${p.manifest.Final_N.toLocaleString('en-US')} proteins.`]];
const widths=[24,12,23,10,13,14,16,16,13];
widths.forEach((w,i)=>results.getRange(`${col(i)}1:${col(i)}38`).format.columnWidth=w);
results.getRange('A6:I6').values=[['Analysis','R² cutoff','Method','N','NED count','NCLE OR','95% CI lower','95% CI upper','P-value']];
header(results,'A6:I6');
const estimates=p.results.filter(r=>r.Term==='NCLE');
results.getRange('A7:I14').values=estimates.map(r=>[
  r.Analysis,r.R2_cutoff,r.Method==='Model 1'?'Length-adjusted':'Crude',r.N,
  p.manifest.Model_diagnostics.find(x=>x.R2_cutoff===r.R2_cutoff).Events,
  r.OR,r.CI_lower,r.CI_upper,r.P_value]);
results.getRange('A7:I8').format.fill=pale;
results.getRange('B7:B14').setNumberFormat('0.00');
results.getRange('D7:E14').setNumberFormat('#,##0');
results.getRange('F7:I14').setNumberFormat('0.000');
section(results,'A17','Contingency tables');
results.getRange('A19:G19').values=[['Analysis','R² cutoff','NCLE status','NED','ED','Total','NED frequency']];
header(results,'A19:G19');
results.getRange('A20:G27').values=p.counts.map(r=>[
  r.Analysis,r.R2_cutoff,r.NCLE===1?'Entangled':'Not entangled',r.NED,r.ED,r.Total,r.NED_fraction]);
for(let i=0;i<p.counts.length;i++) {
  const row=20+i, n=p.counts[i].NCLE;
  const base=`${refs.included},1,${refs.ncle},${n},${refs.r2}`;
  results.getRange(`D${row}`).formulas=[[`=COUNTIFS(${base},"<"&B${row})`]];
  results.getRange(`E${row}`).formulas=[[`=COUNTIFS(${base},">="&B${row})`]];
  results.getRange(`F${row}`).formulas=[[`=SUM(D${row}:E${row})`]];
  results.getRange(`G${row}`).formulas=[[`=D${row}/F${row}`]];
}
results.getRange('A20:G21').format.fill=pale;
results.getRange('B20:B27').setNumberFormat('0.00');
results.getRange('D20:F27').setNumberFormat('#,##0');
results.getRange('G20:G27').setNumberFormat('0.00%');
for(let i=0;i<4;i++) {
  const out=7+2*i, a=20+2*i, b=a+1;
  results.getRange(`D${out}`).formulas=[[`=SUM(F${a}:F${b})`]];
  results.getRange(`E${out}`).formulas=[[`=SUM(D${a}:D${b})`]];
  results.getRange(`F${out}`).formulas=[[`=D${a}*E${b}/(E${a}*D${b})`]];
  const se=`SQRT(1/D${a}+1/E${a}+1/D${b}+1/E${b})`;
  results.getRange(`G${out}`).formulas=[[`=EXP(LN(F${out})-1.95996398454005*${se})`]];
  results.getRange(`H${out}`).formulas=[[`=EXP(LN(F${out})+1.95996398454005*${se})`]];
}
const notes=[
  'Model 1: NED ~ NCLE + z(protein length). OR compares entangled with non-entangled proteins.',
  estimates.every(r=>r.P_value>=.05 && r.CI_lower<=1 && r.CI_upper>=1)
    ? 'No significant NCLE association at any tested cutoff; every NCLE 95% CI includes 1.'
    : 'At least one NCLE estimate reaches P <0.05 or has a 95% CI excluding 1; see all cutoff results.',
  'NED is defined by R² here. Poor fit does not establish bi-exponential kinetics or misfolding.',
  'At cutoff 0.75 there are only nine NED proteins; the regression estimate is imprecise.',
  'Crude P: two-sided Fisher exact. Adjusted P and CI: logistic Wald. No multiplicity adjustment.',
  'Counts and crude OR/CI use Excel formulas. Python regression snapshots require a rerun after edits.',
];
notes.forEach((text,i)=>results.getRange(`A${30+i}`).values=[[text]]);

style(models,'A1:J36');
title(models,'Logistic regression coefficients');
models.getRange('A4').values=[['Model 1: NED ~ NCLE + z(protein length), fitted separately at each R² cutoff.']];
const terms=p.results.filter(r=>r.Method==='Model 1');
models.getRange('A6:I6').values=[['Analysis','R² cutoff','Term','N','Coefficient','Standard error','OR','95% CI lower','95% CI upper']];
header(models,'A6:I6');
models.getRange('A7:I18').values=terms.map(r=>[r.Analysis,r.R2_cutoff,r.Term,r.N,r.Coefficient,r.SE,r.OR,r.CI_lower,r.CI_upper]);
models.getRange('J6').values=[['P-value']];
header(models,'J6:J6');
models.getRange('J7:J18').values=terms.map(r=>[r.P_value]);
const mw=[24,12,23,10,17,17,14,16,16,17];
mw.forEach((w,i)=>models.getRange(`${col(i)}1:${col(i)}36`).format.columnWidth=w);
models.getRange('A7:J9').format.fill=pale;
models.getRange('B7:B18').setNumberFormat('0.00');
models.getRange('D7:D18').setNumberFormat('#,##0');
models.getRange('E7:I18').setNumberFormat('0.000000');
models.getRange('J7:J18').setNumberFormat('0.000E+00');
section(models,'A21','Fit diagnostics');
models.getRange('A23:H23').values=[['R² cutoff','NED count','ED count','Smallest 2×2 cell','Status','Max absolute score','Information condition','Iterations']];
header(models,'A23:H23');
models.getRange('A24:H27').values=p.manifest.Model_diagnostics.map(r=>[
  r.R2_cutoff,r.Events,r.Non_events,r.Smallest_2x2_cell,r.Status,r.Max_abs_score,r.Information_condition_number,r.Iterations]);
models.getRange('A24:A27').setNumberFormat('0.00');
models.getRange('B24:D27').setNumberFormat('#,##0');
models.getRange('F24:F27').setNumberFormat('0.00E+00');
models.getRange('G24:G27').setNumberFormat('0.00');
const diag=p.manifest.Model_diagnostics[0];
models.getRange('A30').values=[[`z(length) = (length − ${diag.Mean_length.toFixed(6)} residues) / ${diag.SD_length_ddof1.toFixed(6)} residues (sample SD, ddof=1).`]];
models.getRange('A31').values=[['Length OR is per one SD increase. Intercept OR is baseline odds, not an exposure effect.']];
models.getRange('A32').values=[['Convergence, design rank and complete/quasi-separation were checked in Python.']];

style(flow,'A1:J21');
title(flow,'Yeast cohort construction','J');
flow.getRange('A4').values=[['Same eligible cohort for the primary analysis and all cutoff-sensitivity analyses.']];
const cutoffs=p.manifest.R2_cutoffs;
flow.getRange('A6:J6').values=[['Sequential step','Remaining','Removed','NCLE=1','NCLE=0','NCLE missing',...cutoffs.map(c=>`NED <${c.toFixed(2)}`)]];
header(flow,'A6:J6');
flow.getRange('A7:J15').values=p.flow.map(r=>[r.Step,r.Remaining,r.Removed_at_step,r.NCLE_1,r.NCLE_0,r.NCLE_missing,...cutoffs.map(c=>r[`NED_${String(c).replace('.','_')}`])]);
flow.getRange('A1:A21').format.columnWidth=78;
flow.getRange('B1:J21').format.columnWidth=15;
flow.getRange('B7:J15').setNumberFormat('#,##0');
flow.getRange('A15:J15').format.fill=pale;
flow.getRange('A17').values=[['Grouped accession records are not expanded into independent proteins. Exact accession matches are used.']];
flow.getRange('A18').values=[['No disorder-fraction filter or complex-membership covariate is used in these statistical analyses.']];
flow.getRange('A19').values=[['Supplied AF length is retained for P40516 (257 residues); current UniProt length is 229 residues.']];

wb.recalculate();
const calculated=results.getRange('D20:G27').values;
for(let i=0;i<p.counts.length;i++) {
  const e=p.counts[i],v=calculated[i];
  if(v[0]!==e.NED || v[1]!==e.ED || v[2]!==e.Total || Math.abs(v[3]-e.NED_fraction)>1e-12)
    throw new Error(`Contingency formula disagreement at row ${20+i}`);
}
for(let i=0;i<4;i++) {
  const v=results.getRange(`F${7+2*i}:H${7+2*i}`).values[0],e=estimates[2*i];
  if(Math.abs(v[0]-e.OR)>1e-10 || Math.abs(v[1]-e.CI_lower)>1e-10 || Math.abs(v[2]-e.CI_upper)>1e-10)
    throw new Error('Crude OR/CI formulas disagree with Python');
}
const err=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:1800});
console.log(err.ndjson);
const inspect=await wb.inspect({kind:'table',range:'Results!A6:I14',include:'values,formulas',tableMaxRows:9,tableMaxCols:9,maxChars:2500});
await fs.mkdir(path.join(root,'verification'),{recursive:true});
await fs.writeFile(path.join(root,'verification/workbook_inspection.ndjson'),inspect.ndjson);
for(const [sheet,range,file] of [
  ['Results','A1:I36','results_preview.png'],
  ['Model coefficients','A1:J33','models_preview.png'],
  ['Cohort flow','A1:J20','flow_preview.png'],
  ['Protein records','A1:K12','records_preview.png'],
  ['Protein records','Z1:AC9','sources_preview.png']]) {
  const blob=await wb.render({sheetName:sheet,range,scale:1.2,format:'png'});
  await fs.writeFile(path.join(root,'verification',file),new Uint8Array(await blob.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(root,'Yeast_NCLE_ED_NED_R2_results.xlsx'));
console.log('Workbook exported; all contingency and crude OR/CI formulas match Python.');
