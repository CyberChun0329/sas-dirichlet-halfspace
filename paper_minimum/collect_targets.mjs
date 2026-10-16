// File identity and literal CSV assembly only. Scientific statistics are Fortran.
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath} from 'node:url';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const [packet,destination]=process.argv.slice(2);
if(!packet||!destination)throw Error('usage: TASKS.csv NEW_AGGREGATE_DIRECTORY');
function local(p){const q=path.resolve(root,p);if(!q.startsWith(root+path.sep))throw Error('outside experiment');return q;}
if(fs.existsSync(local('runs/.unified_run.lock')))throw Error('measurement lock exists; collect only after queue exit');
const lines=fs.readFileSync(local(packet),'utf8').trim().split(/\r?\n/);
if(lines.shift()!=='job,family,case_index,repeats,reference,rows,values,keys,outdir,prep,budget,shared')throw Error('task schema');
const jobs=lines.map(s=>{const f=s.split(',');if(f.length!==12)throw Error('task row');return {id:f[0],family:f[1],n:f[3],dir:f[8],budget:f[10]==='T',shared:f[11]==='T'};});
if(jobs.length!==114||new Set(jobs.map(j=>j.id)).size!==114)throw Error('exact114 required');
const dest=local(destination);if(!destination.startsWith('runs/')||fs.existsSync(dest))throw Error('strict-new runs destination');
const receipt=[];
function csv(p){
  const bytes=fs.readFileSync(local(p));receipt.push({path:p,bytes:bytes.length,sha256:crypto.createHash('sha256').update(bytes).digest('hex')});
  const a=bytes.toString('utf8').trimEnd().split(/\r?\n/);if(!a[0])throw Error('empty '+p);return {header:a[0],rows:a.slice(1)};
}
// Validate each real completion marker before creating an aggregate.
for(const j of jobs)if(!fs.readFileSync(local(j.dir+'/progress.log'),'utf8').trimEnd().endsWith('COMPLETE FINITE_LISTED_POOL'))throw Error('incomplete '+j.id);
fs.mkdirSync(dest);
function assemble(name,kind,selected,families,tagged=false){
  let header,body=[];
  for(const j of selected){
    const c=csv(j.dir+'/'+kind);if(header&&header!==c.header)throw Error('schema changed: '+kind);header=c.header;
    const k=header.split(',').indexOf('family');
    for(const row of c.rows){
      if(!row)continue;
      if(families&&(!families.includes(row.split(',')[k])))continue;
      body.push(tagged?j.id+','+row:row);
    }
  }
  if(!header)throw Error('empty member set: '+name);
  fs.writeFileSync(path.join(dest,name),(tagged?'target_job,':'')+header+'\n'+body.join('\n')+(body.length?'\n':''),{flag:'wx'});
}
for(const file of ['rows.csv','windows.csv','preparation.csv','coverage.csv','cached_correctness_rows.csv'])assemble(file,file,jobs,null,true);
const groups=[
  ['main100',jobs.filter(j=>j.family==='paper'&&j.n==='1'),['DIRECT_2D']],
  ['repeat10',jobs.filter(j=>j.n==='10'),['DIRECT_2D']],
  ['repeat100',jobs.filter(j=>j.n==='100'),['DIRECT_2D']],
  ['uniform_N1',jobs.filter(j=>j.family==='temperature_scale'&&j.n==='1'&&!j.id.includes('0096')),['DIRECT_2D']],
  ['deep_Q',jobs.filter(j=>j.id==='temperature_scale_0096_N1'),['DIRECT_2D']],
  ['deficit',jobs.filter(j=>j.family==='deficit_original'),['DIRECT_2D']],
  ['shared7',jobs.filter(j=>j.shared),['SHARED_1D','ALL_SPACE']],
  ['budget3',jobs.filter(j=>j.budget),['DIRECT_2D']]
];
for(const axis of ['cpu','wall']){
  assemble(axis+'_all_comparisons_tagged.csv',axis+'/comparisons.csv',jobs,null,true);
  for(const [name,members,families]of groups){
    assemble(axis+'_'+name+'_comparisons.csv',axis+'/comparisons.csv',members,families);
    // Only the three predeclared plots have a budget ladder claim.
    const selected=name==='budget3'?members:members.filter(j=>j.budget);
    if(selected.length)assemble(axis+'_'+name+'_budgets.csv',axis+'/budget_frontiers.csv',selected,families);
    else {const c=csv(members[0].dir+'/'+axis+'/budget_frontiers.csv');fs.writeFileSync(path.join(dest,axis+'_'+name+'_budgets.csv'),c.header+'\n',{flag:'wx'});}
  }
}
fs.writeFileSync(path.join(dest,'assembly_receipt.json'),JSON.stringify({observed_utc:new Date().toISOString(),packet,jobs,inputs:receipt,scope:'literal CSV assembly; no numerical statistics or smoke inputs; formal measurement tables use target outputs; cached_correctness_rows retain historical provenance costs, excluded from scientific performance summary inputs'},null,2)+'\n',{flag:'wx'});
console.log(JSON.stringify({status:'ASSEMBLED_NOT_SCIENTIFICALLY_SUMMARIZED',jobs:jobs.length,destination}));
