// Byte identities only; no scientific computation or external project access.
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
const root=path.resolve(import.meta.dirname,'..');
const args=process.argv.slice(2);
if(args.length>1||(args.length===1&&args[0]!=='--with-local-inputs'))throw Error('usage: [--with-local-inputs]');
const sha=p=>crypto.createHash('sha256').update(fs.readFileSync(p)).digest('hex');
function check(r){
 const p=path.resolve(root,r.path);
 if(!p.startsWith(root+path.sep)||fs.lstatSync(p).isSymbolicLink()||fs.statSync(p).size!==r.bytes||sha(p)!==r.sha256)
  throw Error('package identity mismatch '+r.path);
}
const published=JSON.parse(fs.readFileSync(path.join(root,'provenance/package_files.json')));
const sums=fs.readFileSync(path.join(root,'SHA256SUMS'),'utf8').trimEnd().split('\n');
const indexLine=sums.find(l=>l.endsWith('  provenance/package_files.json'));
if(indexLine!==sha(path.join(root,'provenance/package_files.json'))+'  provenance/package_files.json')
 throw Error('package index checksum mismatch');
for(const r of published.files)check(r);
const expectedSums=[...published.files.map(r=>r.sha256+'  '+r.path),indexLine].sort();
if(JSON.stringify([...sums].sort())!==JSON.stringify(expectedSums))throw Error('checksum membership mismatch');
let present=0,missing=0;
const local=JSON.parse(fs.readFileSync(path.join(root,'provenance/local_inputs.json')));
for(const r of local.files){
 if(fs.existsSync(path.join(root,r.path))){check(r);present++;}else missing++;
}
if(args.length&&missing)throw Error('missing excluded local inputs: '+missing);
const expected=new Set([...published.files.map(r=>r.path),'provenance/package_files.json','SHA256SUMS',...local.files.map(r=>r.path)]);
function inspect(relative=''){
 for(const entry of fs.readdirSync(path.join(root,relative),{withFileTypes:true})){
  const name=path.join(relative,entry.name);
  if(name==='build'||name==='.git')continue;
  if(entry.isSymbolicLink())throw Error('unregistered symlink '+name);
  if(entry.isDirectory())inspect(name);
  else if(!expected.has(name))throw Error('unregistered package file '+name);
 }
}
inspect();
console.log(`PASS ${published.files.length} published file identities; ${present} local-only inputs checked, ${missing} absent. No numerical program executed.`);
