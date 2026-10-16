#!/usr/bin/env node
// Filesystem/process orchestration only. All numerical work stays in Fortran.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { spawnSync } from 'node:child_process';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const usage = `Usage (compile first with make -j2):
  node tools/run_paper.mjs validation core|case1|gaussian|table_s2 --out build/NEW_RUN
  node tools/run_paper.mjs validation roots CHUNK_1_TO_8 --out build/NEW_RUN
  node tools/run_paper.mjs precision --out build/NEW_RUN
  node tools/run_paper.mjs periodic unit --out build/NEW_RUN
  node tools/run_paper.mjs periodic grid INDEX_1_TO_8 --out build/NEW_RUN
  node tools/run_paper.mjs periodic invariants INDEX_1_TO_4 --out build/NEW_RUN
  node tools/run_paper.mjs check-timing JOB_ID --out build/NEW_JOB
  node tools/run_paper.mjs timing JOB_ID --out build/NEW_JOB
No arguments or --help only displays this message. Each output must be new.
Validation commands execute one bank/chunk; they do not aggregate a full suite.
Timing commands select exactly one row of paper_minimum/tasks.csv, retain its
saved input paths, and change only outdir in a new build-local task file.`;

function fail(message) { throw new Error(message); }
function read(relative) { return fs.readFileSync(path.join(root, relative), 'utf8'); }
function run(binary, args, cwd = root) {
  const executable = path.join(root, 'build', 'release', binary);
  if (!fs.existsSync(executable)) fail(`Compile first: missing build/release/${binary}`);
  const result = spawnSync(executable, args, {
    cwd, stdio: 'inherit', shell: false,
    env: { ...process.env, VECLIB_MAXIMUM_THREADS: '1', OPENBLAS_NUM_THREADS: '1', OMP_NUM_THREADS: '1' },
  });
  if (result.error) throw result.error;
  if (result.signal) fail(`${binary} terminated by ${result.signal}; partial output retained`);
  if (result.status !== 0) fail(`${binary} exited ${result.status}; partial output retained`);
}
function outputPath(value) {
  if (!/^build\/[A-Za-z0-9_-]+(?:\/[A-Za-z0-9_-]+)*$/.test(value ?? ''))
    fail('--out must be a new root-relative path below build, using letters, digits, _ or -');
  if (value.length > 900) fail('Output path too long for native path fields');
  const absolute = path.join(root, value);
  // Reject symlink ancestors before mkdir so build cannot redirect evidence writes.
  let current = root;
  for (const component of value.split('/')) {
    current = path.join(current, component);
    if (fs.existsSync(current) && fs.lstatSync(current).isSymbolicLink()) fail('Output ancestors must not be symlinks');
  }
  if (fs.existsSync(absolute)) fail(`Output already exists: ${value}`);
  return absolute;
}

try {
  const args = process.argv.slice(2);
  if (args.length === 0 || (args.length === 1 && args[0] === '--help')) {
    console.log(usage);
  } else {
    const outIndex = args.indexOf('--out');
    if (outIndex < 1 || outIndex !== args.length - 2) fail('Supply one final --out build/NEW_RUN');
    const outRelative = args[outIndex + 1];
    const outAbsolute = outputPath(outRelative);
    const [command, ...operands] = args.slice(0, outIndex);
    let binary, nativeArgs;
    if (command === 'timing' || command === 'check-timing') {
      if (operands.length !== 1 || !/^[A-Za-z0-9_-]+$/.test(operands[0])) fail('Exactly one explicit JOB_ID required');
      const lines = read('paper_minimum/tasks.csv').trimEnd().split(/\r?\n/);
      const header = 'job,family,case_index,repeats,reference,rows,values,keys,outdir,prep,budget,shared';
      if (lines[0] !== header) fail('Native task header mismatch');
      const matches = lines.slice(1).filter(line => line.split(',')[0] === operands[0]);
      if (matches.length !== 1) fail('JOB_ID must match exactly one saved paper task');
      const fields = matches[0].split(',');
      if (fields.length !== 12) fail('Native task field count mismatch');
      for (const index of [4, 5, 6, 7]) {
        const relative = fields[index];
        if (!/^[A-Za-z0-9_./-]+$/.test(relative) || path.isAbsolute(relative) || relative.split('/').includes('..'))
          fail('Unsafe saved input path');
        if (!fs.statSync(path.join(root, relative)).isFile()) fail(`Missing saved task input: ${relative}`);
      }
      const release = read('config/source_release.txt').trim();
      if (!/^[a-f0-9]{64}$/.test(release)) fail('Expected frozen 64-character source release identity');
      fields[8] = outRelative;
      const taskRelative = `${outRelative}.task.csv`;
      fs.mkdirSync(path.dirname(outAbsolute), { recursive: true });
      fs.writeFileSync(path.join(root, taskRelative), `${header}\n${fields.join(',')}\n`, { flag: 'wx' });
      run('targeted_driver', [command === 'timing' ? '--run' : '--check', taskRelative,
        operands[0], 'paper_minimum/policy.nml', release]);
    } else {
      if (command === 'validation') {
        const [bank, chunk] = operands;
        if (!['core', 'case1', 'gaussian', 'table_s2', 'roots'].includes(bank)) fail('Choose a named validation bank');
        if (bank === 'roots') {
          if (operands.length !== 2 || !/^[1-8]$/.test(chunk)) fail('roots requires one chunk 1..8');
        } else if (operands.length !== 1) fail('This validation bank takes no chunk');
        binary = 'validation_driver'; nativeArgs = operands;
      } else if (command === 'precision' && operands.length === 0) {
        binary = 'precision_driver'; nativeArgs = ['results'];
      } else if (command === 'periodic') {
        const [mode, index] = operands;
        if (mode === 'unit' && operands.length === 1) nativeArgs = operands;
        else if (mode === 'grid' && operands.length === 2 && /^[1-8]$/.test(index)) nativeArgs = operands;
        else if (mode === 'invariants' && operands.length === 2 && /^[1-4]$/.test(index)) nativeArgs = operands;
        else fail('Choose periodic unit, grid 1..8, or invariants 1..4');
        binary = 'periodic_driver';
      } else fail('Unknown command or operands; use --help');
      fs.mkdirSync(path.join(outAbsolute, 'results'), { recursive: true });
      run(binary, nativeArgs, outAbsolute);
    }
  }
} catch (error) {
  console.error(error.message);
  process.exitCode = 1;
}
