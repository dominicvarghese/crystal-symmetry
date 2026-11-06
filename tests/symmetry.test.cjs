const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const html = fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8');
const scripts = [...html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/g)].map(m => m[1]);
for (const script of scripts) new vm.Script(script); // compile every inline script
const context = vm.createContext({});
vm.runInContext(scripts.find(s => s.includes('function analyzeSingleRep')), context);
const C = context.Crystallography;
const presets = context.PRESET_CRYSTALS;
const I = [[1,0,0],[0,1,0],[0,0,1]];
const close = (a, b, message = '') => assert.ok(C.vecNorm(C.vecSub(a, b)) < 2e-5, `${message}: ${a} != ${b}`);
const expected = {
  zns_216: [8, 96, '-43m'], silicon_227: [8, 192, 'm-3m'],
  rutile_136: [6, 16, '4/mmm'], chalcopyrite_122: [16, 16, '-42m'],
  quartz_152: [9, 6, '32'], baddeleyite_14: [12, 4, '2/m'], rubrene_64: [280, 16, 'mmm'],
};
for (const [key, [nAtoms, nOps, pg]] of Object.entries(expected)) {
  test(`${key}: operation count, group closure, species mapping and animation endpoints`, () => {
    const d = presets[key];
    assert.equal(d.atoms.length, nAtoms);
    assert.equal(d.operations.length, nOps);
    assert.equal(C.identifyCrystalSystemAndProperties(d.operations).pointGroup.hm, pg);
    const opKey = (W, w) => W.flat().join(',')+'|'+w.map(x => ((Math.round(x*1e5)%1e5)+1e5)%1e5).join(',');
    const keys = new Set(d.operations.map(o => opKey(o.W, o.w_raw)));
    for (const op of d.operations) {
      assert.equal(new Set(op.atom_mapping).size, nAtoms);
      for (const a of d.atoms) {
        const b = d.atoms[op.atom_mapping[a.id]];
        assert.equal(a.element, b.element);
        const delta = C.vecSub(C.vecAdd(C.matVec(op.W, a.frac), op.w_raw), b.frac).map(x=>x-Math.round(x));
        assert.ok(C.vecNorm(C.matVec(d.cell_matrix_cols, delta)) <= .02);
      }
      for (const a of d.packed_atoms) for (const style of ['two-stage', 'simultaneous']) {
        close(C.evaluateAnimationStep(a.cart, op, 0, style).posCurrent, a.cart, key+' start');
        const step = C.evaluateAnimationStep(a.cart, op, 1, style);
        const target = C.vecAdd(C.matVec(op.R, a.cart), op.t_cart);
        close(step.posCurrent, target, key+' end'); close(step.posFinal, target);
      }
      for (const other of d.operations) {
        assert.ok(keys.has(opKey(C.matMul(op.W, other.W), C.vecAdd(C.matVec(op.W, other.w_raw), op.w_raw))), key+' closure');
      }
    }
  });
}
test('ZnS offset 2_1: half turn reaches an edge midpoint, then translates along z', () => {
  const d = presets.zns_216;
  const op = d.operations.find(o => o.type === 'screw' && Math.abs(o.point_frac[0]-.5)<1e-6 && Math.abs(o.point_frac[1]-.25)<1e-6 && o.axis_frac[2]===1);
  assert.ok(op);
  const inv = C.matInv(d.cell_matrix_cols);
  close(C.matVec(inv, C.evaluateAnimationStep([0,0,0], op, .5).posCurrent), [1,.5,0]);
  close(C.matVec(inv, C.evaluateAnimationStep([0,0,0], op, 1).posCurrent), [1,.5,.5]);
  const rotated = C.evaluateAnimationStep([0,0,0], op, .5).posCurrent;
  const q = C.vecSub([0,0,0], op.point_cart), r = C.vecSub(rotated, op.point_cart);
  close(C.vecSub(q, C.vecScale(op.axis_cart, C.vecDot(q, op.axis_cart))), C.vecScale(C.vecSub(r, C.vecScale(op.axis_cart, C.vecDot(r, op.axis_cart))), -1));
});
test('Quartz Si has twofold site symmetry and is nonpolar', () => {
  const d = presets.quartz_152;
  assert.equal(C.getAtomSiteSymmetry(d, d.atoms[0]).pointGroup.hm, '2');
  assert.equal(C.identifyCrystalSystemAndProperties(d.operations).polar, false);
});
test('signed screw translation distinguishes 3_1/3_2 and 4_1/4_3', () => {
  const A = [[1,-.5,0],[0,Math.sqrt(3)/2,0],[0,0,1]], W = [[0,-1,0],[1,-1,0],[0,0,1]];
  assert.equal(C.classifySymmetryOperation(1,W,[0,0,1/3],A).symbol,'3_1');
  assert.equal(C.classifySymmetryOperation(1,W,[0,0,2/3],A).symbol,'3_2');
  const W4=[[0,-1,0],[1,0,0],[0,0,1]];
  assert.equal(C.classifySymmetryOperation(1,W4,[0,0,1/4],I).symbol,'4_1');
  assert.equal(C.classifySymmetryOperation(1,W4,[0,0,3/4],I).symbol,'4_3');
});
test('the requested tolerance is respected for every atom, including hydrogen', () => {
  const A = [[4,0,0],[0,5,0],[0,0,6]];
  const atoms = [{element:'C',frac:[0,0,0]}, {element:'C',frac:[.5,.5,.5]}, {element:'H',frac:[.12,.23,.31]}, {element:'H',frac:[.63,.74,.83]}];
  const ops = C.detectSymmetryOperations(A, atoms, .005);
  assert.equal(ops.length,1); // old detector silently used the carbon skeleton
  for(const op of ops) assert.ok(C.validatedAtomMapping(A, atoms, op.W, op.w_raw, .005));
});
test('matching is bijective and rejects a many-to-one atom mapping', () => {
  const atoms = [{element:'C',frac:[.1,0,0]},{element:'C',frac:[.2,0,0]},{element:'C',frac:[.6,0,0]}];
  assert.equal(C.validatedAtomMapping(I,atoms,I,[.5,0,0],.11),null);
});
test('nearby distinct sites are not merged by a fixed 0.06 Å cutoff', () => {
  const d = C.buildCrystalDataset({cellRows:[[10,0,0],[0,11,0],[0,0,12]], rawAtoms:[{element:'C',frac:[.1,.2,.3]},{element:'C',frac:[.104,.2,.3]}],symprec:.005});
  assert.equal(d.atoms.length,2);
});
const cell = `_cell_length_a '4.0(1)'\n_cell_length_b 5\n_cell_length_c 6\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n`;
const atomLoop = 'loop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n';
test('CIF loop rows may span lines or share a line, with inline comments', () => {
  const d=C.parseCIF(cell+atomLoop+'C .1\n.2 .3 # coordinate continuation\nO .4 .5 .6 H .7 .8 .9\n');
  assert.equal(d.atoms.length,3);close(d.atoms[0].frac,[.1,.2,.3]);
});
test('CIF rejects missing cell values, incomplete rows and invalid operators', () => {
  assert.throws(()=>C.parseCIF(atomLoop+'C 0 0 0'),/cell lengths/);
  assert.throws(()=>C.parseCIF(cell+atomLoop+'C 0 0'),/Incomplete/);
  assert.throws(()=>C.parseCIF(cell+'loop_\nC 0 0 0'),/column headers/);
  assert.throws(()=>C.parseCIF(cell+'loop_\n_space_group_symop_operation_xyz\n\'x+oops,y,z\'\n'+atomLoop+'C 0 0 0'),/Invalid CIF symmetry/);
});
test('CIF explicit operations must preserve both cell metric and atom positions', () => {
  assert.throws(()=>C.buildCrystalDataset({cellRows:[[4,0,0],[0,5,0],[0,0,6]],rawAtoms:[{element:'C',frac:[.1,.2,.3]}],explicitSymOps:[{W:I,w:[.5,0,0]}]}),/does not preserve/);
  assert.throws(()=>C.buildCrystalDataset({cellRows:[[4,0,0],[0,5,0],[0,0,6]],rawAtoms:[{element:'C',frac:[0,0,0]}],explicitSymOps:[{W:[[0,-1,0],[1,0,0],[0,0,1]],w:[0,0,0]}]}),/does not preserve/);
});
test('small nonzero strain breaks incompatible operations; animations stay rigid', () => {
  const d=C.applyOrderParameterAndStrain(presets.zns_216,{exx:.0001});
  assert.equal(d.preserved_operations.length,32);
  for(const op of d.operations) {
    for(const a of d.atoms) close(C.evaluateAnimationStep(a.cart,op,1).posCurrent,C.vecAdd(C.matVec(op.R,a.cart),op.t_cart));
    const a=[.3,.2,.8],b=[1.2,.5,.4];
    const delta=C.vecSub(C.evaluateAnimationStep(a,op,1).posCurrent,C.evaluateAnimationStep(b,op,1).posCurrent);
    assert.ok(Math.abs(C.vecNorm(delta)-C.vecNorm(C.vecSub(a,b)))<1e-7);
  }
  const pairOp=d.broken_operations.find(o=>C.findInCellAtomPairs(d,o).some(p=>p.mismatchAngstrom>1e-5));
  assert.ok(pairOp,'broken rigid operation must show a mismatch');
});
test('hydrostatic strain preserves ZnS; uniaxial strain gives a consistent subgroup', () => {
  const h=C.applyOrderParameterAndStrain(presets.zns_216,{exx:.04,eyy:.04,ezz:.04});
  assert.equal(h.preserved_operations.length,96);
  const u=C.applyOrderParameterAndStrain(presets.zns_216,{ezz:.08});
  assert.equal(u.isotropyAnalysis.subgroup.pointGroup.hm,'-42m');
  assert.equal(u.isotropyAnalysis.domainMultiplicity,3);
});
test('known tensor dimensions for ZnS, diamond and quartz', () => {
  for(const [key, counts] of Object.entries({zns_216:[0,1,1,3],silicon_227:[0,1,0,3],quartz_152:[0,2,2,6]})) {
    ['V','[V2]','V[V2]','[[V2]2]'].forEach((sym,i)=>assert.equal(C.computeJahnTensorFormat(presets[key].operations,sym).independentCount,counts[i],key+' '+sym));
  }
});

test('quartz twofold fixes Si and exchanges its four O neighbors in two pairs', () => {
  const d=presets.quartz_152;
  for(const atom of d.atoms.filter(a=>a.element==='Si')) {
    const site=C.getAtomSiteSymmetry(d,atom), op=site.operations.find(o=>o.type==='rotation');
    assert.equal(site.pointGroup.hm,'2'); assert.equal(site.neighbors.length,4);
    assert.equal(op.order,2);
    close(C.evaluateAnimationStep(atom.cart,op,1).posCurrent,atom.cart,'Si center fixed');
    const pairs=C.findSiteNeighborPairs(site,op,d.symprec);
    assert.equal(pairs.length,4); assert.ok(pairs.every(p=>p.moves&&p.atomA.element==='O'&&p.atomB.element==='O'));
    for(const pair of pairs) {
      close(C.evaluateAnimationStep(pair.atomA.cart,op,1).posFinal,pair.atomB.cart,'oxygen neighbor exchanged');
      close(C.evaluateAnimationStep(pair.atomB.cart,op,1).posFinal,pair.atomA.cart,'twofold swaps back');
    }
    // A quarter-turn keeps the center fixed too, but fails to reproduce the shell.
    const quarter=C.findSiteNeighborPairs(site,{...op,angle_deg:90},d.symprec);
    assert.ok(quarter.length<4,'stationary center alone must not imply fourfold symmetry');
  }
});

test('quartz oxygen sites have identity only; not every atom has a twofold site axis', () => {
  for(const a of presets.quartz_152.atoms.filter(a=>a.element==='O')) {
    const site=C.getAtomSiteSymmetry(presets.quartz_152,a);
    assert.equal(site.pointGroup.hm,'1'); assert.equal(site.operations.length,1);
    assert.equal(site.operations[0].type,'identity');
  }
});

test('site membership uses the selected Cartesian tolerance in a large cell', () => {
  const d=C.buildCrystalDataset({cellRows:[[100,0,0],[0,70,0],[0,0,60]],rawAtoms:[
    {element:'C',frac:[.002,.003,.2]},{element:'C',frac:[.998,.997,.2]},
  ],symprec:.005});
  const W=[[-1,0,0],[0,-1,0],[0,0,1]];
  assert.ok(d.operations.some(o=>JSON.stringify(o.W)===JSON.stringify(W)), 'twofold exchanges the two atoms');
  const site=C.getAtomSiteSymmetry(d,d.atoms[0]);
  assert.ok(!site.operations.some(o=>JSON.stringify(o.W)===JSON.stringify(W)), 'exchanging atoms is not site symmetry fixing one atom');
});

test('pair finder maps exact periodic images and reclassifies translated representatives', () => {
  for(const [key,d] of Object.entries(presets)) {
    const atoms=d.packed_atoms.slice(0,20);
    for(const a of atoms) for(const b of atoms) for(const op of C.findOperationsBetweenAtoms(d,a,b)) {
      const end=C.evaluateAnimationStep(a.cart,op,1).posFinal;
      assert.ok(C.vecNorm(C.vecSub(end,b.cart))<=d.symprec+2e-5,key+' pair endpoint');
      if(op.type==='identity') assert.ok(C.vecNorm(op.t_cart)<1e-8);
      if(op.type==='screw') assert.ok(C.vecNorm(op.t_parallel_cart)>1e-8);
      if(op.type==='glide') assert.ok(C.vecNorm(op.t_parallel_cart)>1e-8);
      assert.ok(!op.symbol.endsWith('_0'),'whole axial lattice shift is not a 2_0/3_0 screw');
    }
  }
  const d=presets.zns_216,a=d.atoms[0], b=d.packed_atoms.find(b=>b.element===a.element&&b.frac[0]===0&&b.frac[1]===0&&b.frac[2]===1);
  const translation=C.findOperationsBetweenAtoms(d,a,b).find(o=>o.type==='translation');
  assert.ok(translation); assert.match(translation.label,/Lattice Translation/);
  close(translation.w,[0,0,1]);
});

test('pair finder respects atomic tolerance and excludes broken parent operations', () => {
  const d=presets.zns_216;
  const offset={...d.atoms[0],frac:[.005,0,0],cart:C.matVec(d.cell_matrix_cols,[.005,0,0])};
  assert.equal(C.findOperationsBetweenAtoms({...d,symprec:.005},d.atoms[0],offset).length,0);
  const strained=C.applyOrderParameterAndStrain(d,{ezz:.08});
  for(const op of C.findOperationsBetweenAtoms(strained,strained.atoms[0],strained.atoms[1])) assert.notEqual(op.preserved,false);
});
