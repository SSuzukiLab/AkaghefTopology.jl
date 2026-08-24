import fs from "node:fs";

const root = new URL("../artifacts/baseline/", import.meta.url);
const julia = JSON.parse(fs.readFileSync(new URL("julia_plot_states.json", root)));
const matlabV5 = JSON.parse(fs.readFileSync(new URL("matlab_replay_v5.json", root)));
const matlabV3 = JSON.parse(fs.readFileSync(new URL("matlab_replay_v3.json", root)));
const matlabV2 = JSON.parse(fs.readFileSync(new URL("matlab_replay_v2.json", root)));
const failures = [];
let compared = 0;

const scalarVector = value => Array.isArray(value) ? value : [value];
const findJulia = (group,index) => julia.find(record => record.group === group && record.index === index);
const recordsFor = (records,suffix) => records.filter(record => record.relativePath.includes(suffix));
const normalizeGauss = code => {
  const labels = [...new Set(code.flat().map(Math.abs).filter(Boolean))].sort((a,b) => a-b);
  const map = new Map(labels.map((label,index) => [label,index+1]));
  return code.map(component => component.map(vertex => Math.sign(vertex)*map.get(Math.abs(vertex))));
};

function compare(group,index,matlabObject,label) {
  const candidate = findJulia(group,index);
  const properties = matlabObject.properties;
  const expected = {
    rgauss: normalizeGauss(properties.RGaussCode),
    orientation: scalarVector(properties.ROrientation),
    isWeighted: properties.isWeighted,
    weights: properties.weightRE,
  };
  const actual = {
    rgauss: candidate?.rgauss,
    orientation: candidate?.orientation,
    isWeighted: candidate?.isWeighted,
    weights: candidate?.weights,
  };
  for (const key of ["rgauss","orientation","isWeighted",...(expected.isWeighted ? ["weights"] : [])]) {
    if (JSON.stringify(actual[key]) !== JSON.stringify(expected[key])) {
      failures.push(`${label} ${key}: MATLAB=${JSON.stringify(expected[key])} Julia=${JSON.stringify(actual[key])}`);
    }
  }
  compared += 1;
}

for (const record of recordsFor(matlabV5,"VLExample/call_")) {
  compare("VLExample",record.callId,record.objects.vl,`VLExample#${record.callId}`);
}

const c250818 = recordsFor(matlabV3,"C250818MSTinv_uqsl2/call_");
for (const record of c250818) {
  const objectName = new Map([[4,"vl2"],[6,"vl2"]]).get(record.callId) ?? "vl";
  compare("projects/invariant/C250818MSTinv_uqsl2",record.callId,
          record.objects[objectName],`C250818#${record.callId}`);
}

const c260125 = recordsFor(matlabV3,"C260125MSTinv_uqsl2/call_");
for (const record of c260125) {
  if (record.callId === 1) {
    record.objects.vl2.forEach((object,index) =>
      compare("projects/invariant/C260125MSTinv_uqsl2",index+1,object,`C260125#1[${index+1}]`));
    continue;
  }
  const juliaIndex = record.callId + 2;
  const objectName = new Map([[5,"vl2"],[8,"vl2"]]).get(record.callId) ?? "vl";
  compare("projects/invariant/C260125MSTinv_uqsl2",juliaIndex,
          record.objects[objectName],`C260125#${record.callId}`);
}

const singleMappings = [
  ["C260225M_Oinv/call_001","Docs/experiment/C260225M_Oinv",1,"vl"],
  ["C250811MSTinvUR/call_001","projects/invariant/C250811MSTinvUR",1,"vl"],
  ["C251020KnotDsl2bfromUR/call_001","projects/invariant/C251020KnotDsl2bfromUR",1,"vl"],
  ["C251020KnotDsl2bfromUR/call_002","projects/invariant/C251020KnotDsl2bfromUR",2,"vl"],
  ["C260220LensSpaces/call_001","projects/spine/example/C260220LensSpaces",1,"vl"],
];
for (const [suffix,group,index,objectName] of singleMappings) {
  const record = matlabV2.find(item => item.relativePath.includes(suffix));
  compare(group,index,record.objects[objectName],suffix);
}

if (failures.length) {
  console.error(failures.join("\n"));
  console.error(`FAILED: ${failures.length} field mismatches across ${compared} plot states`);
  process.exit(1);
}
console.log(`PASS: ${compared} MATLAB/Julia plot states matched`);
