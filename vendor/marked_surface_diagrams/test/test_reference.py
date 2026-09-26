"""Regression tests: topology, faithful diagram data, fail-closed geometry.
No dependencies beyond Python 3.11 standard library.
"""
import copy,json,random,sys,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
import reference as R

def load(name): return R.read_document(ROOT/'examples'/f'{name}.json')
class KernelTests(unittest.TestCase):
    def bad(self,d,code):
        with self.assertRaises(R.DiagramError) as cm:R.validate(d)
        self.assertEqual(cm.exception.code,code)
    def test_01_all_examples_json_toml(self):
        for f in (ROOT/'examples').glob('*.json'):
            with self.subTest(f=f.stem):
                d=R.read_document(f);R.validate(d)
                self.assertEqual(d,R.read_document(f.with_suffix('.toml')))
    def test_02_sphere_and_positive_genera(self):
        for g in range(5):
            s=R.validate(load(f'genus_{g}_empty' if g else 'sphere_empty'))
            self.assertEqual(s.surface_report['chi'],2-2*g)
            r=R.regions(s);self.assertEqual(len(r),1);self.assertEqual((r[0]['genus'],r[0]['boundary_components']),(g,0))
    def test_03_fixed_endpoint_annulus_winding(self):
        anchors=[];marks=[];signed=[]
        for tag in ['m1','p0','p1']:
            d=load('annulus_winding_'+tag);s=R.validate(d)
            anchors.append(R.endpoint_anchors(s,'gamma'));marks.append(R.marking_digest(d))
            signed.append(sum(x['direction'] for x in R.seam_word(s,'gamma')))
        self.assertEqual(signed,[-1,0,1]);self.assertEqual(len({repr(a) for a in anchors}),1);self.assertEqual(len(set(marks)),1)
    def test_04_torus_crossing_and_cut_systems(self):
        s=R.validate(load('torus_mer_lon'))
        for c in ('alpha1','beta1'):self.assertEqual([(a['id'],a['sign']) for a in R.intersection_cycles(s)[c]],[('x',1)])
        for role in ('alpha','beta'):self.assertTrue(R.family_cut_report(s,role)['cut_system'])
    def test_05_normalized_disk_is_not_closed_embedded_disk(self):
        r=R.regions(R.validate(load('torus_mer_lon')))[0]
        self.assertTrue(r['normalized_disk']);self.assertFalse(r['embedded_closed_disk'])
    def test_06_parallel_zero_intersection_preserved(self):
        s=R.validate(load('torus_parallel_zero_crossings'))
        self.assertEqual(len(s.curves),2);self.assertEqual(len(s.crossing_at),0)
        self.assertEqual(sorted((r['genus'],r['boundary_components']) for r in R.regions(s)),[(0,2),(0,2)])
        self.assertTrue(R.family_cut_report(s,'alpha')['cut_system']);self.assertTrue(R.family_cut_report(s,'beta')['cut_system'])
    def test_07_nested_contractible_curves_preserved(self):
        s=R.validate(load('nested_zero_crossings'))
        self.assertEqual(sorted((r['genus'],r['boundary_components']) for r in R.regions(s)),[(0,1),(0,2),(0,2)])
        self.assertFalse(R.family_cut_report(s,'alpha')['cut_system'])
    def test_08_three_crossings_both_cycles(self):
        s=R.validate(load('torus_three_crossings'))
        for xs in R.intersection_cycles(s).values():self.assertEqual([(x['id'],x['sign']) for x in xs],[('x1',1),('x2',-1),('x3',1)])
        self.assertEqual(sum(r['embedded_closed_disk'] for r in R.regions(s)),2)
    def test_09_koda_orientation_not_silently_reversed(self):
        s=R.validate(load('torus_three_crossings'));old=copy.deepcopy(s.doc)
        self.assertFalse(any(R.check_joining_disk(s,r['seed'])['koda_joining'] for r in R.regions(s)))
        self.assertEqual(sum(R.check_joining_disk(s,r['seed'],{'alpha1':-1})['koda_joining'] for r in R.regions(s)),1)
        self.assertEqual(old,s.doc)
    def test_10_parser_printer_roundtrip(self):
        for f in (ROOT/'examples').glob('*.json'):
            for c in load(f.stem)['curves']:
                ds=R.parse_word(c['word']);self.assertEqual(ds,R.parse_word(R.print_word(ds)))
        for w in ['', 'P{}','P{a}','P{a b} junk','P{a b} || Q{x y}']:
            with self.subTest(w=w):
                with self.assertRaises(R.DiagramError):R.parse_word(w)
    def test_11_wrong_genus(self):
        d=load('torus_mer_lon');d['surface']['genus']=2;self.bad(d,'SURFACE_GENUS')
    def test_12_wrong_orientation(self):
        d=load('torus_mer_lon');d['patches'][0]['triangles'][0].reverse();self.bad(d,'EDGE_MANIFOLD')
    def test_13_duplicate_side_gluing(self):
        d=load('torus_mer_lon');d['surface']['gluings'].append(d['surface']['gluings'][0]);self.bad(d,'GLUING_DUPLICATE')
    def test_14_missing_route_edge(self):
        d=load('torus_mer_lon');d['curves'][0]['word']='P{v0_2 v4_2}';self.bad(d,'ROUTE_EDGE_MISSING')
    def test_15_missing_crossing(self):
        d=load('torus_mer_lon');d['crossings']=[];self.bad(d,'UNDECLARED_INTERSECTION')
    def test_16_wrong_crossing_sign(self):
        d=load('torus_mer_lon');d['crossings'][0]['sign']=-1;self.bad(d,'CROSSING_SIGN')
    def test_17_same_family_crossing(self):
        d=load('torus_mer_lon');d['curves'][1]['role']='alpha';self.bad(d,'SAME_FAMILY_INTERSECTION')
    def test_18_duplicate_mesh_triangle(self):
        d=load('torus_mer_lon');d['patches'][0]['triangles'].append(d['patches'][0]['triangles'][0]);self.bad(d,'TRIANGLE_DUPLICATE')
    def test_19_boundary_labels(self):
        d=load('annulus_winding_p0');d['surface']['boundary']=[];self.bad(d,'BOUNDARY_LABELS')
    def test_20_route_reversal_recomputes_sign(self):
        d=load('torus_mer_lon');ds=R.parse_word(d['curves'][0]['word']);d['curves'][0]['word']=R.print_word([R.Dart(a.face,a.v,a.u) for a in reversed(ds)])
        self.bad(d,'CROSSING_SIGN');d['crossings'][0]['sign']=-1;R.validate(d)
    def test_21_same_family_edge_overlap(self):
        d=load('torus_mer_lon');c=copy.deepcopy(d['curves'][0]);c['id']='alpha2';d['curves'].append(c);self.bad(d,'SHARED_ROUTE_EDGE')
    def test_22_exact_realization_all(self):
        for f in (ROOT/'examples').glob('*.json'):
            with self.subTest(f=f.stem):
                s=R.validate(load(f.stem));p=R.realize(s);self.assertTrue(R.certify_layout(s,p))
                self.assertTrue(all(isinstance(x,Q) for ps in p.values() for pt in ps.values() for x in pt))
    def test_23_invalid_layout_fails_closed(self):
        s=R.validate(load('torus_mer_lon'));p=R.realize(s);p['P']['v2_2']=p['P']['v1_1']
        with self.assertRaises(R.DiagramError):R.certify_layout(s,p)
    def test_24_certified_one_vertex_motion(self):
        s=R.validate(load('torus_mer_lon'));p=R.realize(s);old=p['P']['v2_2'];q,w=R.move_vertex(s,p,'P','v2_2',[str(old[0]+Q(1,10000)),str(old[1])])
        for t in [Q(k,20) for k in range(21)]:
            z=copy.deepcopy(p);z['P']['v2_2']=tuple((1-t)*a+t*b for a,b in zip(old,q['P']['v2_2']));R.certify_layout(s,z)
        self.assertEqual(w['same_representation'],R.representation_digest(s.doc))
    def test_25_boundary_motion_refused(self):
        s=R.validate(load('torus_mer_lon'));p=R.realize(s)
        with self.assertRaises(R.DiagramError) as cm:R.move_vertex(s,p,'P','v0_2',['0','0'])
        self.assertEqual(cm.exception.code,'MOVE_BOUNDARY')
    def test_26_decimal_svg_roundtrip_certified(self):
        s=R.validate(load('torus_three_crossings'));p=R.realize(s)
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/'x.svg';r=R.write_svg(s,p,f);self.assertTrue(r['certified_centerlines']);self.assertIn('centerlines_certified=true',f.read_text())
    def test_27_stellar_refinements_are_not_physical_punctures(self):
        d=load('torus_mer_lon');original=R.validate(copy.deepcopy(d));random.seed(620)
        for k in range(12):
            ts=d['patches'][0]['triangles'];i=random.randrange(len(ts));a,b,c=ts.pop(i);w=f'refine{k}'
            ts.extend([[a,b,w],[b,c,w],[c,a,w]])
            s=R.validate(d);self.assertEqual(R.intersection_cycles(s),R.intersection_cycles(original))
            self.assertEqual([(r['genus'],r['boundary_components']) for r in R.regions(s)],[(0,1)])
        self.assertEqual(R.marking_digest(d),R.marking_digest(original.doc));self.assertNotEqual(R.representation_digest(d),R.representation_digest(original.doc))
        R.realize(R.validate(d))
    def test_28_exact_parameter_checks(self):
        d=load('annulus_winding_p0');d['patches'][0]['side_parameters']={'east':['0','1/5','1/2','3/4','1']};self.bad(d,'GATE_PARAMETERS')
        d['patches'][0]['side_parameters']['west']=['0','1/4','1/2','4/5','1'];R.validate(d)
    def test_29_physical_endpoint_not_an_unlabelled_free_endpoint(self):
        d=load('annulus_winding_p0');d['curves'][0]['word']='P{v1_1 v1_2 v1_3}';self.bad(d,'ARC_ENDPOINT')
    def test_30_geometry_cannot_change_math_document(self):
        d=load('torus_three_crossings');old=copy.deepcopy(d);s=R.validate(d);R.realize(s);R.regions(s);self.assertEqual(d,old)

    def test_31_saved_region_mark(self):
        d=load('torus_three_crossings');d['region_marks']=[{'id':'D0','seed':['P',68],'claim':'embedded_closed_disk'}]
        R.validate(d)
        d['region_marks'][0]['claim']='koda_joining';self.bad(d,'REGION_CLAIM_FAILED')
    def test_32_mutating_validated_core_invalidates_state(self):
        s=R.validate(load('torus_mer_lon'));s.doc['curves'][0]['word']='P{v0_1 v1_1}'
        with self.assertRaises(R.DiagramError) as cm:R.realize(s)
        self.assertEqual(cm.exception.code,'MATHEMATICAL_STATE_MUTATED')
    def test_33_input_is_defensively_copied(self):
        d=load('torus_mer_lon');s=R.validate(d);d['curves']=[]
        self.assertEqual(len(s.curves),2);R.realize(s)
    def test_34_float_geometry_refused(self):
        s=R.validate(load('torus_mer_lon'));p=R.realize(s);p['P']['v2_2']=(0.0,0.0)
        with self.assertRaises(R.DiagramError) as cm:R.certify_layout(s,p)
        self.assertEqual(cm.exception.code,'LAYOUT_RATIONAL')
    def test_35_unknown_semantic_option_not_silently_ignored(self):
        d=load('torus_mer_lon');d['isotopy_canonical']=True;self.bad(d,'SCHEMA_KEYS')

    def test_36_legacy_projection_ordered_pair_sign(self):
        d=load('torus_mer_lon');a=R.legacy_projection(R.validate(d))
        d['crossings'][0]['curves'].reverse();d['crossings'][0]['sign']=-1
        b=R.legacy_projection(R.validate(d))
        self.assertEqual(a['intersections'],b['intersections'])
        self.assertEqual(a['alpha'],b['alpha']);self.assertEqual(a['beta'],b['beta'])

if __name__=='__main__':unittest.main(verbosity=2)
