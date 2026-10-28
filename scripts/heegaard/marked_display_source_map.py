"""Reconstruct and check whole marked-display source correspondences.

Source: marked_display.jl::prepare_marked_display supplies the oriented native
corner/edge quotient; layout_marked_document.coarse_view retains exact native
barycentric partitions and source paths. These are executable checks of the
full retained object, never trust in saved acceptance booleans.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F


# This bounds only private reconstruction data after decompression. It does
# not change the server's existing 5,000,000-byte request or native limits.
MAX_SOURCE_MAP_DECODED_BYTES=32*1024*1024
MAX_SOURCE_MAP_ENCODED_BYTES=5_000_000


def _source_map_json(value):
    import json
    return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode('utf-8')


def canonical_source_map_recipe(recipe):
    """Prove that JSON normalization changes only tuple/list containers.

    Source: Python json's array encoding. Check every actual node, scalar type,
    dictionary key and value; do not infer semantic equality from a successful
    encoder or silently accept key coercion. The returned recipe is JSON-native.
    """
    import json
    decoded=json.loads(_source_map_json(recipe))
    tuples=0
    def compare(a,b,path):
        nonlocal tuples
        if type(a) in (list,tuple):
            assert type(b) is list and len(a)==len(b),('Changed source-map array',path)
            tuples+=type(a) is tuple
            for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
        elif type(a) is dict:
            assert type(b) is dict and all(type(k) is str for k in a) and a.keys()==b.keys(),('Changed source-map object keys',path)
            for key in a:compare(a[key],b[key],path+'/'+key)
        else:
            assert type(a) in (str,int,float,bool,type(None)) and type(a) is type(b) and a==b,('Changed source-map scalar',path,type(a).__name__,type(b).__name__)
    compare(recipe,decoded,'')
    return decoded,dict(tuple_arrays_normalized=tuples,all_other_types_keys_values_identical=True)


def encode_source_map(recipe):
    """Bounded standard-zlib/Base85 transport; not a mathematical receipt.

    Source: Python base64.b85encode. The denser ASCII carrier preserves the
    exact zlib bytes; old zlib+base64 envelopes remain valid decoder inputs.
    """
    import base64,hashlib,zlib
    if not isinstance(recipe,dict) or recipe.get('schema')!='marked-display-source-map/1':
        raise ValueError('Unknown source-map reconstruction schema')
    raw=_source_map_json(recipe)
    if len(raw)>MAX_SOURCE_MAP_DECODED_BYTES:
        raise ValueError('Source-map reconstruction exceeds the 32 MiB decoded limit')
    data=base64.b85encode(zlib.compress(raw,level=9)).decode('ascii')
    if len(data)>MAX_SOURCE_MAP_ENCODED_BYTES:
        raise ValueError('Source-map envelope exceeds the existing request-size limit')
    return dict(schema='marked-display-source-map-envelope/1',codec='zlib+base85',
                uncompressed_bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),data=data)


def decode_source_map(envelope):
    """Decode one bounded, canonical recipe without accepting its claims.

    Source: Python zlib.decompressobj(max_length), canonical Base64/Base85 and
    JSON object_pairs_hook. Limit the decompressor BEFORE JSON allocation;
    reject trailing streams, duplicate keys, nonfinite numbers and length/hash
    mismatches. The caller must still rerun every source/map/geometric gate.
    """
    import base64,binascii,hashlib,json,re,zlib
    if not isinstance(envelope,dict) or set(envelope)!={'schema','codec','uncompressed_bytes','sha256','data'}:
        raise ValueError('Invalid source-map envelope fields')
    if envelope['schema']!='marked-display-source-map-envelope/1' or envelope['codec'] not in ('zlib+base64','zlib+base85'):
        raise ValueError('Unsupported source-map envelope')
    count=envelope['uncompressed_bytes'];data=envelope['data'];digest=envelope['sha256']
    if type(count) is not int or not 0<count<=MAX_SOURCE_MAP_DECODED_BYTES:
        raise ValueError('Invalid source-map decoded length')
    if not isinstance(data,str) or len(data)>MAX_SOURCE_MAP_ENCODED_BYTES:
        raise ValueError('Invalid source-map encoded length')
    if not isinstance(digest,str) or re.fullmatch('[0-9a-f]{64}',digest) is None:
        raise ValueError('Invalid source-map checksum')
    try:
        if envelope['codec']=='zlib+base64':
            compressed=base64.b64decode(data,validate=True)
            canonical=base64.b64encode(compressed).decode('ascii')
        else:
            compressed=base64.b85decode(data)
            canonical=base64.b85encode(compressed).decode('ascii')
        if canonical!=data:raise ValueError('Noncanonical source-map ASCII carrier')
        decoder=zlib.decompressobj()
        raw=decoder.decompress(compressed,count+1)
    except (ValueError,binascii.Error,zlib.error) as exc:
        raise ValueError('Invalid compressed source-map data') from exc
    if len(raw)!=count or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
        raise ValueError('Source-map length or compressed-stream boundary mismatch')
    if hashlib.sha256(raw).hexdigest()!=digest:
        raise ValueError('Source-map checksum mismatch')
    def object_pairs(pairs):
        result={}
        for key,value in pairs:
            if key in result:raise ValueError('Duplicate source-map JSON key')
            result[key]=value
        return result
    def bad_constant(value):raise ValueError('Nonfinite source-map JSON number')
    try:
        recipe=json.loads(raw,object_pairs_hook=object_pairs,parse_constant=bad_constant)
        if _source_map_json(recipe)!=raw:raise ValueError('Noncanonical source-map JSON')
    except (UnicodeError,RecursionError,OverflowError,TypeError) as exc:
        raise ValueError('Invalid source-map JSON') from exc
    if not isinstance(recipe,dict) or recipe.get('schema')!='marked-display-source-map/1':
        raise ValueError('Unknown source-map reconstruction schema')
    return recipe


def recode_source_map_carrier(envelope):
    """Losslessly update one validated envelope without recompressing geometry.

    Source: decode_source_map and Python base64. A retained accepted result can
    adopt the denser carrier without rerunning layout or mathematical checks;
    the source-map bytes, checksum and decoded-size declaration stay identical.
    This does not grant mathematical validity to an unverified envelope.
    """
    import base64
    decode_source_map(envelope)
    if envelope['codec']=='zlib+base85':return dict(envelope)
    compressed=base64.b64decode(envelope['data'],validate=True)
    return dict(envelope,codec='zlib+base85',data=base64.b85encode(compressed).decode('ascii'))


def check_source_map_saved_request_size(view):
    """Check the whole compact saved-view request, including ordinary fields.

    Source: existing /api/render's {diagram: view} input. The private envelope
    cap alone does not bound the surrounding source/positions/SVG. Julia repeats
    this check on its final returned object after adding the exact cell result.
    """
    size=len(_source_map_json(dict(diagram=view)))
    if not 0<size<=MAX_SOURCE_MAP_ENCODED_BYTES:
        raise ValueError(f'Saved marked-view request exceeds the existing 5000000-byte limit: {size}')
    return size


def verify_native_chart(d, p):
    """Bind a retained chart to a freshly prepared native marked document.

    Every native face is partitioned by a positive manifold disk with its full
    three-edge boundary, including the named infinity cap. Shared edge/corner
    keys, directed curve intervals, full D partition and crossing labels must
    agree. The caller must re-prepare p from the unchanged native source.
    """
    assert d['source_document']==p['source_document']
    assert d['preparation_representation_sha256']==p['source_representation_sha256']
    assert d['source_cut_corner_classes']==p['source_cut_corner_classes'] and d['boundary_vertex_pairs']==p['boundary_vertex_pairs']
    assert d['cuts']==p['cut_curves'] and d['visible_curves']==p['visible_curves'] and d['marking_id']==p['marking_id']
    assert d['infinity_source_triangle']==p['outer_triangle']==d['display_request']['infinity_triangle']
    assert d['disk_region_marks']==p['region_marks']
    native={tuple(t['source']):t for t in p['triangles']};groups=defaultdict(list);keys={};edge_endpoints={}
    orient=lambda a,b,c:(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    for n in native.values():
        for j,e in enumerate(n['edges']):
            u,v=n['vertices'][j],n['vertices'][(j+1)%3]
            if n['edge_signs'][j]<0:u,v=v,u
            assert e not in edge_endpoints or edge_endpoints[e]==(u,v)
            edge_endpoints[e]=(u,v)
    def native_key(n,row):
        nonzero=[i for i,x in enumerate(row) if x]
        if len(nonzero)==1:return ('vertex',n['vertices'][nonzero[0]])
        if len(nonzero)==2:
            j=next(j for j in range(3) if row[(j+2)%3]==0)
            return ('edge',n['edges'][j],row[(j+1)%3] if n['edge_signs'][j]>0 else row[j])
        return ('face',tuple(n['source']),tuple(row))
    for label,rows,vertices in zip(d['reduction_source_triangle_labels'],d['reduction_source_triangle_barycentric'],d['reduction_source_triangles'],strict=True):
        if label is None:assert rows is None;continue
        n=native[tuple(label)];q=[tuple(F(x) for x in row) for row in rows]
        assert all(sum(row)==1 and min(row)>=0 for row in q)
        triangle=[(row[1],row[2]) for row in q];assert orient(*triangle)>0
        groups[tuple(label)].append(triangle)
        for v,row in zip(vertices,q):
            key=native_key(n,row)
            assert v not in keys or keys[v]==key,('native seam point mismatch',v,keys.get(v),key)
            keys[v]=key
    inf=d['infinity_chart'];assert inf['source_triangle']==d['infinity_source_triangle']
    center=tuple(F(x) for x in inf['point_barycentric']);assert sum(center)==1 and min(center)>0
    frame=[tuple(F(x) for x in row) for row in inf['frame_barycentric']]
    assert len(frame)==len(d['chart_outer_boundary']) and all(sum(row)==1 and min(row)>0 for row in frame)
    for v,row in zip(d['chart_outer_boundary'],frame):assert keys[v]==native_key(native[tuple(inf['source_triangle'])],row)
    for a,b in zip(frame,frame[1:]+frame[:1]):
        triangle=[(row[1],row[2]) for row in (center,a,b)];assert orient(*triangle)>0
        groups[tuple(inf['source_triangle'])].append(triangle)
    assert set(groups)==set(native)

    # Each native simplex is covered by one positive manifold disk with precisely
    # its three oriented boundary intervals. This excludes overlaps hidden by an
    # area-only checksum and includes the retained infinity cap.
    partition_triangles=0
    for label,triangles in groups.items():
        edges=defaultdict(list);links=defaultdict(list);vertices={v for t in triangles for v in t}
        for t in triangles:
            assert orient(*t)>0;partition_triangles+=1
            for j,u in enumerate(t):
                v,w=t[(j+1)%3],t[(j+2)%3]
                edges[tuple(sorted((u,v)))].append((u,v));links[u].append((v,w))
        assert all(len(x)==1 or len(x)==2 and x[0]==x[1][::-1] for x in edges.values())
        boundary=[e[0] for e in edges.values() if len(e)==1];bv={v for e in boundary for v in e}
        assert len(vertices)-len(edges)+len(triangles)==1
        for u,link in links.items():
            graph=defaultdict(set)
            for v,w in link:graph[v].add(w);graph[w].add(v)
            visited=set();todo=[next(iter(graph))]
            while todo:
                v=todo.pop()
                if v in visited:continue
                visited.add(v);todo.extend(graph[v]-visited)
            assert visited==set(graph)
            ds=Counter(map(len,graph.values()));assert ds[1]==(2 if u in bv else 0) and all(k in (1,2) for k in ds)
        intervals=defaultdict(list)
        for a,b in boundary:
            if a[1]==b[1]==0:j=0;t,u=a[0],b[0]
            elif sum(a)==sum(b)==1:j=1;t,u=a[1],b[1]
            else:assert a[0]==b[0]==0;j=2;t,u=1-a[1],1-b[1]
            assert 0<=t<u<=1;intervals[j].append((t,u))
        assert set(intervals)=={0,1,2}
        for parts in intervals.values():
            parts.sort();assert parts[0][0]==0 and parts[-1][1]==1
            assert all(a[1]==b[0] for a,b in zip(parts,parts[1:]))
        assert sum(orient(*t) for t in triangles)==1

    paths=d['source_edge_paths'];assert set(map(int,paths))==set(edge_endpoints)
    for e,(u,v) in edge_endpoints.items():
        word=paths[str(e)];assert keys[word[0]]==('vertex',u) and keys[word[-1]]==('vertex',v)
        parameters=[F(0)]+[keys[x][2] for x in word[1:-1]]+[F(1)]
        assert all(keys[x][:2]==('edge',e) for x in word[1:-1])
        assert all(a<b for a,b in zip(parameters,parameters[1:]))
    boundaries={(b['cut'],b['side']):b for b in d['boundaries']}
    for b in p['boundaries']:
        current=boundaries[b['cut'],b['side']]
        for key in ('edges','edge_signs','vertices','source_vertices'):
            expected=b['vertices'] if key=='source_vertices' else b[key]
            assert current[key]==expected
        expected=[paths[str(e)][::s] for e,s in zip(b['edges'],b['edge_signs'])]
        assert current['source_boundary_paths']==expected
    ports={port['source_vertex']:(b,port) for b in d['boundaries'] for port in b['ports']}
    curve_count=0
    for c,b in zip(d['curves'],p['curves'],strict=True):
        for key in ('id','role','closed','word','seam_word'):assert c[key]==b[key]
        for s,t in zip(c['segments'],b['segments'],strict=True):
            for key in ('source','edge','edge_sign'):assert s[key]==t[key]
            expected=paths[str(t['edge'])][::t['edge_sign']]
            assert [expected[0],expected[-1]]==t['vertices']
            u,v=t['vertices']
            if u in ports:
                bd,port=ports[u];expected=[port['vertex']]+bd['source_cap_spokes'][str(u)][1:-1]+expected
            if v in ports:
                bd,port=ports[v];expected=expected+bd['source_cap_spokes'][str(v)][1:-1][::-1]+[port['vertex']]
            assert s['vertices']==expected;curve_count+=1
    patches={x['face']:x for x in d['source_document']['patches']}
    expected_crossings=[]
    for x in p['crossings']:
        f,v=x['at'];nodes={p['source_cut_corner_classes'][f'{f},{i},{j}'] for i,t in enumerate(patches[f]['triangles']) for j,w in enumerate(t) if w==v}
        assert len(nodes)==1;expected_crossings.append(dict(x,vertex=nodes.pop()))
    assert d['crossings']==expected_crossings
    pieces=[dict(mark=mark['id'],source_triangles=piece['triangles'],boundary_cycles=piece['display_boundary_cycles'],
                 contains_infinity=d['infinity_source_triangle'] in piece['triangles']) for mark in p['region_marks'] for piece in mark['cut_pieces']]
    legacy=(len(pieces)==1 and not pieces[0]['contains_infinity'] and all(part['kind']=='curve' for cycle in pieces[0]['boundary_cycles'] for part in cycle))
    if legacy:
        routes={tuple(s['source']):s['vertices'] for c in d['curves'] for s in c['segments']};disk=[]
        for f,_,u,v in p['region_marks'][0]['source_region']['boundary_cycles'][0]:
            route=routes[(f,u,v)] if (f,u,v) in routes else routes[(f,v,u)][::-1];disk.extend(route[:-1])
        assert d['disk']==disk and not d.get('disk_pieces')
    else:assert not d['disk'] and d.get('disk_pieces',[])==pieces

    return dict(native_faces=len(native), exact_partition_triangles=partition_triangles,
                native_edges=len(paths), shared_marked_source_points=len(keys),
                directed_curve_intervals=curve_count, crossings=len(expected_crossings),
                full_native_face_edge_and_seam_binding=True,
                explicit_infinity_cap_covered=True, D_partition_and_directed_boundary_bound=True)

def build_collar_sectors(d, z):
    """Build source/target collar disks without replacing a visible route.

    Source: coarse_view's native/cap labels and retained cap spokes, followed
    by fixed_mesh's explicit clipping. Triangle pCe only proposes each disk
    triangulation; exact full simplicity and both common-refinement coverages
    are verified by marked_display_source_map.jl::_mdsm_collar_sector.
    Boundary vertices are matched by the same ordered rational parameters.
    """
    import numpy as np
    import networkx as nx
    import triangle
    from shapely.geometry import Polygon
    def pt(p):return tuple(F(float(x)) for x in p)
    def blend(a,b,t):return tuple((1-t)*x+t*y for x,y in zip(a,b))
    def center(q):return tuple(sum(p[k] for p in q)/len(q) for k in (0,1))
    def area(q):return sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(q,q[1:]+q[:1]))/2
    def det(q):return 2*area(q)
    def edges(q):return zip(q,q[1:]+q[:1])
    def arc(cycle,a,b):
        i=cycle.index(a);out=[a]
        while out[-1]!=b:
            i=(i+1)%len(cycle);out.append(cycle[i]);assert len(out)<=len(cycle)
        return out
    def matched(a,b):
        grid=sorted({F(i,len(a)-1) for i in range(len(a))}|{F(i,len(b)-1) for i in range(len(b))})
        def sample(q,t):
            if t==1:return q[-1]
            z=t*(len(q)-1);i=z.numerator//z.denominator
            return blend(q[i],q[i+1],z-i)
        return [sample(a,t) for t in grid],[sample(b,t) for t in grid]
    def triangulate(q):
        assert area(q)>0 and len(q)==len(set(q))
        poly=Polygon([[float(x) for x in p] for p in q]);assert poly.is_valid
        P=np.array([[float(x) for x in p] for p in q]);segments=np.array([[i,(i+1)%len(q)] for i in range(len(q))])
        r=triangle.triangulate(dict(vertices=P,segments=segments),'pCe')
        assert np.array_equal(r['vertices'],P)
        T=r['triangles'].tolist();inc=Counter(tuple(sorted(e)) for t in T for e in zip(t,t[1:]+t[:1]))
        assert {e for e,n in inc.items() if n==1}=={tuple(sorted(e)) for e in segments}
        assert len(q)-len(inc)+len(T)==1
        assert all(det([q[v] for v in t])>0 for t in T)
        assert sum(area([q[v] for v in t]) for t in T)==area(q)
        return T
    def serial(points):return [[str(x) for x in p] for p in points]


    P={int(v):pt(q) for v,q in d['positions'].items()}
    T=d['reduction_source_triangles'];labels=d['reduction_source_triangle_labels']
    physical=[t for t,label in zip(T,labels) if label is not None]
    physical_edges={tuple(sorted((a,b))) for t in physical for a,b in edges(t)}
    visible_edges={tuple(sorted((a,b))) for c in d['curves'] for s in c['segments'] for a,b in zip(s['vertices'],s['vertices'][1:])}
    eps=F(1,4);results=[];report=[]
    for bi,b in enumerate(d['boundaries']):
        source_boundary=[v for word in b['source_boundary_paths'] for v in word[:-1]]
        B=set(source_boundary);BE={tuple(sorted(e)) for e in edges(source_boundary)}
        assert area([P[v] for v in source_boundary])<0
        boundary={};piece_area=F(0);available=defaultdict(set)
        def add_polygon(q):
            nonlocal_dummy=None
            assert area(q)>0
            for a,c in edges(q):
                if (c,a) in boundary:del boundary[c,a]
                else:assert (a,c) not in boundary;boundary[a,c]=1
        for t in physical:
            if not B.intersection(t):continue
            cen=center([P[v] for v in t])
            for j,u in enumerate(t):
                v=t[(j+1)%3];mid=blend(P[u],P[v],F(1,2));edge_on_B=tuple(sorted((u,v))) in BE
                for q,flags in (([P[u],mid,cen],[u in B,edge_on_B,False]),([mid,P[v],cen],[edge_on_B,v in B,False])):
                    indices=[i for i,x in enumerate(flags) if x]
                    if not indices:continue
                    if len(indices)==1:
                        i=indices[0];polygon=[q[i],blend(q[i],q[(i+1)%3],eps),blend(q[i],q[(i+2)%3],eps)]
                    else:
                        assert len(indices)==2
                        i=next(i for i in range(3) if flags[i] and flags[(i+1)%3]);j1=(i+1)%3;k=(i+2)%3
                        polygon=[q[i],q[j1],blend(q[j1],q[k],eps),blend(q[i],q[k],eps)]
                    add_polygon(polygon);piece_area+=area(polygon)
            for v in set(t)&B:available[v].add(blend(P[v],cen,eps))
        graph=nx.DiGraph();graph.add_edges_from(boundary)
        assert all(graph.in_degree(v)==graph.out_degree(v)==1 for v in graph)
        cycles=list(nx.simple_cycles(graph));assert len(cycles)==2
        outer=next(c for c in cycles if area(c)>0);inner=next(c for c in cycles if area(c)<0)
        assert area(outer)+area(inner)==piece_area
        assert set(P[v] for v in source_boundary)<=set(inner)
        bridges={};marked=[]
        for v in b['source_vertices']:
            ve=[e for e in visible_edges if v in e and e in physical_edges and e not in BE]
            assert len(ve)<=1,('Multiple visible branches at one cut port',v,ve)
            if ve:
                w=next(w for w in ve[0] if w!=v);q=blend(P[v],P[w],eps/2);marked.append(v)
            else:q=min(available[v])
            assert q in outer
            bridges[v]=q
        ring=[tuple(F(x) for x in z['source_positions'][str(v)]) for v in z['boundary_cycles'][bi]]
        if area(ring)>0:ring.reverse()
        cutpoints={}
        for v in b['source_vertices']:
            spoke=b['source_cap_spokes'][str(v)]
            cut=tuple(F(x) for x in b['target_cutpoints'][str(v)]) if 'target_cutpoints' in b else blend(P[b['cap']],P[spoke[1]],F(1,20))
            assert cut in ring;cutpoints[v]=cut
        source_sum=F(0);target_sum=F(0);sectors=[]
        vs=b['source_vertices']
        for i,v in enumerate(vs):
            w=vs[(i+1)%len(vs)];q,r=bridges[v],bridges[w]
            bpath=arc(inner,P[v],P[w]);cpath=arc(ring,cutpoints[v],cutpoints[w]);outerpath=arc(outer,r,q)
            assert not (set(bpath[1:-1]) & {P[x] for x in vs})
            assert not (set(cpath[1:-1]) & set(cutpoints.values()))
            assert not (set(outerpath[1:-1]) & set(bridges.values()))
            def extension(u):return [bridges[u]]+[P[x] for x in b['source_cap_spokes'][str(u)][:0:-1]]+[cutpoints[u]]
            left,right=extension(v),extension(w)
            source_chains=[[q,P[v]],bpath,[P[w],r],outerpath]
            target_chains=[left,cpath,right[::-1],outerpath]
            src=[];dst=[]
            for a,c in zip(source_chains,target_chains):
                aa,cc=matched(a,c);src.extend(aa[:-1]);dst.extend(cc[:-1])
            assert len(src)==len(dst)
            ts=triangulate(src);tt=triangulate(dst)
            source_sum+=area(src);target_sum+=area(dst)
            sectors.append(dict(boundary_edge=[v,w],source=serial(src),target=serial(dst),source_triangles=ts,target_triangles=tt,
                source_chains=[serial(x) for x in source_chains],target_chains=[serial(x) for x in target_chains]))
        assert source_sum==piece_area and target_sum==area(outer)+area(ring)
        results.append(dict(cut=b['cut'],side=b['side'],native_boundary_vertices=vs,visible_port_vertices=marked,
            source_outer=serial(outer),source_inner=serial(inner),target_inner=serial(ring),sectors=sectors))
        report.append(dict(cut=b['cut'],side=b['side'],sectors=len(sectors),vertices=sum(len(s['source']) for s in sectors),
            source_area=str(source_sum),target_area=str(target_sum),visible_ports=len(marked)))
    output=dict(diagnostic_only=True,layer=d.get('collar_layer','native_to_clipped_source'),exact_polygon_simplicity_pending=True,convex_overlay_pending=True,epsilon=str(eps),boundaries=results,report=report)

    return output


def circle_collar_input(d):
    """Retain the actual fixed-mesh radial extensions for a collar proposal.

    Source: marked_cut_layout.fixed_mesh's radius and directed endpoint/spoke
    construction. Visible ports are reused verbatim; all other radial samples
    use that same formula. These are polygonal domain coordinates, not a claim
    that a binary64 point lies exactly on the circle primitive.
    """
    import copy
    import numpy as np
    a=copy.deepcopy(d)
    a['collar_layer']='final_mesh_polygon_to_circle_polygon'
    a['reduction_source_triangles']=copy.deepcopy(d['triangles'])
    a['reduction_source_triangle_labels']=[['display_mesh',i] for i in range(len(d['triangles']))]
    z=dict(source_positions={},boundary_cycles=[]);n=max(map(int,d['positions']))+1
    for b in a['boundaries']:
        ring=b['mesh_boundary'][:];center=np.array(d['positions'][str(b['cap'])]);P={v:np.array(d['positions'][str(v)]) for v in ring}
        signed=sum(P[u][0]*P[v][1]-P[u][1]*P[v][0] for u,v in zip(ring,ring[1:]+ring[:1]))
        if signed>0:ring.reverse()
        b['source_vertices']=ring;b['source_boundary_paths']=[[u,v] for u,v in zip(ring,ring[1:]+ring[:1])]
        b['source_cap_spokes']={str(v):[b['cap'],v] for v in ring};b['target_cutpoints']={}
        actual={}
        for port in b['ports']:
            v=port['vertex'];neighbors=[]
            for c in d['curves']:
                for s in c['segments']:
                    vs=s['vertices']
                    if vs[0]==v:neighbors.append(vs[1])
                    if vs[-1]==v:neighbors.append(vs[-2])
            assert len(set(neighbors))==1 and neighbors[0] in ring
            actual[neighbors[0]]=d['positions'][str(v)]
        out=[]
        for v in ring:
            q=P[v]-center;point=actual.get(v,(center+b['radius']*q/np.linalg.norm(q)).tolist())
            assert abs(np.linalg.norm(np.array(point)-center)-b['radius'])<1e-12
            exact=[str(F(float(x))) for x in point];b['target_cutpoints'][str(v)]=exact
            z['source_positions'][str(n)]=exact;out.append(n);n+=1
        z['boundary_cycles'].append(out)

    return a,z


def build_port_domain(d, circle, oldseed):
    """Retain the full actual polygonal domain for port-transport verification.

    Source: fixed_mesh's radial port extensions and pair_ports' exact dyadic
    source-edge records. Added unmarked collar vertices may change generated
    indices; actual directed coordinates must still equal the emitted seed.
    The returned numerical inverse maps propose interior coordinates only.
    """
    import copy
    import numpy as np
    from scipy.interpolate import PPoly
    from marked_display_pipeline import pair_ports
    a=copy.deepcopy(d)
    P={int(v):tuple(map(F,q)) for v,q in a['positions'].items()}
    T=[list(t) for t in a['triangles']]
    sub=lambda x,y:(x[0]-y[0],x[1]-y[1])
    cross=lambda x,y:x[0]*y[1]-x[1]*y[0]
    orient=lambda x,y,z:cross(sub(y,x),sub(z,x))
    blend=lambda x,y,t:tuple(u+t*(v-u) for u,v in zip(x,y))
    nextid=max(P)+1;holes=[]
    for b,c in zip(a['boundaries'],circle['boundaries'],strict=True):
        assert (b['cut'],b['side'])==(c['cut'],c['side']), 'Circle collar label/side mismatch'
        ring=b['mesh_boundary'][:]
        actual={}
        for port in b['ports']:
            v=port['vertex']
            adjacent=[s['vertices'][1] if s['vertices'][0]==v else s['vertices'][-2] for curve in a['curves'] for s in curve['segments'] if v in (s['vertices'][0],s['vertices'][-1])]
            assert len(set(adjacent))==1
            actual[adjacent[0]]=v
        inner=[]
        for v in ring:
            point=tuple(F(x) for x in c['target_cutpoints'][str(v)])
            if v in actual:
                w=actual[v];assert P[w]==point
            else:
                w=nextid;nextid+=1;P[w]=point;a['positions'][str(w)]=list(map(float,point))
            inner.append(w)
        for u,v,U,V in zip(ring,ring[1:]+ring[:1],inner,inner[1:]+inner[:1]):
            for t in ([u,v,V],[u,V,U]):
                if orient(*(P[k] for k in t))<0:t.reverse()
                assert orient(*(P[k] for k in t))>0
                T.append(t)
        holes.append((b['cap'],inner))

    a['triangles']=T
    seed=pair_ports(a);Q={int(v):tuple(map(F,q)) for v,q in seed['positions'].items()}
    maps=[(np.array(m['center_svg']),m['inner_radius_svg'],m['outer_radius_svg'],[m['angle_origin']],PPoly(np.array(m['spline_coefficients']),np.array(m['spline_breaks']))) for m in seed['port_transport']['maps']]
    # Added unmarked collar vertices change generated IDs, never the drawn paths.
    for ca,cb in zip(seed['curves'],oldseed['curves'],strict=True):
        assert ca['id']==cb['id']
        for sa,sb in zip(ca['segments'],cb['segments'],strict=True):
            assert [seed['positions'][str(v)] for v in sa['vertices']]==[oldseed['positions'][str(v)] for v in sb['vertices']], 'Actual common seed geometry changed'


    for v,record in seed['port_transport']['vertex_sources'].items():
        u,w=record['edge'];P[int(v)]=blend(P[u],P[w],F(record['parameter']))
    domain=dict(source_positions={str(v):list(map(str,q)) for v,q in P.items()},
                triangles=T,holes=[dict(cap=cap,vertices=ring) for cap,ring in holes],seed=seed)
    return domain,(a,seed,oldseed,P,Q,T,holes,maps)

def build_port_map(data, seed, oldseed, P, Q, T, holes, maps):
    """Propose a CDT and certify its exact whole-step polygonal realization.

    Source: Triangle constrained refinement, SciPy brentq inverse proposals,
    and pair_ports' retained ordered source-edge parameters. Refined constraint
    vertices are exact on both paired edges. Positive local determinants and
    the full boundary gates below do not yet assert closed-quotient gluing.
    """
    import math
    import numpy as np
    import triangle
    from scipy.optimize import brentq
    sub=lambda a,b:(a[0]-b[0],a[1]-b[1])
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    orient=lambda a,b,c:cross(sub(b,a),sub(c,a))
    blend=lambda a,b,t:tuple(x+t*(y-x) for x,y in zip(a,b))
    for v,r in seed['port_transport']['vertex_sources'].items():
        u,w=r['edge'];P[int(v)]=blend(P[u],P[w],F(r['parameter']))
    incidence=defaultdict(list)
    for t in T:
        for u,v in zip(t,t[1:]+t[:1]):incidence[tuple(sorted((u,v)))].append((u,v))
    he={tuple(sorted((u,v))) for _,ring in holes for u,v in zip(ring,ring[1:]+ring[:1])}
    exterior=[e[0] for k,e in incidence.items() if len(e)==1 and k not in he]
    constraints=set(tuple(sorted(e)) for e in exterior)|he
    for c in seed['curves']:
        for s in c['segments']:constraints.update(tuple(sorted(e)) for e in zip(s['vertices'],s['vertices'][1:]))
    disk=seed.get('disk',[]);constraints.update(tuple(sorted(e)) for e in zip(disk,disk[1:]+disk[:1]))
    assert [seed['positions'][str(v)] for v in disk]==[oldseed['positions'][str(v)] for v in oldseed.get('disk',[])]
    for key in ('disk_pieces','disk_region_marks','crossings'):
        assert seed.get(key)==oldseed.get(key),('Changed marked feature',key)
    for x in seed['crossings']:
        assert seed['positions'][str(x['vertex'])]==oldseed['positions'][str(x['vertex'])]
    ids=sorted({v for e in constraints for v in e});idx={v:i for i,v in enumerate(ids)};edges=sorted(constraints)
    q=[[float(x) for x in Q[v]] for v in ids]
    inp=dict(vertices=q,segments=[[idx[u],idx[v]] for u,v in edges],segment_markers=[[i+1] for i in range(len(edges))])
    # Source: Triangle 20250106 triangle.tri.triangulate documents holes as an
    # optional two-dimensional array. Its Cython bridge indexes row zero even
    # for an empty array; omit that optional field for the no-cut domain.
    if holes:inp['holes']=[[float(x) for x in Q[c]] for c,_ in holes]
    result=triangle.triangulate(inp,'pqa0.002')
    assert np.array_equal(result['vertices'][:len(ids)],np.array(q))
    print('CDT proposal',len(result['vertices']),len(result['triangles']),flush=True)
    points={i:tuple(map(F,r)) for i,r in enumerate(result['vertices'])}
    source={i:P[v] for i,v in enumerate(ids)}
    # Every boundary/curve Steiner point uses the SAME exact edge parameter on
    # both sides. Project only the internal proof vertex to the actual straight
    # target constraint, so no rounded Triangle intersection replaces a route.
    constrained={};segments=[]
    for edge,marker in zip(result['segments'],result['segment_markers']):
        u,v=edges[int(marker[0])-1];a,b=Q[u],Q[v];axis=max(range(2),key=lambda k:abs(b[k]-a[k]))
        for i in map(int,edge):
            t=(points[i][axis]-a[axis])/(b[axis]-a[axis]);assert 0<=t<=1
            src,dst=blend(P[u],P[v],t),blend(a,b,t)
            if i in constrained:assert constrained[i]==(src,dst)
            constrained[i]=(src,dst)
        segments.append((tuple(map(int,edge)),int(marker[0])))
    for i,(src,dst) in constrained.items():
        if i in source:assert source[i]==src and points[i]==dst
        source[i]=src;points[i]=dst
    def inverse(q):
        q=np.array(list(map(float,q)))*600
        for center,r,support,src,f in maps:
            v=q-center;rho=np.linalg.norm(v)
            if rho<1e-10 or rho>=support:continue
            z=np.clip((rho-r)/(support-r),0,1);w=1-z*z*(3-2*z);a0=src[0];tau=2*math.pi
            fun=lambda a:a+w*(float(f(a))-a)
            lo=fun(a0);angle=lo+(math.atan2(v[1],v[0])-lo)%tau
            a=brentq(lambda a:fun(a)-angle,a0,a0+tau,xtol=1e-14)
            return tuple(map(F,(center+rho*np.array([math.cos(a),math.sin(a)]))/600))
        return tuple(map(F,q/600))
    for i in points:
        if i not in source:source[i]=inverse(points[i])
    tris=[list(map(int,t)) for t in result['triangles']]
    def low(u,v,w):
        e,f=sub(source[v],source[u]),sub(source[w],source[u])
        de,df=sub(sub(points[v],points[u]),e),sub(sub(points[w],points[u]),f)
        c,b,a=cross(e,f),cross(de,f)+cross(e,df),cross(de,df)
        candidates=[c,a+b+c]
        if a>0 and 0 < -b/(2*a) < 1:candidates.append(c-b*b/(4*a))
        return min(candidates)
    areas=[low(*t) for t in tris]
    report=dict(vertices=len(points),faces=len(tris),nonpositive=sum(x<=0 for x in areas),minimum=float(min(areas)),
                source_nonpositive=sum(orient(*(source[v] for v in t))<=0 for t in tris),target_nonpositive=sum(orient(*(points[v] for v in t))<=0 for t in tris))
    print(report,flush=True)
    assert min(areas)>0,('CDT whole-step failure',report)
    assert len(set(source.values()))==len(source) and len(set(points.values()))==len(points)
    # Verify connected oriented 2-manifold combinatorics, including every link.
    occurrences=defaultdict(list);links=defaultdict(list)
    for t in tris:
        for j,u in enumerate(t):
            v,w=t[(j+1)%3],t[(j+2)%3];occurrences[tuple(sorted((u,v)))].append((u,v));links[u].append((v,w))
    assert all(len(e)==1 or len(e)==2 and e[0]==e[1][::-1] for e in occurrences.values())
    boundary={k for k,e in occurrences.items() if len(e)==1};boundary_vertices={v for e in boundary for v in e}
    for u,link in links.items():
        graph=defaultdict(set)
        for v,w in link:graph[v].add(w);graph[w].add(v)
        seen=set();todo=[next(iter(graph))]
        while todo:
            v=todo.pop()
            if v in seen:continue
            seen.add(v);todo.extend(graph[v]-seen)
        assert seen==set(graph)
        degrees=Counter(map(len,graph.values()))
        assert degrees[1]==(2 if u in boundary_vertices else 0) and all(k in (1,2) for k in degrees)
    assert len(points)-len(occurrences)+len(tris)==1-len(holes)
    # Triangle subdivides boundary constraints with the SAME parameter on each
    # side; the whole-time geometric cycles equal their unrefined versions.
    boundary_markers={m for e,m in segments if tuple(sorted(e)) in boundary}
    expected_boundary=he|{tuple(sorted(e)) for e in exterior}
    assert {edges[m-1] for m in boundary_markers}==expected_boundary
    # Exact source/target position equality above binds the source domain to the
    # already verified fixed mesh plus polygon annuli; same source winding=1.
    def oldlow(u,v,w):
        ii=[idx[x] for x in (u,v,w)]
        return low(*ii)
    # All outside edges retain a convex enclosing polygon throughout the step.
    nextedge={u:v for u,v in exterior};start=next(iter(nextedge));outer=[start]
    while nextedge[outer[-1]]!=start:outer.append(nextedge[outer[-1]])
    assert len(outer)==len(exterior)
    for i,u in enumerate(outer):assert oldlow(outer[i-1],u,outer[(i+1)%len(outer)])>0
    boxes=[]
    for cap,ring in holes:
        assert P[cap]==Q[cap]
        # Add fixed centers only to the orientation evaluator, not the mesh.
        center_id=max(points)+1;points[center_id]=Q[cap];source[center_id]=P[cap]
        sign=1 if orient(P[cap],P[ring[0]],P[ring[1]])>0 else -1
        for u,v in zip(ring,ring[1:]+ring[:1]):assert low(center_id,idx[u],idx[v])>0 if sign==1 else low(center_id,idx[v],idx[u])>0
        del points[center_id];del source[center_id]
        pp=[chart[v] for chart in (P,Q) for v in ring];boxes.append([(min(q[k] for q in pp),max(q[k] for q in pp)) for k in range(2)])
        for u,v in exterior:
            for w in ring:assert oldlow(u,v,w)>0
    for i,b in enumerate(boxes):
        for c in boxes[i+1:]:assert any(b[k][1]<c[k][0] or c[k][1]<b[k][0] for k in range(2))
    report.update(whole_step_positive=True,initial_manifold_verified=True,whole_step_boundary_embedded=True,
                  all_actual_curve_and_D_constraints_retained=True,closed_surface_isotopy_complete=False)
    output=dict(report=report,source_positions={str(i):list(map(str,q)) for i,q in source.items()},
                target_positions={str(i):list(map(str,q)) for i,q in points.items()},triangles=tris,
                original_vertex_ids=ids,constraints=edges,refined_constraint_edges=segments)

    return output


def build_exact_circle_anchors(m, d, *, boundary_cycles):
    """Construct witness-only exact anchors for the decoded circle quotient.

    Source: the actual CDT, fixed boundary pairing and the renderer's binary64
    circle attributes. Rational Pythagorean points are only proposals until the
    complete determinant/boundary checks below pass. These exact points must
    NEVER be rounded into the emitted graph and called exact circle incidence.
    """
    import math
    P={int(v):tuple(F(x) for x in q) for v,q in m['target_positions'].items()};Q=dict(P)
    T=m['triangles'];idx={v:i for i,v in enumerate(m['original_vertex_ids'])}
    sub=lambda a,b:(a[0]-b[0],a[1]-b[1])
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    orient=lambda a,b,c:cross(sub(b,a),sub(c,a))

    edges=defaultdict(list)
    for t in T:
        for u,v in zip(t,t[1:]+t[:1]):edges[tuple(sorted((u,v)))].append((u,v))
    nextedge={e[0][0]:e[0][1] for e in edges.values() if len(e)==1}
    unseen=set(nextedge);cycles=[]
    while unseen:
        first=next(iter(unseen));cycle=[first]
        while nextedge[cycle[-1]]!=first:cycle.append(nextedge[cycle[-1]])
        unseen.difference_update(cycle);cycles.append(cycle)
    area=lambda c:sum(cross(P[u],P[v]) for u,v in zip(c,c[1:]+c[:1]))
    outer=next(c for c in cycles if area(c)>0);holes=[c for c in cycles if area(c)<0]
    assert len(holes)==len(d['boundaries'])
    records=[];used=set()
    assert len(boundary_cycles)==len(d['boundaries'])
    by_cap={int(b['cap']):b['vertices'] for b in boundary_cycles}
    assert len(by_cap)==len(boundary_cycles)
    for b in d['boundaries']:
        ports={idx[v['vertex']] for v in b['ports']}
        # Source: build_port_domain's actual ordered inner cycle, retained in
        # the CDT constraint set even when it contains no visible port.
        original_ring=[idx[v] for v in by_cap[b['cap']]]
        # fixed_mesh/build_port_domain order this cycle as the positive cap
        # disk boundary. The punctured surface induces the opposite direction.
        # Change only the comparison convention, preserving its first vertex
        # and the complete order; native edge lambda and mesh data stay fixed.
        assert len(original_ring)>=3 and len(set(original_ring))==len(original_ring) and area(original_ring)>0
        original_ring=[original_ring[0]]+original_ring[:0:-1]
        matches=[c for c in holes if set(original_ring)<=set(c)]
        assert len(matches)==1;ring=matches[0]
        actual=[v for v in ring if v in set(original_ring)]
        offset=actual.index(original_ring[0])
        assert actual[offset:]+actual[:offset]==original_ring
        assert ports<=set(ring) and not(set(ring)&used);used.update(ring)
        # Exact interpretation follows the ACTUAL decoded SVG primitive, including
        # its Float64 multiplication by600, rather than exact(normalized)*600.
        center=tuple(F(float(x)*600)/600 for x in d['positions'][str(b['cap'])]);r=F(float(b['radius'])*600)/600
        for v in ring:
            q=sub(P[v],center);t=F(math.tan(math.atan2(float(q[1]),float(q[0]))/2))
            a=r*(1-t*t)/(1+t*t);z=r*2*t/(1+t*t);Q[v]=(center[0]+a,center[1]+z)
            assert dot(sub(Q[v],center),sub(Q[v],center))==r*r
        records.append(dict(cut=b['cut'],side=b['side'],center=center,radius=r,ring=ring,ports=ports))
    pairs={u:v for a,b in d['boundary_vertex_pairs'] for u,v in ((a,b),(b,a))}
    port_pairs=[]
    for cut in d['cuts']:
        top,bottom=sorted([b for b in d['boundaries'] if b['cut']==cut],key=lambda b:b['row'])
        ca,cb=[tuple(F(float(x)*600)/600 for x in d['positions'][str(b['cap'])]) for b in (top,bottom)]
        bp={v['source_vertex']:v for v in bottom['ports']}
        for port in top['ports']:
            u=idx[port['vertex']];v=idx[bp[pairs[port['source_vertex']]]['vertex']]
            q=sub(Q[u],ca);Q[v]=(cb[0]+q[0],cb[1]-q[1]);port_pairs.append((u,v))
            assert sub(Q[u],ca)==(Q[v][0]-cb[0],cb[1]-Q[v][1])

    def low(u,v,w):
        e,f=sub(P[v],P[u]),sub(P[w],P[u]);de,df=sub(sub(Q[v],Q[u]),e),sub(sub(Q[w],Q[u]),f)
        c,b,a=cross(e,f),cross(de,f)+cross(e,df),cross(de,df);vs=[c,a+b+c]
        if a>0 and 0 < -b/(2*a) < 1:vs.append(c-b*b/(4*a))
        return min(vs)
    minimum=min(low(*t) for t in T);assert minimum>0,('Exact-circle PL correction folds',float(minimum))
    assert all(P[v]==Q[v] for v in outer)
    for b in records:
        center_id=max(P)+1;P[center_id]=Q[center_id]=b['center']
        for u,v in zip(b['ring'],b['ring'][1:]+b['ring'][:1]):assert low(center_id,v,u)>0
        del P[center_id];del Q[center_id]
        assert all(dot(sub(Q[v],b['center']),sub(Q[v],b['center']))==b['radius']**2 for v in b['ring'])
        for u,v in zip(outer,outer[1:]+outer[:1]):
            for w in b['ring']:assert low(u,v,w)>0
    boxes=[]
    for b in records:
        points=[chart[v] for chart in (P,Q) for v in b['ring']]
        boxes.append([(min(q[k] for q in points),max(q[k] for q in points)) for k in range(2)])
    for i,b in enumerate(boxes):
        for c in boxes[i+1:]:assert any(b[k][1]<c[k][0] or c[k][1]<b[k][0] for k in range(2))
    maxport=max((max(abs(x-y) for x,y in zip(P[v],Q[v])) for pair in port_pairs for v in pair),default=F(0))
    report=dict(exact_paired_ports=len(port_pairs),boundary_vertices=len(used),
                target_circle_numeric_contract='exact rational of decoded Float64 SVG center/radius, divided by600 for normalized coordinates',
                whole_step_minimum_double_area=str(minimum),whole_step_positive=True,
                maximum_port_coordinate_change_svg=float(maxport*600),
                closed_surface_quotient_complete=False,renderer_changed=False)
    output=dict(report=report,corrected_boundary_positions={str(v):list(map(str,Q[v])) for v in used},
                holes=[dict(cut=b['cut'],side=b['side'],center=list(map(str,b['center'])),radius=str(b['radius']),vertices=b['ring'],ports=sorted(b['ports'])) for b in records],
                port_pairs=port_pairs,outer=outer)

    return output


def build_graph_free_collars(m, a):
    """Construct and verify collar supports fixing the actual marked graph.

    Source: SciPy HalfspaceIntersection proposes facets; exact rational H/V
    equality, circle support and all constrained-edge intersections verify the
    proposal. The final Bezier graph is checked separately against this SAME
    support; fixing only a polyline is not enough for final output binding.
    """
    import itertools
    import numpy as np
    from scipy.spatial import HalfspaceIntersection
    P={int(v):tuple(F(x) for x in q) for v,q in m['target_positions'].items()}
    P.update({int(v):tuple(F(x) for x in q) for v,q in a['corrected_boundary_positions'].items()})
    sub=lambda a,b:(a[0]-b[0],a[1]-b[1])
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    count=defaultdict(int)
    for t in m['triangles']:
        for u,v in zip(t,t[1:]+t[:1]):count[tuple(sorted((u,v)))]+=1
    graph=[e for e,marker in m['refined_constraint_edges'] if count[tuple(sorted(e))]==2]
    reports=[];collars=[]
    for b in a['holes']:
        center=tuple(F(x) for x in b['center']);R=F(b['radius']);ports=set(b['ports'])
        normals=[sub(P[v],center) for v in b['vertices']]
        bounds=[R*R+(F(0) if v in ports else R*R/F(10**7)) for v in b['vertices']]
        half=np.array([[float(n[0]),float(n[1]),-float(z)] for n,z in zip(normals,bounds)])
        proposal=HalfspaceIntersection(half,np.zeros(2))
        vertices=set()
        for facet in proposal.dual_facets:
            for i,j in itertools.combinations(facet,2):
                n,q=normals[i],normals[j];den=cross(n,q)
                if not den:continue
                v=((bounds[i]*q[1]-bounds[j]*n[1])/den,(n[0]*bounds[j]-q[0]*bounds[i])/den)
                if all(dot(v,k)<=z for k,z in zip(normals,bounds)):vertices.add(v)
        vertices=sorted(vertices,key=lambda v:np.arctan2(float(v[1]),float(v[0])))
        assert len(vertices)>=3
        for u,v in zip(vertices,vertices[1:]+vertices[:1]):
            assert cross(u,v)>0
            # Every actual polygon edge is a supplied supporting halfspace. Together
            # with all-vertex feasibility this proves the exact H/V equality.
            assert any(dot(n,u)==z==dot(n,v) for n,z in zip(normals,bounds))
        # Circle containment follows from |n|=R and offset>=R^2. Equality occurs
        # only at the declared visible port, never an untracked boundary point.
        assert all(dot(n,n)==R*R and z>=R*R for n,z in zip(normals,bounds))
        box=[(min(v[k] for v in vertices),max(v[k] for v in vertices)) for k in range(2)]
        contacts=set();tested=0
        for u,v in graph:
            x,y=sub(P[u],center),sub(P[v],center)
            if any(max(x[k],y[k])<box[k][0] or min(x[k],y[k])>box[k][1] for k in range(2)):continue
            tested+=1;lo,hi=F(0),F(1);delta=sub(y,x)
            for n,z in zip(normals,bounds):
                slope=dot(n,delta);room=z-dot(n,x)
                if slope>0:hi=min(hi,room/slope)
                elif slope<0:lo=max(lo,room/slope)
                elif room<0:lo,hi=F(1),F(0);break
            if lo>hi:continue
            assert lo==hi and ((lo==0 and u in ports) or (lo==1 and v in ports)),('Visible graph meets collar support',b['cut'],b['side'],u,v,float(lo),float(hi))
            contacts.add(u if lo==0 else v)
        assert contacts==ports
        reports.append(dict(cut=b['cut'],side=b['side'],outer_vertices=len(vertices),exact_contacts=len(contacts),tested_graph_edges=tested,graph_fixed=True))
        collars.append(dict(cut=b['cut'],side=b['side'],center=b['center'],radius=b['radius'],ports=sorted(ports),
            outer_positions=[[str(x+o) for x,o in zip(v,center)] for v in vertices],
            halfspace_normals=[list(map(str,n)) for n in normals],halfspace_offsets=list(map(str,bounds))))
    out=dict(report=reports,collars=collars,
        radial_projection='On each ray, map [r_inner(theta), r_outer(theta)] affinely onto [R, r_outer(theta)]; identity beyond. At tangent contact rays both inner/outer equal R and the map is identity. Positive radial slope away from contacts, continuous extension at contacts.',
        quotient_correction='Use the recorded continuous degree-one lifts in coordinates (theta,s), rho=R+s*(r_outer(theta)-R). Send theta to (1-w(s))*theta+w(s)*g(theta), s unchanged, with w(0)=1,w(1)=0. Visible touching port rays are fixed. If there are no visible ports the support has strictly positive width on the whole circle and the same periodic formula applies, using the explicit native lambda-zero and psi lifts. Extend by identity outside. Strictly increasing angular lifts give positive orientation.',
        whole_visible_graph_pointwise_fixed=True,all_time_quotient_recipe_pending=True)

    return out


def build_boundary_trace(d, z, fixed, native, circle):
    """Compose complete native-edge parameters through both collar boundaries.

    Source: build_collar_sectors retains the ordered inner-chain correspondence;
    fixed_mesh retains the same source boundary cycle. Each interval is traced
    exactly and full target boundary coverage is required without reversal.
    The interiors of these SAME collar sectors are verified separately in Julia.
    """
    def pt(x):return tuple(F(a) for a in x)
    def fl(x):return tuple(F(float(a)) for a in x)
    def blend(a,b,t):return tuple((1-t)*x+t*y for x,y in zip(a,b))
    def cross(a,b):return a[0]*b[1]-a[1]*b[0]
    def sub(a,b):return tuple(x-y for x,y in zip(a,b))
    def parameter(p,a,b):
        v=sub(b,a);w=sub(p,a)
        if cross(v,w):return None
        k=0 if v[0] else 1
        return w[k]/v[k]
    def matched(a,b):
        grid=sorted({F(i,len(a)-1) for i in range(len(a))}|{F(i,len(b)-1) for i in range(len(b))})
        def sample(q,t):
            if t==1:return q[-1]
            z=t*(len(q)-1);i=z.numerator//z.denominator
            return blend(q[i],q[i+1],z-i)
        return [sample(a,t) for t in grid],[sample(b,t) for t in grid]
    def transport(chain,maps):
        out=[]
        for (s,p),(t,q) in zip(chain,chain[1:]):
            parts=[]
            for a,b,c,d in maps:
                u=parameter(a,p,q);v=parameter(b,p,q)
                if u is None or v is None:continue
                lo=max(F(0),min(u,v));hi=min(F(1),max(u,v))
                if lo>=hi:continue
                parts.append((lo,hi,blend(c,d,(lo-u)/(v-u)),blend(c,d,(hi-u)/(v-u))))
            parts.sort(key=lambda x:x[0]);assert parts and parts[0][0]==0 and parts[-1][1]==1
            for left,right in zip(parts,parts[1:]):assert left[1]==right[0] and left[3]==right[2]
            for lo,hi,a,b in parts:
                row=(s+(t-s)*lo,a)
                if out and out[-1][0]==row[0]:assert out[-1]==row
                else:out.append(row)
                out.append((s+(t-s)*hi,b))
        return out

    def maps_from_boundary(b):
        out=[]
        for s in b['sectors']:
            a,c=matched(list(map(pt,s['source_chains'][1])),list(map(pt,s['target_chains'][1])))
            out.extend((a[i],a[i+1],c[i],c[i+1]) for i in range(len(a)-1))
        return out

    P={int(k):fl(v) for k,v in d['positions'].items()};S={int(k):pt(v) for k,v in z['source_positions'].items()};R={int(k):fl(v) for k,v in fixed['positions'].items()}
    reports=[];result=[]
    expected_sides=[(b['cut'],b['side']) for b in d['boundaries']]
    assert len(set(expected_sides))==len(expected_sides)
    assert [(b['cut'],b['side']) for b in native['boundaries']]==expected_sides
    assert [(b['cut'],b['side']) for b in circle['boundaries']]==expected_sides
    for i,b in enumerate(d['boundaries']):
        ring=z['boundary_cycles'][i]
        fixedmap=[(S[u],S[v],R[u],R[v]) for u,v in zip(ring,ring[1:]+ring[:1])]
        native_map=maps_from_boundary(native['boundaries'][i]);circle_map=maps_from_boundary(circle['boundaries'][i]);words=[]
        for word in b['source_boundary_paths']:
            chain=[(F(j,len(word)-1),P[v]) for j,v in enumerate(word)]
            for mapping in (native_map,fixedmap,circle_map):chain=transport(chain,mapping)
            assert chain[0][0]==0 and chain[-1][0]==1 and all(a[0]<c[0] for a,c in zip(chain,chain[1:]))
            words.append(dict(source_vertices=[word[0],word[-1]],parameters=[str(t) for t,p in chain],circle_polygon_points=[[str(x) for x in p] for t,p in chain]))
        for a,c in zip(words,words[1:]+words[:1]):assert a['circle_polygon_points'][-1]==c['circle_polygon_points'][0]
        # Entire final polygon is covered, without reversal, by the traced words.
        target_edges=[(m[2],m[3]) for m in circle_map];coverage=[F(0) for _ in target_edges]
        for word in words:
            points=list(map(pt,word['circle_polygon_points']))
            for p,q in zip(points,points[1:]):
                matches=[]
                for j,(a,c) in enumerate(target_edges):
                    u=parameter(p,a,c);v=parameter(q,a,c)
                    if u is not None and v is not None and 0<=u<v<=1:matches.append((j,v-u))
                assert len(matches)==1,(b['cut'],p,q,matches)
                j,length=matches[0];coverage[j]+=length
        assert all(x==1 for x in coverage)
        report=dict(cut=b['cut'],side=b['side'],native_edges=len(words),trace_segments=sum(len(w['parameters'])-1 for w in words),exact_full_circle_polygon_coverage=True)
        result.append(dict(report=report,words=words));reports.append(report)
    output=dict(diagnostic_only=True,exact_boundary_trace=True,full_gluing_parameterization_pending=True,boundaries=result)

    return output


def build_boundary_gluing(d, fixed, trace):
    """Bind every paired native cut interval to a shared rational lambda.

    Source: prepare_marked_display's boundary_vertex_pairs and the exact full
    boundary trace. Float64 port replay is a diagnostic check only; the fixed
    quotient family uses the exact lambda reversal and normalized lifts.
    """
    import math
    import numpy as np
    from marked_display_pipeline import pair_ports
    seed=pair_ports(fixed)
    pairs={v:w for a,b in d['boundary_vertex_pairs'] for v,w in ((a,b),(b,a))}
    P={int(v):np.array(x) for v,x in fixed['positions'].items()};S={int(v):np.array(x) for v,x in seed['positions'].items()}
    output=[];reports=[]
    for cut in d['cuts']:
        indices=sorted([i for i,b in enumerate(d['boundaries']) if b['cut']==cut],key=lambda i:d['boundaries'][i]['row']);ti,bi=indices
        top,bottom=[d['boundaries'][i] for i in indices];vs=top['source_vertices'];n=len(vs)
        ports={x['source_vertex']:x for x in fixed['boundaries'][ti]['ports']}
        anchors=[j for j,v in enumerate(vs) if v in ports];visible_pairs=len(anchors)
        # Source: retained directed native boundary. With no marked ports its
        # first source vertex defines lambda zero; this is not a visible curve.
        assert vs
        if not anchors:anchors=[0]
        m=len(anchors)
        start=anchors[0];order=vs[start:]+vs[:start]
        anchors=[j for j,v in enumerate(order) if v in ports] if ports else [0]
        lam={};weights={}
        for j,a in enumerate(anchors):
            b=anchors[j+1] if j+1<m else n
            for k in range(a,b):lam[order[k]]=F(j,m)+F(k-a,m*(b-a));weights[order[k]]=F(1,m*(b-a))
        assert sum(weights.values())==1
        assert all(pairs[v] in bottom['source_vertices'] for v in order)
        bv=bottom['source_vertices'];assert all(order[(order.index(pairs[a])-1)%n]==pairs[b] for a,b in zip(bv,bv[1:]+bv[:1]))
        phase=None
        if ports:
            first=ports[order[0]]['vertex'];phase=math.atan2(*(S[first]-S[top['cap']])[::-1])
        maxerror=0.
        boundary_records=[]
        for i,b in zip(indices,(top,bottom)):
            side=1 if i==ti else -1;words=[];star_dets=[]
            cap=tuple(F(float(x)) for x in P[b['cap']])
            for w in trace['boundaries'][i]['words']:
                a,c=w['source_vertices']
                if side==1:l0=lam[a];delta=weights[a]
                else:
                    # Reversed native edge has precisely lambda(t)=lambda_top(1-t).
                    top_a=pairs[c];delta=-weights[top_a];l0=lam[top_a]+weights[top_a]
                    assert (l0-lam[pairs[a]])%1==0
                pars=[F(x) for x in w['parameters']];ls=[l0+t*delta for t in pars]
                points=[tuple(F(x) for x in q) for q in w['circle_polygon_points']]
                for a0,b0 in zip(points,points[1:]):
                    aa=tuple(x-y for x,y in zip(a0,cap));bb=tuple(x-y for x,y in zip(b0,cap));det=aa[0]*bb[1]-aa[1]*bb[0];assert det<0;star_dets.append(det)
                words.append(dict(source_vertices=w['source_vertices'],native_parameters=w['parameters'],common_lambda=[str(x) for x in ls],circle_polygon_points=w['circle_polygon_points']))
            for port in fixed['boundaries'][i]['ports']:
                native=port['source_vertex'];l=lam[native if side==1 else pairs[native]];angle=side*(phase-2*math.pi*float(l));expected=np.array([math.cos(angle),math.sin(angle)])*fixed['boundaries'][i]['radius']
                maxerror=max(maxerror,float(np.linalg.norm(S[port['vertex']]-S[b['cap']]-expected)))
            boundary_records.append(dict(side=b['side'],orientation_sign=side,words=words,strict_radial_angle_monotonicity_exact=True))
        assert maxerror<1e-12
        report=dict(cut=cut,native_paired_edges=n,visible_pairs=visible_pairs,all_edge_parameters_reflect_exactly=True,visible_port_float_replay_max=maxerror)
        reports.append(report);output.append(dict(report=report,phase=phase,angle_formula='orientation_sign * (phase - 2*pi*common_lambda)',boundaries=boundary_records))
    # Rounded PCHIP pieces are numerical proposals and need not form an exact
    # continuous map at their knots. The old diagnostic inverse-PCHIP formula
    # is deliberately not a production proof. build_fixed_quotient instead
    # uses the certified full PL steps and its explicit normalized lifts.
    result=dict(cuts=output,port_transport=seed['port_transport'],
        exact_scope='full native rational lambda, paired edge reversal and polygon-ray order',
        float_port_replay_is_diagnostic_only=True,
        fixed_quotient_family_required=True)

    return result


# Source: normalized SVG numeric readback; retain the actual binary64 scale step.
def build_rendered_mesh(m,anchors):
    S={int(i):tuple(F(x)*600 for x in q) for i,q in m['target_positions'].items()}
    R=dict(S);ids=m['original_vertex_ids'];index={v:i for i,v in enumerate(ids)}
    for i in range(len(ids)):R[i]=tuple(F(float(x/600)*600) for x in S[i])
    blend=lambda x,y,t:tuple(u+t*(v-u) for u,v in zip(x,y))
    assigned={}
    for edge,marker in m['refined_constraint_edges']:
        u,v=(index[x] for x in m['constraints'][marker-1]);axis=max(range(2),key=lambda k:abs(S[v][k]-S[u][k]))
        for w in edge:
            t=(S[w][axis]-S[u][axis])/(S[v][axis]-S[u][axis]);assert 0<=t<=1 and blend(S[u],S[v],t)==S[w]
            q=blend(R[u],R[v],t)
            if w in assigned:assert assigned[w]==q
            assigned[w]=q
    R.update(assigned);A=dict(R)
    for i,q in anchors['corrected_boundary_positions'].items():A[int(i)]=tuple(F(x)*600 for x in q)
    return S,R,A


def build_infinity_extension(d,mesh,fixed):
    """Attach the specified native infinity cap to the actual exterior map.

    Source: coarse_view.infinity_chart and fixed_mesh's homothetic exterior
    collar; this is the retained exact ha4_infinity_extension construction.
    On each native fan use reciprocal radius, the existing fixed PL map in
    the finite annulus, then rays through its ACTUAL exterior target vertices.
    The named native barycentric center maps to the compactification point.
    """
    pt=lambda q:tuple(F(x) for x in q)
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    orient=lambda a,b,c:cross(sub(b,a),sub(c,a))
    area=lambda q:sum(cross(a,b) for a,b in zip(q,q[1:]+q[:1]))
    P={int(v):tuple(F(float(x)) for x in q) for v,q in d['positions'].items()}
    S={int(v):pt(q) for v,q in mesh['source_positions'].items()}
    R={int(v):tuple(F(float(x)) for x in q) for v,q in fixed['positions'].items()}
    frame=d['chart_outer_boundary'];ordered=frame if area([P[v] for v in frame])>0 else frame[::-1]
    ext=mesh['outer'];assert len(ext)==len(ordered)>=3
    direct=not d['boundaries']
    assert bool(d['cuts'])==bool(d['boundaries'])
    if direct:
        assert d['genus']==0 and ext==ordered
        assert all(S[v]==P[v]==R[v] for v in frame)
        # Source: existing no-cut fixed_mesh identity and the common quotient's
        # zero-centered exterior fan. Its exact strict-star test is retained.
        origin=(F(0),F(0))
    else:
        origins=[tuple(2*x-y for x,y in zip(P[v],S[w])) for v,w in zip(ordered,ext)]
        assert len(set(origins))==1;origin=origins[0]
    target_origin=(F(0),F(0))
    assert all(orient(origin,S[a],S[b])>0 and orient(target_origin,R[a],R[b])>0 for a,b in zip(ext,ext[1:]+ext[:1]))
    inf=d['infinity_chart'];U=pt(inf['point_barycentric'][1:]);B=[pt(q[1:]) for q in inf['frame_barycentric']]
    assert len(B)==len(frame) and all(sum(map(F,q))==1 and min(map(F,q))>0 for q in [inf['point_barycentric']]+inf['frame_barycentric'])
    factors=[]
    for a,b,u,v in zip(B,B[1:]+B[:1],frame,frame[1:]+frame[:1]):
        native=orient(U,a,b);source=orient(origin,P[u],P[v]);assert native>0 and source<0
        factors.append(-source/native)
    assert inf['source_triangle']==d['display_request']['infinity_triangle']==fixed['display_request']['infinity_triangle']
    return dict(attachment_kind='direct_no_cut' if direct else 'homothetic_collar',source_triangle=inf['source_triangle'],
        infinity_point_native_barycentric=inf['point_barycentric'],native_frame_barycentric=inf['frame_barycentric'],
        source_frame_vertices=frame,source_plane_origin=list(map(str,origin)),
        native_frame_exterior_vertices=[ext[ordered.index(v)] for v in frame],
        source_exterior_vertices=ext,source_exterior_positions=[list(map(str,S[v])) for v in ext],
        target_exterior_positions=[list(map(str,R[v])) for v in ext],target_plane_origin=['0','0'],
        positive_jacobian_factors=list(map(str,factors)),
        fan_map='U+t*(B_i(s)-U) -> O+(Q_i(s)-O)/t; t=0 -> infinity',
        exterior_map='0<t<=1: O_target+(R_frame_i(s)-O_target)/t' if direct else '1/2<=t<=1: existing fixed PL map; 0<t<=1/2: O_target+(R_exterior_i(s)-O_target)/(2*t)',
        finite_frame_join_parameter='1',outer_collar_join_parameter='1' if direct else '1/2')


def build_fixed_quotient(m, a, g, support, d):
    """Check full native lambda and construct one fixed-quotient family.

    Source: exact boundary trace, CDT, actual SVG rounding mesh and graph-free
    sector collars. The chart base is C0=K0 R0 B, including generally nonidentity
    K0. Continuous lifts and every paired boundary interval are retained;
    endpoint gluing alone is never substituted for all-time compatibility.
    """
    from functools import cmp_to_key
    pt=lambda q:tuple(F(x) for x in q)
    P={int(v):pt(q) for v,q in m['source_positions'].items()};Q={int(v):pt(q) for v,q in m['target_positions'].items()}
    _,rendered,anchored=build_rendered_mesh(m,a)
    U={v:tuple(x/600 for x in q) for v,q in rendered.items()};R={v:tuple(x/600 for x in q) for v,q in anchored.items()}
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    sub=lambda a,b:(a[0]-b[0],a[1]-b[1])
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    blend=lambda a,b,t:tuple(x+t*(y-x) for x,y in zip(a,b))
    serial=lambda q:list(map(str,q))
    idx={v:i for i,v in enumerate(m['original_vertex_ids'])}

    def parameter(p,a,b):
        v,w=sub(b,a),sub(p,a)
        if cross(v,w):return None
        return w[0]/v[0] if v[0] else w[1]/v[1]

    def trace(lambdas,points,edges):
        out=[]
        for l0,l1,x,y in zip(lambdas,lambdas[1:],points,points[1:]):
            pieces=[]
            for u,v in edges:
                r,s=parameter(P[u],x,y),parameter(P[v],x,y)
                if r is None or s is None:continue
                lo,hi=max(F(0),min(r,s)),min(F(1),max(r,s))
                if lo>=hi:continue
                records=[]
                for t in (lo,hi):
                    z=(t-r)/(s-r)
                    records.append((l0+t*(l1-l0),blend(P[u],P[v],z),blend(Q[u],Q[v],z),blend(U[u],U[v],z),blend(R[u],R[v],z)))
                pieces.append((lo,hi,*records))
            pieces.sort(key=lambda z:z[0]);assert pieces and pieces[0][0]==0 and pieces[-1][1]==1
            for c,e in zip(pieces,pieces[1:]):assert c[1]==e[0] and c[3]==e[2]
            for lo,hi,c,e in pieces:
                if out and out[-1][0]==c[0]:assert out[-1]==c
                else:out.append(c)
                out.append(e)
        return out

    def clockwise(points,first):
        def local(q):return (dot(first,q),-cross(first,q))
        def upper(q):return q[1]>0 or q[1]==0 and q[0]>=0
        def compare(a,b):
            x,y=local(a),local(b)
            if upper(x)!=upper(y):return -1 if upper(x) else 1
            c=cross(x,y)
            return -1 if c>0 else 1 if c<0 else 0
        return sorted(points,key=cmp_to_key(compare))

    output=[]
    for cut in g['cuts']:
        name=cut['report']['cut'];boundaries=[];native_anchor_lambdas={}
        for record in cut['boundaries']:
            b=next(x for x in a['holes'] if x['cut']==name and x['side']==record['side']);ring=b['vertices'];edges=list(zip(ring,ring[1:]+ring[:1]));words=[]
            for word in record['words']:
                lambdas=list(map(F,word['common_lambda']));pars=list(map(F,word['native_parameters']))
                assert pars[0]==0 and pars[-1]==1
                assert all(l==lambdas[0]+s*(lambdas[-1]-lambdas[0]) for s,l in zip(pars,lambdas))
                chain=trace(lambdas,list(map(pt,word['circle_polygon_points'])),edges)
                words.append(dict(source_vertices=word['source_vertices'],lambda_values=[str(x[0]) for x in chain],source_points=[serial(x[1]) for x in chain],seed_points=[serial(x[2]) for x in chain],rendered_points=[serial(x[3]) for x in chain],anchor_points=[serial(x[4]) for x in chain]))
                for native,l in zip(word['source_vertices'],(lambdas[0],lambdas[-1])):
                    key=(record['side'],native)
                    if key in native_anchor_lambdas:assert native_anchor_lambdas[key]==l%1
                    native_anchor_lambdas[key]=l%1
            for x,y in zip(words,words[1:]+words[:1]):
                for field in ('source_points','seed_points','rendered_points','anchor_points'):assert x[field][-1]==y[field][0]
            center=pt(b['center']);radius=F(b['radius']);big=radius*F(6,5)
            assert all(dot(sub(chart[v],center),sub(chart[v],center))<big*big for chart in (P,Q,U,R) for v in ring)
            # The semantic circle center is the decoded SVG center. Recheck the
            # whole radial graph about THIS center through all three PL stages,
            # rather than inheriting star-shapedness about a nearby scaled center.
            for X,Y in ((P,Q),(Q,U),(U,R)):
                for u,v in edges:
                    x,y=sub(X[v],center),sub(X[u],center)
                    dx,dy=sub(sub(Y[v],center),x),sub(sub(Y[u],center),y)
                    c0,b0,a0=cross(x,y),cross(dx,y)+cross(x,dy),cross(dx,dy)
                    values=[c0,a0+b0+c0]
                    if a0>0 and 0 < -b0/(2*a0) < 1:values.append(c0-b0*b0/(4*a0))
                    assert min(values)>0,'Decoded-circle center is not a strict radial kernel throughout the PL family'
            collar=next(c for c in support['collars'] if c['cut']==name and c['side']==b['side'])
            assert all(dot(sub(pt(v),center),sub(pt(v),center))<big*big for v in collar['outer_positions'])
            boundaries.append(dict(cut=name,side=b['side'],center=b['center'],radius=b['radius'],initial_support_radius=str(big),words=words))
        top,bottom=sorted([b for b in d['boundaries'] if b['cut']==name],key=lambda b:b['row'])
        center=pt(next(h for h in a['holes'] if h['cut']==name and h['side']==top['side'])['center']);topanchors={}
        for port in top['ports']:
            l=native_anchor_lambdas[(top['side'],port['source_vertex'])];v=idx[port['vertex']];topanchors[sub(R[v],center)]=l
        visible_anchor=bool(topanchors)
        if not visible_anchor:
            radius=F(next(h for h in a['holes'] if h['cut']==name and h['side']==top['side'])['radius'])
            topanchors[(radius,F(0))]=F(0)
        first=next(q for q,l in topanchors.items() if l==0)
        # Quarter-turn auxiliary points make every radial-chord interval < pi;
        # no transcendental angle or unrecorded integer winding is introduced.
        extra={first,(-first[1],first[0]),(-first[0],-first[1]),(first[1],-first[0])}
        order=clockwise(set(topanchors)|extra,first);assert order[0]==first
        lam=dict(topanchors);anchor_indices=[i for i,q in enumerate(order) if q in topanchors]
        for j,i in enumerate(anchor_indices):
            end=anchor_indices[j+1] if j+1<len(anchor_indices) else len(order)
            l0=topanchors[order[i]];l1=topanchors[order[end]] if end<len(order) else F(1)
            assert l1>l0
            for k in range(i+1,end):lam[order[k]]=l0+(l1-l0)*F(k-i,end-i)
        assert all(cross(u,v)<0 for u,v in zip(order,order[1:]+order[:1]))
        psi=dict(lambda_values=[str(lam[q]) for q in order]+['1'],top_relative_points=[serial(q) for q in order]+[serial(first)],
                 interpolation='R * normalize((1-s)*A+s*B), s affine in the recorded lambda interval; bottom is exact y reflection of top.')
        # The first anchor remains in one open semicircle during BOTH linear PL
        # steps. This supplies the unique time-continuous lift ending at angle0.
        pairing={u:v for e in d['boundary_vertex_pairs'] for u,v in (e,e[::-1])}
        byside={b['side']:b for b in boundaries}
        topwords={tuple(w['source_vertices']):w for w in byside[top['side']]['words']}
        bottomwords={tuple(w['source_vertices']):w for w in byside[bottom['side']]['words']}
        for (u,v),word in topwords.items():
            reverse=bottomwords[pairing[v],pairing[u]]
            x,y=map(F,(word['lambda_values'][0],word['lambda_values'][-1]))
            a0,b0=map(F,(reverse['lambda_values'][0],reverse['lambda_values'][-1]))
            assert (x-b0)%1==0 and (y-a0)%1==0 and y-x==a0-b0
        for b in (top,bottom):
            record=byside[b['side']]
            c=pt(next(h for h in a['holes'] if h['cut']==name and h['side']==b['side'])['center'])
            if visible_anchor:
                port=next(v for v in b['ports'] if native_anchor_lambdas[(b['side'],v['source_vertex'])]==0)
                v=idx[port['vertex']];points=[chart[v] for chart in (P,Q,U,R)]
                identity=dict(vertex=v)
            else:
                candidates=[(w,j) for w in record['words'] for j,l in enumerate(w['lambda_values']) if F(l)%1==0]
                assert candidates
                points=[pt(candidates[0][0][field][candidates[0][1]]) for field in ('source_points','seed_points','rendered_points','anchor_points')]
                assert all([pt(w[field][j]) for field in ('source_points','seed_points','rendered_points','anchor_points')]==points for w,j in candidates)
                identity=dict(native_lambda='0',visible_port=False)
            final=sub(points[-1],c)
            assert all(dot(final,sub(q,c))>0 for q in points)
            target=first if b['side']==top['side'] else (first[0],-first[1])
            record['lift_anchor']=dict(**identity,source=serial(points[0]),seed=serial(points[1]),rendered=serial(points[2]),final=serial(points[3]),
                psi_relative_ray=serial(target),psi_lift_from_final_ray=dict(det=str(cross(final,target)),dot=str(dot(final,target))),
                rule='The native lambda-zero ray uses atan2(det(final,ray),dot(final,ray)) in (-pi/2,pi/2), ending at0. The psi lambda-zero ray uses atan2(recorded det,recorded dot) in (-pi,pi], choosing +pi on the negative axis. Extend each lift once around its recorded oriented boundary; no additional integer turn is permitted.')
        output.append(dict(cut=name,boundaries=boundaries,canonical_paired_boundary=psi))

    # Disjoint supports throughout: all variable radial profiles lie below these
    # fixed big radii. They also leave the fixed infinity/exterior neighborhood.
    bs=[b for c in output for b in c['boundaries']]
    for i,b in enumerate(bs):
        for c in bs[i+1:]:assert dot(sub(pt(b['center']),pt(c['center'])),sub(pt(b['center']),pt(c['center'])))>(F(b['initial_support_radius'])+F(c['initial_support_radius']))**2
    def quadratic_min(c,b,a):
        values=[c,a+b+c]
        if a>0 and 0 < -b/(2*a) < 1:values.append(c-b*b/(4*a))
        return min(values)
    outer=a['outer'];outer_minimum=None
    markers={tuple(sorted(e)):marker for e,marker in m['refined_constraint_edges']}
    for u,v in zip(outer,outer[1:]+outer[:1]):
        for X,Y in ((P,Q),(Q,U),(U,R)):
            e,de=sub(X[v],X[u]),sub(sub(Y[v],Y[u]),sub(X[v],X[u]));du=sub(Y[u],X[u])
            for b in bs:
                h=sub(pt(b['center']),X[u]);dh=(-du[0],-du[1])
                lower=quadratic_min(cross(e,h),cross(de,h)+cross(e,dh),cross(de,dh))
                length2=max(dot(e,e),dot(sub(Y[v],Y[u]),sub(Y[v],Y[u])))
                assert lower>0 and lower*lower>F(b['initial_support_radius'])**2*length2
            # Fixed origin supports an explicit homogeneous extension to infinity.
            lower=quadratic_min(cross(X[u],X[v]),cross(du,X[v])+cross(X[u],sub(Y[v],X[v])),cross(du,sub(Y[v],X[v])))
            assert lower>0
            outer_minimum=lower if outer_minimum is None else min(outer_minimum,lower)
        original_edge=m['constraints'][markers[tuple(sorted((u,v)))]-1]
        old_a,old_b=[tuple(map(F,d['positions'][str(w)])) for w in original_edge]
        for w in (u,v):
            fraction=parameter(P[w],old_a,old_b)
            assert fraction is not None and 0<=fraction<=1
    recipe=dict(cuts=output,
        numeric_decoder='IEEE binary64 round-trip SVG parse; exact binary rational interpretation of decoded circle centers/radii and Bezier controls',
        parameter_time='Three certified straight vertex steps: CDT source to seed, seed to actual SVG-rounded mesh, then rounded mesh to exact anchors. The retained boundary traces include every stage.',
        radial_support='Use the recorded big circle during the first two steps. During the anchor step interpolate its radial profile linearly to the graph-free final outer polygon. Convexity of norm and exact final radius ensure this stays outside both moving inner polygon and exact circle, except the declared final port contacts.',
        R_t='Map each radial interval [r_inner(t,theta),r_outer(t,theta)] affinely onto [R,r_outer(t,theta)], identity outside. Extend continuously by identity at final touching port rays.',
        phi_t='Full recorded native-lambda boundary trace through the PL steps and R_t; use the recorded continuous anchor lift.',
        g_t='psi composed with inverse(phi_t), using those degree-one lifts. Both paired sides use the SAME recorded psi with exact reflection.',
        K_t='In coordinates rho=R+s*(r_outer(t,theta)-R), send (theta,s) to ((1-w(s))*theta+w(s)*g_t(theta),s), w(s)=1-s. Identity outside. At final touching rays g_1 fixes the anchor, so the formula extends continuously.',
        F_t='K_t composed with R_t composed with f_t composed with B, where B is the retained full native chart/collar/infinity map.',
        compatibility='For every native lambda and every time, F_t boundary equals the SAME psi(lambda) on top and its exact reflection on bottom. Thus F_t descends to one fixed quotient, and F_t composed with inverse(F_0) is an ambient isotopy there.',
        final_graph='R_1 and K_1 are supported in the verified graph-free collar and fix the exact-anchor Bezier graph pointwise. The NEW emitted endpoint-normalized Float64 Beziers are connected by the certified control/boundary homotopy. Old endpoint-gapped output is not claimed isotopic to a joined graph.',
        base_reference='C0=K_0 composed with R_0 composed with B (f_0=id). Generally g_0 is NOT identity. Retain this complete initial collar action, native cell/edge/gluing marking and normalized lifts. F_t composed with inverse(C0) starts at identity on the same psi quotient. C0 is not automatically a previously chosen reference chart.',
        exterior_extension=dict(origin=['0','0'],outer_stages=[[serial(chart[v]) for v in outer] for chart in (P,Q,U,R)],whole_step_positive_wedge_bound=str(outer_minimum),
            formula='For each exterior wedge, x=r*((1-s)*A+s*B), r>=1, maps to r*((1-s)*A_t+s*B_t). Positive star wedges give a homogeneous orientation-preserving extension fixing infinity; circle-collar supports stay strictly inside the moving exterior polygon.'),
        status='analytic family construction with finite rational boundary/lift gates; final normalized Bezier transfer certified separately; full native composition/common-entry replay pending')

    return recipe


def build_literal_boundary_paths(m, anchors, fixed, family):
    """Bind renderer L chains to the same complete directed native intervals.

    Source: the actual SVG-rounding mesh and fixed-quotient native trace. Only
    forward collinear subdivisions can disappear. Full wraps and all native
    darts are retained so the renderer cannot choose a different D complement.
    """
    _,R,_=build_rendered_mesh(m,anchors);ids=m['original_vertex_ids'];index={v:i for i,v in enumerate(ids)}
    sub=lambda x,y:(x[0]-y[0],x[1]-y[1])
    cross=lambda x,y:x[0]*y[1]-x[1]*y[0]
    dot=lambda x,y:x[0]*y[0]+x[1]*y[1]
    def simplify(q):
        out=[]
        for v in q:
            if out and out[-1]==v:continue
            while len(out)>1 and cross(sub(out[-1],out[-2]),sub(v,out[-1]))==0:
                assert dot(sub(out[-1],out[-2]),sub(v,out[-1]))>0;out.pop()
            out.append(v)
        return out
    def arc(ring,a,b,full=False):
        i=ring.index(a);out=[a]
        while out[-1]!=b or full and len(out)==1:
            i=(i+1)%len(ring);out.append(ring[i]);assert len(out)<=len(ring)+1
        return out
    records=[]
    for hole in anchors['holes']:
        cut,side=hole['cut'],hole['side'];ring=hole['vertices']
        original=[v for v in ring if v<len(ids)];assert set(hole['ports'])<=set(original)
        # All retained coordinates serialize exactly under the declared decoder.
        assert all(F(float(x))==x for v in original for x in R[v])
        for u,v in zip(original,original[1:]+original[:1]):
            assert simplify([R[w] for w in arc(ring,u,v)])==[R[u],R[v]]
        bd=next(b for b in fixed['boundaries'] if (b['cut'],b['side'])==(cut,side))
        trace=next(b for c in family['cuts'] for b in c['boundaries'] if (b['cut'],b['side'])==(cut,side))
        words={tuple(w['source_vertices']):w for w in trace['words']}
        native=bd['source_vertices'];ports={v['source_vertex']:index[v['vertex']] for v in bd['ports']}
        source_ports=[v for v in native if v in ports]
        interval_records=[]
        for u,v in zip(source_ports,source_ports[1:]+source_ports[:1]):
            source_word=arc(native,u,v,full=u==v);points=[]
            for x,y in zip(source_word,source_word[1:]):
                q=[tuple(F(t)*600 for t in z) for z in words[x,y]['rendered_points']]
                assert not points or points[-1]==q[0];points.extend(q if not points else q[1:])
            literal=arc(original,ports[u],ports[v],full=u==v)
            assert simplify(points)==simplify([R[w] for w in literal]),('Native cut interval does not bind literal boundary',cut,side,u,v)
            interval_records.append(dict(from_source_vertex=u,to_source_vertex=v,native_vertices=source_word,
                                         polygon_vertices=[original.index(w) for w in literal]))
        records.append(dict(cut=cut,side=side,orientation='clockwise',
            points_svg=[[float(x) for x in R[v]] for v in original],
            source_ports={str(v):original.index(w) for v,w in ports.items()},
            native_boundary_vertices=native,native_boundary_darts=bd['source'],consecutive_port_intervals=interval_records))
    output=dict(schema='marked-display-attachment-boundaries/1',numeric_decoder='IEEE binary64 round-trip',
                role='literal polygonal cut-surface domain, mapped by the retained chart to decoded-circle quotient; not a claim of circle-glyph point-set incidence',
                boundaries=records)

    return output


def build_exact_anchor_beziers(original, oldseed, m, a, *, minimum_separation=12., puncture_clearance=12., crossing_radius=20.):
    """Retain exact witness controls and check their numerical proposal.

    Source: the exact decoded-circle anchors, the source-directed Bezier end
    controls, and existing smooth/fairing validators. The original emitted
    graph is never overwritten. Exact rational replacements are consumed by
    verify_bezier_attachment, not rounded and claimed to lie on a circle.
    """
    import copy
    import numpy as np
    from check_smooth_marked import check,leaves
    from check_fairing_isotopy import check as check_motion
    result=copy.deepcopy(original);seed=copy.deepcopy(oldseed)
    anchor={};centers={};radii={}
    for hole in a['holes']:
        for i in hole['ports']:
            v=m['original_vertex_ids'][i]
            anchor[v]=tuple(F(x)*600 for x in a['corrected_boundary_positions'][str(i)])
            centers[v]=tuple(F(x)*600 for x in hole['center']);radii[v]=F(hole['radius'])*600
            seed['positions'][str(v)]=[float(x/600) for x in anchor[v]]
    replacements=[];maximum=F(0)
    for curve in result['curves']:
        for segment in curve['segments']:
            for is_start,v in ((True,segment['vertices'][0]),(False,segment['vertices'][-1])):
                if v not in anchor:continue
                bi=0 if is_start else len(segment['beziers'])-1;ei,ci=(0,1) if is_start else (3,2)
                bez=segment['beziers'][bi];endpoint,control=[tuple(map(F,bez[k])) for k in (ei,ci)]
                normal=tuple((x-y)/radii[v] for x,y in zip(anchor[v],centers[v]));assert sum(x*x for x in normal)==1
                distance=sum((x-y)*n for x,y,n in zip(control,endpoint,normal));assert distance>0
                corrected_control=tuple(x+distance*n for x,n in zip(anchor[v],normal))
                for k,q in ((ei,anchor[v]),(ci,corrected_control)):
                    maximum=max(maximum,max(abs(x-F(y)) for x,y in zip(q,bez[k])))
                    bez[k]=list(map(float,q))
                    replacements.append(dict(curve=curve['id'],source=segment['source'],piece=bi,control=k,exact_svg=list(map(str,q))))
            pieces=[np.array(q) for q in segment['beziers']]
            points=[pieces[0][0]]+[leaf[-1] for q in pieces for leaf in leaves(q)]
            assert len(points)==len(segment['vertices'])
            for v,q in zip(segment['vertices'],points):result['positions'][str(v)]=(q/600).tolist()
    assert maximum<F(1,10**7), 'Exact endpoint correction exceeds the existing isotopy tolerance'
    result=check(result,seed,minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,crossing_radius=crossing_radius)
    motion=check_motion(seed,result)
    assert result['source_document']==original['source_document'] and result['display_request']==original['display_request']
    return dict(replacements=replacements,report=dict(exact_endpoint_replacements=len(anchor),
        control_replacements=len(replacements),maximum_control_coordinate_change_svg=str(maximum),
        same_geometry_thresholds=True,smooth_check=result['smooth_check'],layout_isotopy_check=motion,
        rounded_candidate_is_diagnostic_only=True,emitted_output_unchanged=True))


def verify_bezier_attachment(d, seed, m, a, replacement, collar_data):
    """Check the actual normalized control graph through its attachment chart.

    Source: existing dyadic Bezier subdivision and exact quadratic orientation,
    swept control hulls, moving port planes and the SAME graph-free support.
    This certifies the encoded polygonal-domain graph, not point-set incidence
    with a separately displayed circle glyph or hardware rasterization.
    """
    import math,time
    import numpy as np
    from shapely.geometry import box
    from shapely.strtree import STRtree
    started=time.monotonic()
    normalization=d['source_endpoint_normalization'];shared=normalization['shared_svg_endpoints']
    fixed_vertices={x['vertex'] for x in seed['crossings']}|{v['vertex'] for b in seed['boundaries'] for v in b['ports']}
    for curve,original in zip(d['curves'],seed['curves'],strict=True):
        for key in ('id','role','closed','word','seam_word'):assert curve[key]==original[key]
        for segment,oldsegment in zip(curve['segments'],original['segments'],strict=True):
            for key in ('source','edge','edge_sign'):assert segment[key]==oldsegment[key]
            assert [segment['vertices'][0],segment['vertices'][-1]]==[oldsegment['vertices'][0],oldsegment['vertices'][-1]]
            for first,v in ((True,segment['vertices'][0]),(False,segment['vertices'][-1])):
                endpoint=segment['beziers'][0 if first else -1][0 if first else -1]
                assert endpoint==shared[str(v)]
                if v in fixed_vertices:assert endpoint==[float(x)*600 for x in seed['positions'][str(v)]]
    sub=lambda x,y:(x[0]-y[0],x[1]-y[1])
    add=lambda x,y:(x[0]+y[0],x[1]+y[1])
    cross=lambda x,y:x[0]*y[1]-x[1]*y[0]
    dot=lambda x,y:x[0]*y[0]+x[1]*y[1]
    blend=lambda x,y,t:tuple(u+t*(v-u) for u,v in zip(x,y))
    orient=lambda x,y,z:cross(sub(y,x),sub(z,x))
    def minimum(c,b,a):
        values=[c,a+b+c]
        if a>0 and 0 < -b/(2*a) < 1:values.append(c-b*b/(4*a))
        return min(values)
    def moving_det(P,Q,u,v,w):
        e,f=sub(P[v],P[u]),sub(P[w],P[u]);de,df=sub(sub(Q[v],Q[u]),e),sub(sub(Q[w],Q[u]),f)
        return minimum(cross(e,f),cross(de,f)+cross(e,df),cross(de,df))
    def moving_dot(n0,n1,q0,q1):
        dn,dq=sub(n1,n0),sub(q1,q0)
        return minimum(dot(n0,q0),dot(dn,q0)+dot(n0,dq),dot(dn,dq))
    def hull(points):
        points=sorted(set(points))
        if len(points)<3:return points
        lower=[];upper=[]
        for chain,seq in ((lower,points),(upper,reversed(points))):
            for q in seq:
                while len(chain)>1 and orient(chain[-2],chain[-1],q)<=0:chain.pop()
                chain.append(q)
        return lower[:-1]+upper[:-1]
    def separate(A,B):
        for polygon in (A,B):
            for u,v in zip(polygon,polygon[1:]+polygon[:1]):
                edge=sub(v,u);n=(-edge[1],edge[0]);aa=[dot(n,q) for q in A];bb=[dot(n,q) for q in B]
                if max(aa)<min(bb) or max(bb)<min(aa):return True
        return False
    def bbox(points):
        lo=[min(q[k] for q in points) for k in range(2)];hi=[max(q[k] for q in points) for k in range(2)]
        return box(*(math.nextafter(float(x),-math.inf) for x in lo),*(math.nextafter(float(x),math.inf) for x in hi))

    # The SVG stores Float64(x*600), which differs from exact rational x*600.
    # Bind every original constrained vertex to the actual seed SVG coordinate,
    # and every constrained Steiner vertex to the SAME original-edge parameter.

    S,R,A=build_rendered_mesh(m,a);ids=m['original_vertex_ids'];index={v:i for i,v in enumerate(ids)}
    # Unmarked collar vertices shift later generated IDs. Bind by the already
    # retained exact directed geometry, never confuse the two numeric ID spaces.
    by_position={tuple(F(x)*600 for x in m['target_positions'][str(i)]):i for i in range(len(ids))}
    for curve in seed['curves']:
        for segment in curve['segments']:
            for v in segment['vertices']:
                exact=tuple(F(x)*600 for x in seed['positions'][str(v)])
                assert exact in by_position and R[by_position[exact]]==tuple(F(float(x)*600) for x in seed['positions'][str(v)])
    mins=[]
    for start,end in ((S,R),(R,A)):
        smallest=min(moving_det(start,end,*t) for t in m['triangles']);assert smallest>0
        mins.append(str(smallest))
        print('whole mesh step positive',float(smallest),round(time.monotonic()-started,2),flush=True)

    # Convex outer boundary and clockwise star-shaped holes are checked on this
    # ACTUAL refined mesh, not on unrefined source cycles.
    outer=a['outer'];holes=a['holes'];outer_triples=0
    for i,u in enumerate(outer):
        for start,end in ((S,R),(R,A)):
            # Collinear subdivisions are permitted, with nonzero edge lengths.
            assert moving_det(start,end,outer[i-1],u,outer[(i+1)%len(outer)])>=0
        outer_triples+=1
    for hole in holes:
        c=tuple(F(x)*600 for x in hole['center']);ring=hole['vertices']
        temporary=max(R)+1
        for start,end in ((S,R),(R,A)):
            start[temporary]=end[temporary]=c
            for u,v in zip(ring,ring[1:]+ring[:1]):assert moving_det(start,end,temporary,v,u)>0
            for u,v in zip(outer,outer[1:]+outer[:1]):
                for w in ring:assert moving_det(start,end,u,v,w)>0
            del start[temporary];del end[temporary]
    boxes=[]
    for h in holes:
        points=[chart[v] for chart in (S,R,A) for v in h['vertices']]
        boxes.append([(min(x[k] for x in points),max(x[k] for x in points)) for k in range(2)])
    for i,x in enumerate(boxes):
        for y in boxes[i+1:]:assert any(x[k][1]<y[k][0] or y[k][1]<x[k][0] for k in range(2))

    repl={(r['curve'],tuple(r['source']),r['piece'],r['control']):tuple(F(x) for x in r['exact_svg']) for r in replacement['replacements']}
    def split(q):
        ab=[blend(x,y,F(1,2)) for x,y in zip(q,q[1:])];bc=[blend(x,y,F(1,2)) for x,y in zip(ab,ab[1:])];z=blend(*bc,F(1,2))
        return [q[0],ab[0],bc[0],z],[z,bc[1],ab[2],q[-1]]
    def small(q):
        return max(math.hypot(*(float(x-y) for x,y in zip(q[k],blend(q[0],q[-1],F(k,3))))) for k in (1,2))<=.003
    leaves=[];joins={};actual_segments=0
    def emit(q,r,left,right,tag,depth=0):
        if small(q) and small(r):
            # A fixed positive derivative projection implies each intermediate
            # cubic is embedded and retains its source direction.
            direction=sub(q[-1],q[0])
            assert all(dot(sub(y,x),direction)>0 for curve in (q,r) for x,y in zip(curve,curve[1:]))
            leaves.append(dict(q=q,r=r,left=left,right=right,hull=None));return
        assert depth<25
        q0,q1=split(q);r0,r1=split(r);middle=tag+('middle',depth)
        emit(q0,r0,left,middle,tag+(0,),depth+1);emit(q1,r1,middle,right,tag+(1,),depth+1)
    for curve in d['curves']:
        for si,segment in enumerate(curve['segments']):
            actual_segments+=1;pieces=segment['beziers'];source=tuple(segment['source'])
            for bi,piece in enumerate(pieces):
                q=[tuple(map(F,x)) for x in piece];r=[repl.get((curve['id'],source,bi,k),x) for k,x in enumerate(q)]
                left=('vertex',segment['vertices'][0]) if bi==0 else ('piece',curve['id'],si,bi)
                right=('vertex',segment['vertices'][-1]) if bi+1==len(pieces) else ('piece',curve['id'],si,bi+1)
                for key,point0,point1 in ((left,q[0],r[0]),(right,q[-1],r[-1])):
                    if key in joins:assert joins[key]==(point0,point1),('literal shared endpoint mismatch',key)
                    joins[key]=(point0,point1)
                emit(q,r,left,right,(curve['id'],si,bi))
    print('exact directed leaves',len(leaves),'segments',actual_segments,round(time.monotonic()-started,2),flush=True)

    def geth(i):
        if leaves[i]['hull'] is None:leaves[i]['hull']=hull(leaves[i]['q']+leaves[i]['r'])
        return leaves[i]['hull']
    def joint_separate(x,y,key):
        i=0 if x['left']==key else -1;j=0 if y['left']==key else -1
        assert x['q'][i]==y['q'][j] and x['r'][i]==y['r'][j]
        xx=[sub(v,z[i]) for z in (x['q'],x['r']) for v in z if v!=z[i]]
        yy=[sub(v,z[j]) for z in (y['q'],y['r']) for v in z if v!=z[j]]
        # Search separating line through the exact moving shared endpoint. The
        # rays are only proposals; every control's sign is subsequently exact.
        ex=np.array(list(map(float,xx[-1])));ey=np.array(list(map(float,yy[-1])))
        guesses=[tuple(map(F,ex/np.linalg.norm(ex)-ey/np.linalg.norm(ey)))]
        guesses.extend(((-v[1],v[0]) for v in xx+yy))
        for n in guesses:
            vx=[dot(n,v) for v in xx];vy=[dot(n,v) for v in yy]
            if min(vx)>0 and max(vy)<0 or max(vx)<0 and min(vy)>0:return True
        return False

    tree=STRtree([bbox(x['q']+x['r']) for x in leaves]);pairs=shared=0
    for i,x in enumerate(leaves):
        for j in tree.query(bbox(x['q']+x['r'])):
            j=int(j)
            if j<=i:continue
            pairs+=1;y=leaves[j];common=set((x['left'],x['right']))&set((y['left'],y['right']))
            if common:
                assert len(common)==1 and joint_separate(x,y,next(iter(common))),('shared-join separation failed',i,j,common)
                shared+=1
            else:assert separate(geth(i),geth(j)),('swept hull separation failed',i,j)
        if i and i%3000==0:print('curve hulls',i,pairs,round(time.monotonic()-started,2),flush=True)

    port_by_vertex={ids[v]:(h,v) for h in holes for v in h['ports']}
    hole_checks=port_checks=0
    for h in holes:
        ring=h['vertices'];H=hull([chart[v] for chart in (R,A) for v in ring]);hb=bbox(H)
        center=tuple(F(x)*600 for x in h['center'])
        for i in tree.query(hb):
            i=int(i);x=leaves[i];ends=[k for k in (x['left'],x['right']) if len(k)==2 and k[0]=='vertex' and k[1] in port_by_vertex and port_by_vertex[k[1]][0] is h]
            if ends:
                assert len(ends)==1;key=ends[0];v=port_by_vertex[key[1]][1];ei=0 if x['left']==key else -1
                assert x['q'][ei]==R[v] and x['r'][ei]==A[v]
                n0,n1=sub(R[v],center),sub(A[v],center)
                for w in ring:
                    assert moving_dot(n0,n1,sub(R[v],R[w]),sub(A[v],A[w]))>=0,('moving port support fails',v,w)
                for k in range(4):
                    if k==ei%4:continue
                    assert moving_dot(n0,n1,sub(x['q'][k],R[v]),sub(x['r'][k],A[v]))>0,('port Bezier enters moving hole',i,k)
                port_checks+=1
            else:assert separate(geth(i),H),('Bezier swept hull meets moving hole',i,h['cut'],h['side'])
            hole_checks+=1

    # At the exact endpoint the analytic radial/quotient collars must fix the
    # actual BEZIER graph as well as the seed polyline. Check the same retained
    # outer halfspace polygons against every exact corrected control hull.
    collar_checks=[]
    for collar in collar_data['collars']:
        hole=next(h for h in holes if (h['cut'],h['side'])==(collar['cut'],collar['side']))
        center=tuple(F(x)*600 for x in collar['center'])
        normals=[tuple(F(x) for x in n) for n in collar['halfspace_normals']]
        offsets=[F(x)*600 for x in collar['halfspace_offsets']]
        polygon=[tuple(F(x)*600 for x in q) for q in collar['outer_positions']]
        contacts=set();tested=0
        for i in tree.query(bbox(polygon)):
            i=int(i);x=leaves[i];tested+=1
            values=[[dot(n,sub(q,center))-z for q in x['r']] for n,z in zip(normals,offsets)]
            if any(min(v)>0 for v in values):continue
            ends=[k for k in (x['left'],x['right']) if len(k)==2 and k[0]=='vertex' and k[1] in port_by_vertex and port_by_vertex[k[1]][0] is hole]
            assert len(ends)==1,('corrected Bezier intersects collar support',i,collar['cut'],collar['side'])
            key=ends[0];endpoint=0 if x['left']==key else 3
            assert any(v[endpoint]==0 and all(v[j]>0 for j in range(4) if j!=endpoint) for v in values),('collar tangent-plane support failed',i)
            contacts.add(port_by_vertex[key[1]][1])
        assert contacts==set(hole['ports'])
        collar_checks.append(dict(cut=hole['cut'],side=hole['side'],tested_leaf_hulls=tested,exact_port_contacts=len(contacts),corrected_bezier_graph_fixed=True))

    allpoints=[q for x in leaves for k in ('q','r') for q in x[k]]
    corners=[]
    if allpoints:
        lo=[min(q[k] for q in allpoints) for k in range(2)];hi=[max(q[k] for q in allpoints) for k in range(2)]
        corners=[(x,y) for x in (lo[0],hi[0]) for y in (lo[1],hi[1])]
    else:
        assert not d['curves'] and not seed['curves']
        assert not d.get('disk') and not d.get('disk_pieces') and not d.get('disk_region_marks')
    for chart in (R,A):
        for u,v in zip(outer,outer[1:]+outer[:1]):assert all(orient(chart[u],chart[v],q)>0 for q in corners)

    # D is encoded by these same directed source intervals; keep its word and all
    # crossing labels unchanged. No independent approximate D route is introduced.
    for key in ('disk_pieces','disk_region_marks','crossings','source_document','display_request'):
        assert d.get(key)==seed.get(key),('marked metadata changed',key)
    report=dict(directed_segments=actual_segments,exact_bezier_leaves=len(leaves),candidate_pairs=pairs,
                shared_join_pairs=shared,moving_hole_checks=hole_checks,exact_port_support_checks=port_checks,
                mesh_step_minimum_double_areas_svg=mins,outer_vertices=len(outer),
                final_graph_free_collar_checks=collar_checks,
                literal_normalized_beziers_to_exact_anchor_isotopy=True,
                source_to_CDT_overlay_complete=False,full_native_marked_surface_complete=False,
                elapsed_seconds=time.monotonic()-started)

    return report


def verify_svg_coordinates(d):
    """Read back actual SVG numbers under the declared binary64 decoder.

    Source: draw_marked_cut.render's source-tagged curve/D paths. This checks
    serialized controls, implicit starts, shared vertices, circle attributes
    and retained polygonal cut chains. Native provenance and attachment-chart
    validity are separate mandatory checks; no pixel/decimal-real claim is made.
    """
    import hashlib,re,xml.etree.ElementTree as ET
    from draw_marked_cut import disk_boundary_cycles,attachment_boundary_path
    assert d.get('svg_coordinate_encoding')=='float64-roundtrip'
    svg=d['svg'];root=ET.fromstring(svg);ns={'s':'http://www.w3.org/2000/svg'}
    assert root.tag=='{http://www.w3.org/2000/svg}svg'
    number=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
    token=r'[MCLHVZmhvz]|'+number
    def commands(path):
        assert not re.sub(r'[\s,]','',re.sub(token,'',path)), 'Unsupported SVG path token'
        tokens=re.findall(token,path);out=[];i=0
        counts={'M':2,'C':6,'L':2,'H':1,'V':1,'Z':0,'m':2,'h':1,'v':1,'z':0}
        while i<len(tokens):
            command=tokens[i];i+=1;count=counts[command]
            values=list(map(float,tokens[i:i+count]));assert len(values)==count
            out.append((command,values));i+=count
        return out
    xy=lambda q:[float(q[0]),-float(q[1])]
    cubic=lambda q:[x for point in q[1:] for x in xy(point)]
    groups=[g for g in root.findall('s:g',ns) if 'data-curve' in g.attrib and not g.get('data-orientation')]
    assert len(groups)==len(d['curves']);endpoints={};pieces=0
    for group,curve in zip(groups,d['curves'],strict=True):
        assert group.get('data-curve')==curve['id']
        paths=group.findall('s:path',ns);assert len(paths)==len(curve['segments'])
        for node,segment in zip(paths,curve['segments'],strict=True):
            decoded=commands(node.get('d'));assert decoded[0]==('M',xy(segment['beziers'][0][0]))
            assert len(decoded)==len(segment['beziers'])+1
            point=decoded[0][1]
            for command,q in zip(decoded[1:],segment['beziers'],strict=True):
                assert point==xy(q[0]) and command==('C',cubic(q))
                point=command[1][-2:];pieces+=1
            for v,q in ((segment['vertices'][0],decoded[0][1]),(segment['vertices'][-1],point)):
                if v in endpoints:assert endpoints[v]==q
                endpoints[v]=q
                assert q==xy(d['source_endpoint_normalization']['shared_svg_endpoints'][str(v)])
    holes=root.findall('s:circle[@data-puncture]',ns);assert len(holes)==len(d['boundaries'])
    for node,b in zip(holes,d['boundaries'],strict=True):
        x,y=d['positions'][str(b['cap'])]
        assert [float(node.get(k)) for k in ('cx','cy','r')]==[600*x,-600*y,600*b['radius']]
    crossings=root.findall('s:circle[@data-crossing]',ns);assert len(crossings)==len(d['crossings'])
    for node,crossing in zip(crossings,d['crossings'],strict=True):
        assert node.get('data-crossing')==crossing['id']
        assert [float(node.get(k)) for k in ('cx','cy')]==endpoints[crossing['vertex']]
    fills=root.findall('s:path[@data-disk]',ns)
    if d.get('disk'):
        assert len(fills)==1
        expected=[('M',xy(d['disk_beziers'][0][0]))]+[('C',cubic(q)) for q in d['disk_beziers']]+[('Z',[])]
        assert commands(fills[0].get('d'))==expected
    else:
        assert len(fills)==len(d.get('disk_pieces',[]))
        segments={(c['id'],tuple(s['source'])):s for c in d['curves'] for s in c['segments']}
        for node,piece in zip(fills,d.get('disk_pieces',[]),strict=True):
            assert node.get('data-disk')==str(piece['mark']) and node.get('fill-rule')=='evenodd'
            expected=[]
            if piece['contains_infinity']:
                left,top,width,height=map(float,root.get('viewBox').split())
                expected=[('M',[left,top]),('h',[width]),('v',[height]),('h',[-width]),('Z',[])]
            for parts in disk_boundary_cycles(d,piece):
                first=parts[0];segment=segments[first['curve'],tuple(first['source'])]
                v=segment['vertices'][0 if first['direction']==1 else -1]
                start=endpoints[v];expected.append(('M',start));point=start
                for part in parts:
                    if part['kind']=='curve':
                        segment=segments[part['curve'],tuple(part['source'])]
                        controls=segment['beziers'] if part['direction']==1 else [q[::-1] for q in segment['beziers'][::-1]]
                        for q in controls:
                            assert point==xy(q[0]);expected.append(('C',cubic(q)));point=xy(q[-1])
                    else:
                        chain=attachment_boundary_path(d,part);assert point==xy(chain[0])
                        for q in chain[1:]:expected.append(('L',xy(q)))
                        point=xy(chain[-1])
                assert point==start;expected.append(('Z',[]))
            assert commands(node.get('d'))==expected
    return dict(numeric_decoder='IEEE binary64 round-trip parse, then exact binary rational interpretation',
        cubic_pieces=pieces,control_points=4*pieces,shared_endpoint_vertices=len(endpoints),
        puncture_circle_primitives=len(holes),crossings=len(crossings),D_fill_paths=len(fills),
        all_emitted_coordinates_recovered=True,svg_sha256=hashlib.sha256(svg.encode()).hexdigest(),
        literal_circle_glyph_incidence_claim=False,hardware_rasterization_claim=False)
