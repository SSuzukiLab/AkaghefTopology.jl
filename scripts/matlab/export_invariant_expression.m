function expression = export_invariant_expression(link)
%EXPORT_INVARIANT_EXPRESSION Reproduce VirtualLink.calcInvariant's tensor expression.
isWeighted = link.isWeighted;
link.calcREdgeTable;
[~,orientation] = link.getRGaussCode();
edgeTable = link.REdgeTable;
antipodeFactors = repmat("(S^%d){%d,%d}",1,height(edgeTable));
indices = abs(edgeTable.Crossing)*4+(edgeTable.Crossing>0)-[1,3];
if isWeighted
    powers = 2*edgeTable.Weight;
else
    powers = zeros(size(edgeTable.Weight));
end
zeroPower = powers==0;
for index=1:height(edgeTable)
    antipodeFactors(index)=sprintf(antipodeFactors(index),powers(index), ...
        indices(index,1),indices(index,2));
end
antipodeFactors=antipodeFactors(~zeroPower);
indexMap=dictionary( ...
    [indices(zeroPower,1);indices(~zeroPower,1);indices(:,2)], ...
    [indices(zeroPower,2);indices(~zeroPower,1);indices(:,2)]);
vertexFactors=repmat("Tp",1,length(orientation));
vertexFactors(orientation<0)="Tm";
for index=1:length(orientation)
    vertexFactors(index)=vertexFactors(index)+sprintf("{%s}", ...
        join(string(indexMap(4*(index-1)+(1:4))),","));
end
expression=join([vertexFactors,antipodeFactors]);
end
