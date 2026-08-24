function topology_replay_record(callId, lineNumber, phase)
%TOPOLOGY_REPLAY_RECORD Save topology object state at a plot boundary.

outputDirectory = string(getappdata(groot,"TopologyReplayOutputDirectory"));
variables = evalin("caller","whos");
state = struct("callId",callId,"lineNumber",lineNumber,"phase",phase, ...
    "objects",struct());
for variable = variables(:).'
    if ~ismember(variable.class,["VirtualLink","URDiagram"])
        continue
    end
    value = evalin("caller",variable.name);
    state.objects.(variable.name) = serializeTopologyObject(value);
end
save(fullfile(outputDirectory,sprintf("call_%03d_%s.mat",callId,phase)),"state");
end

function serialized = serializeTopologyObject(value)
serialized = repmat(struct("class","","size",[],"properties",struct()),size(value));
for objectIndex = 1:numel(value)
    object = value(objectIndex);
    serialized(objectIndex).class = class(object);
    serialized(objectIndex).size = size(object);
    names = properties(object);
    for propertyIndex = 1:numel(names)
        name = names{propertyIndex};
        try
            serialized(objectIndex).properties.(name) = object.(name);
        catch exception
            serialized(objectIndex).properties.(name) = ...
                "<unavailable: "+string(exception.identifier)+">";
        end
    end
end
end
