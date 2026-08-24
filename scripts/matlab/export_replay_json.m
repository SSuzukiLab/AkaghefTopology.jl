function export_replay_json(replayRoot,outputFile)
%EXPORT_REPLAY_JSON Export stable topology fields from replay MAT snapshots.
arguments
    replayRoot (1,1) string
    outputFile (1,1) string
end

files = dir(fullfile(replayRoot,"**","call_*_after.mat"));
records = repmat(struct("relativePath","","callId",0,"lineNumber",0, ...
    "objects",struct()),numel(files),1);
for fileIndex = 1:numel(files)
    inputFile = fullfile(files(fileIndex).folder,files(fileIndex).name);
    loaded = load(inputFile,"state");
    records(fileIndex).relativePath = erase(string(inputFile),replayRoot+filesep);
    records(fileIndex).callId = loaded.state.callId;
    records(fileIndex).lineNumber = loaded.state.lineNumber;
    objectNames = fieldnames(loaded.state.objects);
    for nameIndex = 1:numel(objectNames)
        name = objectNames{nameIndex};
        source = loaded.state.objects.(name);
        target = repmat(struct("class","","size",[],"properties",struct()),size(source));
        for objectIndex = 1:numel(source)
            target(objectIndex).class = source(objectIndex).class;
            target(objectIndex).size = source(objectIndex).size;
            target(objectIndex).properties = selectedProperties(source(objectIndex).properties);
        end
        records(fileIndex).objects.(name) = target;
    end
end

parent = fileparts(outputFile);
if ~isfolder(parent), mkdir(parent); end
writelines(jsonencode(records,PrettyPrint=true),outputFile);
end

function target = selectedProperties(source)
target = struct();
names = ["GaussCode","RGaussCode","orientation","ROrientation", ...
    "PDCode","RPDCode","Ncircle","isWeighted","weightRE","isCut"];
for name = names
    if isfield(source,name)
        target.(name) = source.(name);
    end
end
for name = ["EdgeTable","REdgeTable"]
    if ~isfield(source,name) || ~istable(source.(name))
        continue
    end
    tableValue = source.(name);
    edge = struct();
    % `Position` is the Sage/allocate_pos polyline used by VirtualLink.plot.
    % It is required for geometric Julia-vs-MATLAB verification; exporting
    % only Crossing/Arc can prove combinatorics but cannot detect a different
    % valid orthogonal embedding.
    for column = ["ID","Crossing","Weight","StrandID","Arc","Position"]
        if ismember(column,string(tableValue.Properties.VariableNames))
            if column == "Position"
                % JSON has no complex scalar.  Preserve each MATLAB complex
                % polyline as N-by-2 [real imag] rows in REdgeTable order.
                edge.Position = positionRows(tableValue.Position);
            else
                edge.(column) = tableValue.(column);
            end
        end
    end
    target.(name) = edge;
end
end

function rows = positionRows(value)
rows = cell(size(value));
for index = 1:numel(value)
    point = value{index};
    if isempty(point)
        rows{index} = zeros(0,2);
    else
        point = point(:);
        rows{index} = [real(point),imag(point)];
    end
end
end
