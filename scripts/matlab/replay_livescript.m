function result = replay_livescript(sourceFile, outputDirectory, topologyRoot)
%REPLAY_LIVESCRIPT Replay every topology plot in an extracted Live Script.
%   The sourceFile must be the plain .m produced by MATLAB's Live Script
%   converter. Plot calls are isolated so that one failing diagram does not
%   silently suppress later diagrams in the same Live Script.

arguments
    sourceFile (1,1) string {mustBeFile}
    outputDirectory (1,1) string
    topologyRoot (1,1) string {mustBeFolder}
end

addpath(genpath(topologyRoot));
addpath(fileparts(mfilename("fullpath")));
outputDirectory = string(java.io.File(char(outputDirectory)).getCanonicalPath());
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
close all force;

source = readlines(sourceFile);
plotLine = ~startsWith(strtrim(source), "%") & ...
    ~cellfun(@isempty, regexp(cellstr(source), ...
    '\.(plot|plotSage)([^A-Za-z0-9_]|$)', "once"));
closeLine = ~startsWith(strtrim(source), "%") & ...
    ~cellfun(@isempty, regexp(cellstr(source),'^\s*close\s+all',"once"));
callIds = zeros(size(source));
callIds(plotLine) = 1:nnz(plotLine);

instrumented = strings(0,1);
for lineNumber = 1:numel(source)
    line = source(lineNumber);
    if plotLine(lineNumber) || closeLine(lineNumber)
        instrumented(end+1) = "topology_replay_capture;"; %#ok<AGROW>
    end
    if ~plotLine(lineNumber)
        instrumented(end+1) = line; %#ok<AGROW>
        continue
    end
    callId = callIds(lineNumber);
    instrumented(end+1) = sprintf("topology_replay_record(%d,%d,'before');",callId,lineNumber); %#ok<AGROW>
    instrumented(end+1) = "try"; %#ok<AGROW>
    instrumented(end+1) = line; %#ok<AGROW>
    instrumented(end+1) = sprintf("topology_replay_record(%d,%d,'after');",callId,lineNumber); %#ok<AGROW>
    instrumented(end+1) = "catch topologyReplayException"; %#ok<AGROW>
    instrumented(end+1) = sprintf("topology_replay_error(%d,%d,topologyReplayException);",callId,lineNumber); %#ok<AGROW>
    instrumented(end+1) = "end"; %#ok<AGROW>
end

setappdata(groot,"TopologyReplayOutputDirectory",char(outputDirectory));
setappdata(groot,"TopologyReplayErrors",struct( ...
    "callId",{},"lineNumber",{},"identifier",{},"message",{}));
setappdata(groot,"TopologyReplayFigureCount",0);

temporaryDirectory = string(tempname);
mkdir(temporaryDirectory);
cleanup = onCleanup(@()rmdir(temporaryDirectory,"s"));
[~,stem] = fileparts(sourceFile);
instrumentedFile = fullfile(temporaryDirectory,stem+".m");
writelines(instrumented,instrumentedFile);

scriptError = [];
try
    run(instrumentedFile);
catch exception
    scriptError = exception;
end
topology_replay_capture;

errors = getappdata(groot,"TopologyReplayErrors");
exportedFigureCount = getappdata(groot,"TopologyReplayFigureCount");
result = struct( ...
    "source",sourceFile, ...
    "plotCalls",nnz(plotLine), ...
    "figures",exportedFigureCount, ...
    "plotErrors",errors, ...
    "scriptError",scriptError);
save(fullfile(outputDirectory,"replay_result.mat"),"result");

summary = [ ...
    "source="+sourceFile
    "plot_calls="+nnz(plotLine)
    "figures="+exportedFigureCount
    "plot_errors="+numel(errors)
    "script_error="+string(~isempty(scriptError))];
writelines(summary,fullfile(outputDirectory,"replay_summary.txt"));

close all force;
clear cleanup;
end
