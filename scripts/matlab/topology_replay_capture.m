function topology_replay_capture
%TOPOLOGY_REPLAY_CAPTURE Export every not-yet-captured replay figure.

outputDirectory = string(getappdata(groot,"TopologyReplayOutputDirectory"));
figureCount = getappdata(groot,"TopologyReplayFigureCount");
figures = findall(groot,"Type","figure");
figures = figures(isgraphics(figures));
if ~isempty(figures)
    numbers = arrayfun(@(figureHandle)figureHandle.Number,figures);
    [~,order] = sort(numbers);
    figures = figures(order);
end
for figureHandle = figures(:).'
    if isappdata(figureHandle,"TopologyReplayCaptured")
        continue
    end
    axesHandles = findall(figureHandle,"Type","axes");
    if isempty(axesHandles)
        continue
    end
    figureCount = figureCount+1;
    try
        exportgraphics(axesHandles(end), ...
            fullfile(outputDirectory,sprintf("figure_%03d.png",figureCount)), ...
            "Resolution",144);
        setappdata(figureHandle,"TopologyReplayCaptured",true);
    catch exception
        fprintf(2,"EXPORT_ERROR figure=%d id=%s message=%s\n", ...
            figureCount,exception.identifier,exception.message);
    end
end
setappdata(groot,"TopologyReplayFigureCount",figureCount);
end
