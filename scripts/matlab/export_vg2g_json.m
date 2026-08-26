function export_vg2g_json(inputPath,outputPath)
% export_vg2g_json Run only the MATLAB vg2g stage and write JSON results.
%
% This is the Sage-independent reference path used by the Julia parity
% harness. It deliberately calls VirtualLink.debugVG2G rather than calcVG2G,
% because calcVG2G continues into setSageLink.
    arguments
        inputPath (1,1) string
        outputPath (1,1) string
    end

    inputs=jsondecode(fileread(inputPath));
    results=repmat(struct("name","","gauss",{{}},"orientation",[]),size(inputs));
    for index=1:numel(inputs)
        rgc=toGaussCell(inputs(index).rgauss);
        orientation=double(inputs(index).orientation(:).');
        [gauss,planarOrientation]=VirtualLink.debugVG2G(rgc,orientation);
        results(index).name=string(inputs(index).name);
        results(index).gauss=gauss;
        results(index).orientation=planarOrientation;
    end

    payload=struct("schema_version",1,"matlab_version",version,"results",results);
    fileID=fopen(outputPath,"w");
    assert(fileID~=-1,"export_vg2g_json:OpenFailed", ...
        "Cannot open output file: %s",outputPath)
    cleanup=onCleanup(@()fclose(fileID));
    fwrite(fileID,jsonencode(payload,PrettyPrint=true),"char");
end

function gauss=toGaussCell(value)
% JSON nested arrays decode as a numeric matrix when rectangular and as a
% cell array when ragged. Normalize both forms to VirtualLink's cell input.
    if iscell(value)
        gauss=cellfun(@(component)double(component(:).'),value,UniformOutput=false);
    elseif isvector(value)
        gauss={double(value(:).')};
    else
        gauss=arrayfun(@(row){double(value(row,:))},1:size(value,1));
    end
end
