function [messageOut] = buildReactionIDlayout_diff(mapReactionIdsA, mapReactionIdsB, filename)
% Builds a layout for MINERVA from two sets of reaction IDs.
%
% USAGE:
%
%    [messageOut] = buildReactionIDlayout_diff(mapReactionIdsA, mapReactionIdsA, filename)
%
% INPUTS:
%    mapReactionIdsA:       cell vector of reaction IDs (e.g. recon22.rxns(1:200))
%    mapReactionIdsB:       cell vector of reaction IDs (e.g. recon22.rxns(101:300))
%    filename:              name for your layout text file (e.g. 'DCM')
% 
% OUTPUT:
%    messageOut:            message!
%
% .. Author: Michiel Adriaens (17-01-2020)

% Create overlay content reaction by reaction
content = 'name\treactionIdentifier\tlineWidth\tcolor\n';
allMapReactionIds = unique([mapReactionIdsA; mapReactionIdsB]);
for i=1:length(allMapReactionIds)
    mapReactionId = allMapReactionIds{i};
    presentInA = any(strcmp(mapReactionIdsA, mapReactionId) == 1);
    presentInB = any(strcmp(mapReactionIdsB, mapReactionId) == 1);
    % Present in both
    if presentInA && presentInB
        line = strcat('\t', mapReactionId, '\t', 5, '\t', '#7c6c5a', '\n'); % taupe
        content = strcat(content, line);
    end
    
    % Present in A only
    if presentInA && ~presentInB
        line = strcat('\t', mapReactionId, '\t', 5, '\t', '#FFD700', '\n'); % gold
        content = strcat(content, line);
    end
    
    % Present in B only
    if ~presentInA && presentInB
        line = strcat('\t', mapReactionId, '\t', 5, '\t', '#1E90FF', '\n'); % blue
        content = strcat(content, line);
    end

end

% Export overlay to upload manually:
fileID = fopen(strcat(filename, '_reactionIDs_Difference_overlay.txt'), 'w');
fprintf(fileID, content);
fclose(fileID);
messageOut = "done";