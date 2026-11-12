function pp = aidaInit(pp, sc, j)
% aidaInit  Write AIDA initialization parameters to JSON.
%
%   pp = aidaInit(pp, sc, j)
%     sc ∈ {'primary','secondary'}
%     j  is required when sc == 'secondary' (index of the secondary)
%
% Output JSON path:
%   ./data_sharing/AIDA_init.json
%
% JSON schema (as read by read_aida_params_json in C++):
% {
%   "flag1": <int>,
%   "flag2": <int>,
%   "flag3": <int>,
%   "gravOrd": <int>,
%   "mass": <double>,
%   "A_drag": <double>,
%   "Cd": <double>,
%   "A_srp": <double>,
%   "Cr": <double>
% }% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
    % Ensure data_sharing exists
    outDir = "./data_sharing";
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % Validate pp.aida fields exist
    requiredAida = ["flag1","flag2","flag3","gravOrd"];
    for f = requiredAida
        if ~isfield(pp,'aida') || ~isfield(pp.aida, f)
            error('aidaInit:MissingAidaField', ...
                  'pp.aida.%s is missing.', f);
        end
    end

    % Select spacecraft parameters
    switch lower(sc)
        case "primary"
            if ~isfield(pp,'primary') || j > numel(pp.primary)
                error('aidaInit:BadIndex',...
                      'Primary index j=%d is invalid.', j);
            end
            S = pp.primary(j);
        case "secondary"
            if ~isfield(pp,'secondary') || j > numel(pp.secondary)
                error('aidaInit:BadIndex',...
                      'Secondary index j=%d is invalid.', j);
            end
            S = pp.secondary(j);
        otherwise
            error('aidaInit:BadSC','sc must be ''primary'' or ''secondary''.');
    end

    % Validate required spacecraft fields
    requiredSC = "mass"; %["mass","A_drag","Cd","A_srp","Cr"];
    for f = requiredSC
        if ~isfield(S, f)
            error('aidaInit:MissingSCField', ...
                  'Field %s is missing for %s.', f, sc);
        end
    end

    % Build JSON payload matching the C++ parser (read_aida_params_json)
    payload = struct( ...
        'flag1',   int64(pp.aida.flag1), ...
        'flag2',   int64(pp.aida.flag2), ...
        'flag3',   int64(pp.aida.flag3), ...
        'gravOrd', int64(pp.aida.gravOrd), ...
        'mass',    S.mass); % 'A_drag',  S.A_drag,  'Cd',      S.Cd, 'A_srp',   S.A_srp,  'Cr',      S.Cr
    

    % Write pretty JSON (easy to inspect)
    jsonText = jsonencode(payload, 'PrettyPrint', true);

    % Save
    outPath = fullfile(outDir, 'AIDA_init.json');
    fid = fopen(outPath,'w');
    if fid < 0
        error('aidaInit:IO','Cannot open %s for writing.', outPath);
    end
    fwrite(fid, jsonText, 'char');
    fclose(fid);
end