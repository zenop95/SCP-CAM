function [dynMaps, STM, cartMaps, state, cart, ...
          nu, pp] = propagateDA(iter, oldTraj, validation, pp)
% Per-primary DA propagation, with common control history.
N  = pp.N;
Np = numel(pp.primary);
t  = pp.t;

if validation
    order = 1;
    u = reshape(pp.u, 3, N)';
else
    order = 2;
    u = zeros(N, 3);
end

% Each output has an added dimension for primary's mixand
dynMaps   = zeros(6, 9, N, Np);
STM       = zeros(6, 6, N, Np);
cartMaps  = zeros(6, 6, N, Np);
state     = zeros(6, N, Np);
cart      = zeros(6, N, Np);
nu=[];

for ip = 1:Np

    pp = aidaInit(pp, 'primary', ip);
    input=struct('N',N,'t',t,'et',pp.et,'order',order,...
        'uMax',pp.uMax,'scaling',pp.scaling(1),...
        'flagRtn',pp.flagRtn,'u',U);

    if iter == 1 || validation
        input.flagProp = 1;
        input.x0 = pp.primary(ip).x0(:)';

    else
        input.flagProp = 2;

        % oldTraj has one trajectory per primary mixand:
        % size(oldTraj) = [6,N,Np].
        if size(oldTraj,3) ~= Np
            error('propagateDA:OldTrajSize', ...
                ['For iter > 1, oldTraj must be 6 x N x Np. ', ...
                 'Received size %s; Np = %d.'], ...
                 mat2str(size(oldTraj)), Np);
        end

        input.oldTraj = reshape(oldTraj(:,:,ip), 6, [])';

        % Same optimized control history for every mixand.
        input.u = reshape(pp.majorIter(end).u, 3, [])';
        % u=reshape(pp.majorIter(end).u,3,N)';
    end

    % Write input for this mixand.
    fid = fopen('./data_sharing/initial_state.json', 'w');

    if fid == -1
        error('propagateDA:FileOpen', ...
              'Could not open initial_state.json.');
    end

    fwrite(fid, jsonencode(input), 'char');
    fclose(fid);

    % Run C++ propagation synchronously.
    [status, cmdOut] = system('wsl ./build/bin/propAida');

    if status ~= 0
        error('propagateDA:CppFailure', ...
              'C++ propagation failed for primary mixand %d failed: %s', ip, cmdOut);
    end

    output = jsondecode(fileread('./data_sharing/out_prop.json'));

    % Check C++ output dimensions.
    state_ip = output.state';
    cart_ip  = state_ip;
   

    if ~isequal(size(state_ip), [6,N])
        error('propagateDA:StateSize', ...
            'Expected state size 6 x %d; received %s.', ...
            N, mat2str(size(state_ip)));
    end

    state(:,:,ip) = state_ip;
    cart(:,:,ip)  = cart_ip;

    % Dynamics maps.
    dm = permute(output.dynMaps, [2,3,1]);
    dynMaps(:,:,:,ip) = dm;
    STM(:,:,:,ip)     = dm(1:6,1:6,:);
   
    % Cartesian map is identity in the current implementation.
    cartMaps(:,:,:,ip) = repmat(eye(6), [1,1,N]); % cartMaps(:,:,:,a)=repmat(eye(6),1,1,N);

    if isfield(output, 'trustRegion')
        nu_ip = output.trustRegion';

        if isempty(nu)
            nu = zeros(size(nu_ip,1), size(nu_ip,2), Np);
        end

        nu(:,:,ip) = nu_ip;
    end
    
    if isfield(output,'timeMs')
        if ~isfield(pp,'timeSubtr')
            pp.timeSubtr=0;
        end
        pp.timeSubtr = pp.timeSubtr + output.timeMs/1000;
    end
end
end
