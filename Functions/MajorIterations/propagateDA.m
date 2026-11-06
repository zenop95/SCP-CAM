% function [dynMaps,STM,llMaps,cartMaps,state,cart,lla_const,nu,pp] = ...
%                                   propagateDA(iter,oldTraj,validation,pp)
% % PropagateDA calls the relevant C++ function to perform the first or
% % second order DA propagation that gives the constant part of the
% % propagation, the linear maps and the nonlinearity parameter nu.
% % 
% % INPUT:  t          = [-] (N,1) Time vector
% %         iter       = [-] (1,1) Current iteration counter
% %         oldTraj    = [-] (6N,1) Relative trajectory from previous iteration
% %         validation = [bool] (1,1) Flag to perform validation propagation
% %         pp         = postprocess structure
% %
% % OUTPUT: dynMaps  = [-] (6,9,N) Linear maps for the dynamics
% %         STM      = [-] (6,6,N) State transition matrices
% %         llMaps   = [-] (2,3,N) Linear maps from state to geodetic   
% %         cartMaps = [-] (6,6,N) Linear maps from state to cartesian   
% %         state    = [-] (6,N)   Constant part of the propagated state
% %         cart     = [-] (6,N)   Constant part of the cartesian coordinates
% %         ll_const = [-] (2,N)   Constant part of the geodetic coordinates
% %         nu       = [-] (6,N)   Nonlinearity parameter nu
% %         pp       = [struct]    Postprocess structure. 
% % 
% % Author: Zeno Pavanello, 2022
% % E-mail: zpav176@aucklanduni.ac.nz
% %--------------------------------------------------------------------------
% N=pp.N; t=pp.t;
% if validation 
%     order = 1; 
%     u     = reshape(pp.u,3*N,1); 
% else 
%     order = 2; 
%     u     = zeros(3*N,1); 
% end
% 
% input         = struct();
% input.N       = N; 
% input.t       = t; 
% input.et      = pp.et;
% input.order   = order; 
% input.uMax    = pp.uMax;
% input.scaling = pp.scaling(1);
% input.flagRtn = pp.flagRtn;
% 
% if iter==1||validation
%     input.flagProp = 1;
%     input.x0       = pp.primary.x0(:)';
%     input.u        = reshape(u,3,[])';
% else
%     input.flagProp = 2;
%     input.oldTraj  = reshape(oldTraj,6,[])';
%     input.u        = reshape(pp.majorIter(end).u,3,[])';
% end
% 
% fid = fopen('./data_sharing/initial_state.json','w'); 
% fwrite(fid,jsonencode(input),'char'); 
% fclose(fid);
% !wsl ./build/bin/propAida
% output = jsondecode(fileread('./data_sharing/out_prop.json'));
% 
% state     = output.state';
% cart      = state;
% lla_const = output.latlonConst';
% 
% % build dynMaps (9x6xN), cartMaps(6x6), llMaps(9x3)
% dynMaps  = permute(output.dynMaps,[2,3,1]);          % 9x6xN
% llMaps   = permute(output.latlonMaps,[2,3,1]);       % 9x3xN
% STM      = dynMaps(1:6,1:6,:);
% llMaps   = llMaps(:,1:6,:);
% cartMaps = repmat(eye(6),[1,1,N]);
% if isfield(output,'trustRegion')
%     nu = output.trustRegion'; 
% else 
%     nu = []; 
% end
% 
% pp.timeSubtr = pp.timeSubtr+output.timeMs/1000;
% end

% function [dynMaps, STM, llMaps, cartMaps, state, cart, lla_const, nu, pp] = propagateDA(iter, oldTraj, validation, pp)
% 
% N  = pp.N; % discretization (forward + backward)
% t  = pp.t;
% Np = numel(pp.primary);
% 
% if validation
%     order = 1;
%     u = reshape(pp.u, 3*N, 1);
% else
%     order = 2;
%     u = zeros(3*N, 1);
% end
% 
% % Allocate outputs with a final dimension for the primary mixand.
% dynMaps   = zeros(9, 6, N, Np);
% STM       = zeros(6, 6, N, Np);
% llMaps    = zeros(9, 6, N, Np);
% cartMaps  = zeros(6, 6, N, Np);
% state     = zeros(6, N, Np);
% cart      = zeros(6, N, Np);
% lla_const = zeros(2, N, Np);
% nu        = [];
% 
% for ip = 1:Np
% 
%     input       = struct();
%     input.N     = N;
%     input.t     = t;
%     input.et    = pp.et;
%     input.order = order;
%     input.uMax  = pp.uMax;
%     input.scaling = pp.scaling(1);
%     input.flagRtn = pp.flagRtn;
% 
%     if iter == 1 || validation
%         input.flagProp = 1;
% 
%         % One initial state for this particular primary mixand.
%         input.x0 = pp.primary(ip).x0(:)';
% 
%         % Identical control history for every primary mixand.
%         input.u = reshape(u, 3, [])';
% 
%     else
%         input.flagProp = 2;
% 
%         % Each mixand must have its OWN previous trajectory.
%         % oldTraj must be 6 x N x Np.
%         input.oldTraj = reshape(oldTraj(:,:,ip), 6, [])';
% 
%         % All primary mixands use the same optimized control.
%         input.u = reshape(pp.majorIter(end).u, 3, [])';
%     end
% 
%     % Write input and run the existing C++ propagator.
%     fid = fopen('./data_sharing/initial_state.json', 'w');
%     if fid == -1
%         error('propagateDA:FileOpen', ...
%               'Could not open initial_state.json for writing.');
%     end
% 
%     fwrite(fid, jsonencode(input), 'char');
%     fclose(fid);
% 
%     status = system('wsl ./build/bin/propAida');
%     if status ~= 0
%         error('propagateDA:CppFailure', ...
%               'C++ propagation failed for primary mixand %d.', ip);
%     end
% 
%     output = jsondecode(fileread('./data_sharing/out_prop.json'));
% 
%     % Store this mixand's propagation.
%     state(:,:,ip) = output.state';
%     cart(:,:,ip)  = state(:,:,ip);
%     lla_const(:,:,ip) = output.latlonConst';
% 
%     dm = permute(output.dynMaps, [2, 3, 1]);
%     lm = permute(output.latlonMaps, [2, 3, 1]);
% 
%     dynMaps(:,:,:,ip) = dm;
%     STM(:,:,:,ip)     = dm(1:6,1:6,:);
%     llMaps(:,:,:,ip)  = lm(:,1:6,:);
% 
%     cartMaps(:,:,:,ip) = repmat(eye(6), [1,1,N]);
% 
%     if isfield(output, 'trustRegion')
%         nu_ip = output.trustRegion';
%         if isempty(nu)
%             nu = zeros(size(nu_ip,1), size(nu_ip,2), Np);
%         end
%         nu(:,:,ip) = nu_ip;
%     end
% 
%     pp.timeSubtr = pp.timeSubtr + output.timeMs/1000;
% end
% 
% end

function [dynMaps, STM, llMaps, cartMaps, state, cart, ...
          lla_const, nu, pp] = propagateDA(iter, oldTraj, validation, pp)

N  = pp.N;
t  = pp.t;
Np = numel(pp.primary);

if validation
    order = 1;
    u = reshape(pp.u, 3*N, 1);
else
    order = 2;
    u = zeros(3*N, 1);
end

% Each output has a final dimension for the primary mixand.
dynMaps   = zeros(6, 9, N, Np);
STM       = zeros(6, 6, N, Np);
llMaps    = zeros(6, 9, N, Np);
cartMaps  = zeros(6, 6, N, Np);
state     = zeros(6, N, Np);
cart      = zeros(6, N, Np);

% C++ latlonConst contains THREE quantities.
lla_const = zeros(3, N, Np);

nu = [];

for ip = 1:Np

    input = struct();

    input.N       = N;
    input.t       = t;
    input.et      = pp.et;
    input.order   = order;
    input.uMax    = pp.uMax;
    input.scaling = pp.scaling(1);
    input.flagRtn = pp.flagRtn;

    if iter == 1 || validation

        input.flagProp = 1;
        input.x0 = pp.primary(ip).x0(:)';

        % Common control for all primary mixands.
        input.u = reshape(u, 3, [])';

    else

        input.flagProp = 2;

        % oldTraj must contain one trajectory per primary mixand:
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
    status = system('wsl ./build/bin/propAida');

    if status ~= 0
        error('propagateDA:CppFailure', ...
              'C++ propagation failed for primary mixand %d.', ip);
    end

    output = jsondecode(fileread('./data_sharing/out_prop.json'));

    % Check C++ output dimensions.
    state_ip = output.state';
    cart_ip  = state_ip;
    lla_ip   = output.latlonConst';

    if ~isequal(size(state_ip), [6,N])
        error('propagateDA:StateSize', ...
            'Expected state size 6 x %d; received %s.', ...
            N, mat2str(size(state_ip)));
    end

    if ~isequal(size(lla_ip), [3,N])
        error('propagateDA:LLASize', ...
            'Expected latlonConst size 3 x %d; received %s.', ...
            N, mat2str(size(lla_ip)));
    end

    state(:,:,ip) = state_ip;
    cart(:,:,ip)  = cart_ip;
    lla_const(:,:,ip) = lla_ip;

    % Dynamics maps.
    dm = permute(output.dynMaps, [2,3,1]);
    lm = permute(output.latlonMaps, [2,3,1]);

    dynMaps(:,:,:,ip) = dm;
    STM(:,:,:,ip)     = dm(1:6,1:6,:);
    %llMaps(:,:,:,ip)  = lm(:,1:6,:);

    % Cartesian map is identity in the current implementation.
    cartMaps(:,:,:,ip) = repmat(eye(6), [1,1,N]);

    if isfield(output, 'trustRegion')
        nu_ip = output.trustRegion';

        if isempty(nu)
            nu = zeros(size(nu_ip,1), size(nu_ip,2), Np);
        end

        nu(:,:,ip) = nu_ip;
    end

    pp.timeSubtr = pp.timeSubtr + output.timeMs/1000;
end

end
