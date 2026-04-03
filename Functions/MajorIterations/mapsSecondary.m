function [x, C] = mapsSecondary(Ntca,x0,C0,pp)
% mapsSecondary Propagates the secondary orbit for the required time 
% and it propagates the s/c covariance with the associatedd linear maps. 
% 
% INPUT:  N     = [-] (1,1) Number of propagation nodes
%         t     = [-] (*,1) Vector of times
%         et    = [-] (1,1) Ephemeris time at conjunction
%         x0    = [-] (6,1) State of the secondary spacecraft at conjunction
%         C0    = [-] (6,6) Covariance of the secondary spacecraft at conjunction
%         pp    = [struct]  Postprocess structure
%
% OUTPUT: x  = [-] (6,N)   Propagated state
%         C  = [-] (6,6,N) Propagated covariance
% 
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------
%% Fix forward and backward times

Dt=(Ntca-1)*pp.dt;
et=pp.et+Dt;
t=pp.t-Dt;
N=pp.N;

input=struct();
input.N=N;
input.scaling=pp.scaling(1);
input.t=t(:)';
input.et=et;
input.x0=x0(:)';

fid=fopen('./data_sharing/initial_state.json','w'); fwrite(fid,jsonencode(input),'char'); fclose(fid);
!wsl ./build/bin/propAidaDebris
output=jsondecode(fileread('./data_sharing/out_prop_sec.json'));

x = output.constPart';
stm = permute(output.stm,[2,3,1]);

% propagate covariance
C=nan(6,6,N);
N_back=pp.N_back+Ntca-1;
C(:,:,N_back+1)=C0;
for i=N_back+1:-1:2
    invSTM=inv(stm(:,:,i)); C(:,:,i-1)=invSTM*C(:,:,i)*invSTM';
end
for i=N_back+2:N
    C(:,:,i)=stm(:,:,i)*C(:,:,i-1)*stm(:,:,i)';
end
pp.timeSubtr=pp.timeSubtr+output.timeMs/1000;
end
