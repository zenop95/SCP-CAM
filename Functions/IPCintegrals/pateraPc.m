function Pc = pateraPc(r,v,P,R)
% pateraPc computes the non linear collision probability as taking into
% account the finite dimension of the ellipsoid on the direction of the
% relative velocity, according to the method developed by Patera.

% INPUT: r  = [m]   3xN relative position array
%        v  = [m/s] 3xN relative velocity array
%        P  = [m^2] 3x3xN Covariance Matrix of the relative position
%        R  = [m]   Hard Body Radius (HBR) of the collision

% OUTPUT: Ipc = [-] Computed Instantaneous Collision Probability


% Bibliography: Patera, R. P. (2003). Satellite Collision Probability for 
% Nonlinear Relative Motion. Journal of Guidance, Control, and Dynamics, 
% 26(5), 728–733. https://doi.org/10.2514/2.5127

% Author: Zeno Pavanello, 2022
%--------------------------------------------------------------------------
N = size(r,2);
validateattributes(r,{'double'},{'2d','nrows',3,'ncols',N})
validateattributes(v,{'double'},{'2d','nrows',3,'ncols',N})
validateattributes(R,{'double'},{'scalar','positive'})
n   = 1000;
Pci = nan(N,1);
for i = 1:N
%     try chol(P(:,:,i));
        [V,D] = eig(P(:,:,i));
        if det(V) < 0; V(:,1) = -V(:,1); end
        [a,b] = sort(diag(D),'ascend');
        D = diag(a);
        V = V(:,b);
        a = sqrt(D(1,1)); b = sqrt(D(2,2)); c = sqrt(D(3,3));
        S = diag([1 a/b a/c]);
        A = S*V';
        rN = A*r(:,i);
        vN = A*v(:,i);
        % Encounter frame
        z = vN/norm(vN);
        x = cross(vN,rN); x = x/norm(x);
        y = cross(z,x); y = y/norm(y);
        ToEncPlane = [x'; y'; z'];
    
        % Rotate vectors and matrix        
        rEnc     = ToEncPlane*rN;
        vEnc     = ToEncPlane*vN; 
        r2d      = rEnc(1:2);
       
        %Define the vectors for the integration on the 2D plane
        sphere    = R*eye(3); %Hard Body sphere
        stretched = ToEncPlane*A*sphere*A'*ToEncPlane'; %stretch and rotate to bring in the same frame as covariance
        hbr2d     = stretched(1:2,1:2);         %Take the 2d part excluding the component in the velocity direction
        [F,G]     = eig(hbr2d);
        zeroIsIn  = isInBox(0,0,G(1,1),G(2,2),-r2d); % check if the ellipse contains the origin (before transformation because it is easier)
        E         = linspace(0,2*pi,n);      % Discretize the ellipise in its own reference frame
        disc      = [G(1,1)*cos(E); G(2,2)*sin(E)];
        for j = 1:n
            discRot(:,j) = F*disc(:,j); %Bring the discretized ellipse back to the covariance frame
        end
        X = discRot + r2d; %Add the position as an offset to the HBR sphere

        %Probability of collision linear and non linear
        int      = 0;
        for j = 1:n-1
            int = int + exp(-(norm(X(:,j))^2+norm(X(:,j+1))^2)/(4*a^2))*...
                        asin(norm(cross([X(:,j);0],[X(:,j+1);0])) ...
                        /norm(X(:,j)*norm(X(:,j+1))));
        end
        if zeroIsIn
            Pc2d   = 1 - int/(2*pi);
        else
            Pc2d   = -int/(2*pi);
        end
        Pci(i) = vEnc(3)/(a*sqrt(2*pi))*exp(-rEnc(3)^2/(2*a^2))*Pc2d;
%     catch
%     error('the covariance matrix is not positive semi-definite')
%     end
end
Pc = sum(Pci);
end
