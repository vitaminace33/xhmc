function [V, DV] = linearAlkanePE(q,params,varargin)
%LINEARALKANEPE Potential energy of a linear alkane molecule.
%   Consider a linear alkane molecule with n=m+2 carbon atoms.
%
%       CH_3---(CH_2)_m---CH_3
%
%   For a given set of system parameters, the potential is the sum of
%   several multi-body potentials: namely, the 2-body potential, (if m>0)
%   the 3-bodypotential, (if m>1) the 4-body potential and (if m>2) the
%   Lennard-Jones potential. Refer to [1] for specific details.
%
%   LINEARALKANEPE computes the value and the gradient of the potential
%   energy of a linear alkanel molecule (of arbitrary length n) given its
%   spatial configuration. The gradient is computed partially by regressive
%   accumulation.
%
%   Inputs:
%
%           q   Spatial configuration of each carbon atom. Array of size
%               [3,n], for wich q(:,i) is the space configuration of the
%               i-th carbon atom.
%      params   System parameters. Structure that, if not empty, must
%               include the following scalar fields: k0, d0, kth, th0, c1,
%               c2, c3, sig33, sig32, sig22, eps33, eps32, eps22. Per
%               default, if empty, it is set to the values given in [1].
%          lj   LJ potential switch (optional). Boolean that toggles the
%               activation of the LJ potential. Per default, it is set to
%               true (1).
%
%   Outputs:
%
%       V   Potential energy (scalar).
%      DV   Gradient of the potential energy. Array of size [3,n], n=m+2,
%           for which DV(j,i) is the partial derivative of V with respect
%           to q(j,i) evaluated at q.
%
%   Notes:
%
%       * The 2D input/ouput arrays can be easily transformed from/to 1D
%         arrays by the commands:
%             q2D = reshape(q1D,[3 n]);
%            DV1D = reshape(VD2D,[1 3*n]);
%
%   [1] E. Cances, F. Legoll, G. Stoltz: "Theoretical and Numerical
%   Comparison of some Sampling Methods for Molecular Dynamics", ESAIM:
%   M2AM (2007)
%
%   URL: http://google.com/+cedricmartinezcampos
%
%   Copyleft 2014 Cedric M. Campos, Dr., joint work with
%                 Mari Paz Calvo Cabrero, Prof. Dr.
%                 J. M. Sanz-Serna, Prof. Dr.
%
%   Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
%   http://creativecommons.org/licenses/by-nc-sa/4.0/
%
%   Supported by the European Union under the European Social Fund and the
%   Junta de Castilla y Leon (Spain).
%
%   $Revision: 0.2.0.00 $  $Date: 2014/07/25 10:45:34 $

narginchk(2, 3)
if nargin < 3
    lj = true;
else
    lj = varargin{1};
end

if ~isempty(params)
       k0 = params.k0;
       d0 = params.d0;
      kth = params.kth;
      th0 = params.th0;
       c1 = params.c1;
       c2 = params.c2;
       c3 = params.c3;
    eps33 = params.eps33;
    eps32 = params.eps32;
    eps22 = params.eps22;
    sig33 = params.sig33;
    sig32 = params.sig32;
    sig22 = params.sig22;
    if ~any([eps33 eps32 eps22])
        lj = false;
    end
else
       k0 = 1000;
       d0 = 1;
      kth = 208;
      th0 = 1.187;
       c1 = 1.18;
       c2 = -.23;
       c3 = 2.64;
    eps33 = .294;
    eps32 = .241;
    eps22 = .198;
    sig33 = 2.55;
    sig32 = 2.55;
    sig22 = 2.55;
end

dim = 3;
n = size(q,2);
N = n*(n-1)/2;

persistent index ij
if isempty(index) || size(index,2)~=n;
    index = zeros([n-1 n]); % hack: lowers computation, however
    ij = zeros([1,n-1]);    %       increases slightly memory load.
    for i = 1:n-1
        for j = i+1:n
            index(i,j) = index_(i,j,n);
        end
        ij(i) = index(i,i+1);
    end
end

I = eye(dim);
z = zeros([dim,1]);

 r = zeros([dim N]);
 d = zeros([1 N]);
Dd = zeros([dim N]);
if lj
    for i = 1:n-1
        for j = i+1:n
            k = index(i,j);
            r(:,k) = q(:,j)-q(:,i);
            d(k) = norm(r(:,k));
            Dd(:,k) = r(:,k)/d(k);
        end
    end
else
    for i = 1:n-1
        j = i+1;
        k = index(i,j);
        r(:,k) = q(:,j)-q(:,i);
        d(k) = norm(r(:,k));
        Dd(:,k) = r(:,k)/d(k);
    end
end
    

 un = zeros([dim n-1]);
Dun = zeros([dim dim n-1]);
for i = 1:n-1
    k = ij(i);
     un(:,i) = r(:,k)/d(k);
    Dun(:,:,i) = (I-un(:,i)*un(:,i)')/d(k);
end

 dotu = dot(un(:,1:n-2),un(:,2:n-1)); % bug: any(abs(dotu)>1) == 1
 dotu( dotu >  1) =  1; % bug workaround
 dotu( dotu < -1) = -1; % bug workaround

 theta = acos( dotu );
Dtheta = -1./sqrt(1-dotu.^2);
 
%  theta = zeros([1 n-2]);
% Dtheta = zeros([1 n-2]);
% for i = 1:n-2
%      theta(i) = acos( dotu(i) );
%     Dtheta(i) = -1/sqrt(1-dotu(i)^2);
% end

 V2 = 0;
DV2 = zeros([3 n-1]);
for i = 1:n-1
    [v,Dv] = vdos(d(ij(i)),d0,k0);
     V2 =  V2 +  v;
    DV2(:,i) = Dv*Dd(:,ij(i));
end

 V3 = 0;
DV3 = zeros([1 n-2]);
for i = 1:n-2
    [v,DV3(i)] = vtres(theta(i),th0,kth);
     V3 =  V3 +  v;
end
DthV3 = DV3.*Dtheta;
DV3 = zeros([dim n-1]);
    DV3(:,1)   =                        DthV3(1)*un(:,2);
for i = 2:n-2
    DV3(:,i)   = DthV3(i-1)*un(:,i-1) + DthV3(i)*un(:,i+1);
end
    DV3(:,n-1) = DthV3(n-2)*un(:,n-2);
for i = 1:n-1
    DV3(:,i) = Dun(:,:,i)*DV3(:,i);
end

crossr = zeros([dim n-2]);
hatr = zeros([dim dim n-1]);
hatr(:,:,1) = hat(r(:,1));
for i = 1:n-2
    hatr(:,:,i+1) = hat(r(:,ij(i+1)));
    crossr(:,i) = hatr(:,:,i)*r(:,ij(i+1));
end

%  un = zeros([dim n-2]);
% Dun = zeros([dim dim n-2]);
for i = 1:n-2
     ncrossri = norm(crossr(:,i));
     un(:,i) = crossr(:,i)/ncrossri;
    Dun(:,:,i) = (I-un(:,i)*un(:,i)')/ncrossri;
end

cosphi = -dot(un(:,1:n-3),un(:,2:n-2));

V4 = 0;
DV4 = zeros([1 n-3]);
for i = 1:n-3
    [v,DV4(i)] = vcuatro(cosphi(i),c1,c2,c3);
     V4 =  V4 +  v;
end
% DV4 = -DV4;              % minus sign is already taken into conideration
DuV4 = zeros([dim n-2]);   % in the block loop below
    DuV4(:,1)   =                      DV4(1)*un(:,2);
for i = 2:n-3
    DuV4(:,i)   = DV4(i-1)*un(:,i-1) + DV4(i)*un(:,i+1);
end
    DuV4(:,n-2) = DV4(n-3)*un(:,n-3);
for i = 1:n-2
    DuV4(:,i) = Dun(:,:,i)*DuV4(:,i);
end
DV4 = zeros([dim n-1]);
    DV4(:,1)   =                           - hatr(:,:,2)*DuV4(:,1);
for i = 2:n-2
    DV4(:,i)   = hatr(:,:,i-1)*DuV4(:,i-1) - hatr(:,:,i+1)*DuV4(:,i);
end
    DV4(:,n-1) = hatr(:,:,n-2)*DuV4(:,n-2);

 V =  V2 +  V3 +  V4;
DV = DV2 + DV3 + DV4;
DV = [z,DV] - [DV,z];

if lj
    VLJ = 0;
    DVLJ = zeros([dim n]);
    sigma = sig22;
    epsilon = eps22;
    for j = 4:n
        if j == n-1
            sigma = sig32;
            epsilon = eps32;
        elseif j == n
            sigma = sig33;
            epsilon = eps33;
        end
        for i = 1:n-(j-1)
            k = i+(j-1);
            ik = index(i,k);
            [v,Dv] = vlj(d(ik),sigma,epsilon);
            VLJ =  VLJ +  v;
            Dvik = Dv*Dd(:,ik);
            DVLJ(:,i) = DVLJ(:,i) - Dvik;
        DVLJ(:,k) = DVLJ(:,k) + Dvik;
        end
    end
     V =  V +  VLJ;
    DV = DV + DVLJ;
end
end

function [V, DV] = vdos(d, d0, k0)
 V = k0*(d-d0).^2/2;
DV = k0*(d-d0);
end

function [V, DV] = vtres(theta, th0, kth)
 V = kth*(theta-th0).^2/2;
DV = kth*(theta-th0);
end

function [V, DV] = vcuatro(x, c1, c2, c3)
 V = (1-x)*(c1+2*c2*(1+x)+c3*(1+4*x*(1+x)));
DV = -c1+3*c3-4*x*(c2+3*c3*x);
end

function [V, DV] = vlj(d,sigma,epsilon)
sd = (sigma/d).^6;
 V = 4*epsilon*sd.*(sd-1);
DV = 24*epsilon./d.*sd.*(1-2*sd);
end

function k = index_(i,j,n)
k = (i-1)*n-i*(i-1)/2+(j-i);
end

function a = hat(v)
a = [    0, -v(3),  v(2); ...
      v(3),     0, -v(1); ...
     -v(2),  v(1),     0];
end
