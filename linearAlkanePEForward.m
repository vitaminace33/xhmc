function [V, DV] = linearAlkanePEForward(q,params,varargin)
%LINEARALKANEPEFORWARD Potential energy of a linear alkane molecule.
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
%   spatial configuration. The gradient is computed completely by forward
%   accumulation. Slow but clear explicit implementation.
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
%   $Revision: 0.3.0.00 $  $Date: 2014/07/25 10:45:47 $

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

index = zeros([n-1 n]); % hack: lowers computation, however
for i = 1:n-1           %       increases slightly memory load.
    for j = i+1:n
        index(i,j) = index_(i,j,n);
    end
end

 r = zeros([dim N]);
Dr = zeros([N n]); % hack: should be [N dim dim n], but then the
for i = 1:n-1      %       +/-1 below should be +/-eye(dim). This
    for j = i+1:n  %       saves some memory and cpu load.
        ij = index(i,j);
         r(:,ij) = q(:,j)-q(:,i);
        Dr(ij,i) = -1;
        Dr(ij,j) =  1;
    end
end

 d = zeros([1 N]);
Dd = zeros([N dim n]);
for i = 1:n-1
    for j = i+1:n
        ij = index(i,j);
        d(ij) = norm(r(:,ij));
        Dd(ij,:,i) = r(:,ij)/d(ij)*Dr(ij,i);
        Dd(ij,:,j) = r(:,ij)/d(ij)*Dr(ij,j);
    end
end

 un = zeros([dim n-1]);
Dun = zeros([n-1 dim dim n]);
 Id = eye(dim);
for i = 1:n-1
    ij = index(i,i+1);
     un(:,i) = r(:,ij)/d(ij);
    Duni = (Id-un(:,i)*un(:,i)')/d(ij);
    for k = i:i+1
        Dun(i,:,:,k) = Duni*Dr(ij,k);
    end
end

 dotu = dot(un(:,1:n-2),un(:,2:n-1));
Ddotu = zeros([n-2 dim n]);
for i = 1:n-2
    for j = 1:dim
        for k = i:i+2
            Ddotu(i,j,k) = Dun(i  ,:,j,k)*un(:,i+1) + ...
                           Dun(i+1,:,j,k)*un(:,i  );
        end
    end
end
clear un Dun

 theta = zeros([1 n-2]);
Dtheta = zeros([n-2 dim n]);
for i = 1:n-2
     theta(i) = acos( dotu(i) );
    Dtheta(i,:,:) = -1/sqrt(1-dotu(i)^2)*Ddotu(i,:,:);
end
clear dotu Ddotu

 crossr = zeros([dim n-2]);
Dcrossr = zeros([n-2 dim dim n]);
% Id = eye(dim);
jk = 1;
crossrjk = hat(r(:,jk));
for i = 1:n-2
    ij = jk;
    jk = index(i+1,i+2);
%     crossr(:,i) = cross(r(:,ij),r(:,jk));
%     for j = 1:dim
%         for k = 1:n
%             Dcrossr(i,:,j,k) = Dr(ij,k)*cross(Id(:,j),r(:,jk)) ...
%                                + Dr(jk,k)*cross(r(:,ij),Id(:,j));
%         end
%     end
    crossrij = crossrjk;
    crossrjk = hat(r(:,jk));
    crossr(:,i) = crossrij*r(:,jk);
    for k = i:i+2
        Dcrossr(i,:,:,k) = Dr(jk,k)*crossrij - Dr(ij,k)*crossrjk;
    end
end
clear r Dr

 un = zeros([dim n-2]);
Dun = zeros([n-2 dim dim n]);
 Id = eye(dim);
for i = 1:n-2
     ncrossri = norm(crossr(:,i));
     un(:,i) = crossr(:,i)/ncrossri;
    Duni = (Id-un(:,i)*un(:,i)')/ncrossri;
    for k = i:i+2
      Dun(i,:,:,k) = Duni*reshape(Dcrossr(i,:,:,k),[dim dim]);
    end
end
clear crossr Dcrossr

 dotu = dot(un(:,1:n-3),un(:,2:n-2));
Ddotu = zeros([n-3 dim n]);
for i = 1:n-3
    for j = 1:dim
        for k = i:i+3
            Ddotu(i,j,k) = Dun(i  ,:,j,k)*un(:,i+1) + ...
                           Dun(i+1,:,j,k)*un(:,i  );
        end
    end
end
clear un Dun

 cosphi =  -dotu;
Dcosphi = -Ddotu;
clear dotu Ddotu

 V2 = 0;
DV2 = zeros([dim n]);
for i = 1:n-1
    ij = index(i,i+1);
    [v,Dv] = vdos(d(ij),d0,k0);
     V2 =  V2 +  v;
    DV2 = DV2 + Dv*reshape(Dd(ij,:,:),[dim n]);
end

 V3 = 0;
DV3 = zeros([dim n]);
for i = 1:n-2
    [v,Dv] = vtres(theta(i),th0,kth);
     V3 =  V3 +  v;
    DV3 = DV3 + Dv*reshape(Dtheta(i,:,:),[dim n]);
end

 V4 = 0;
DV4 = zeros([dim n]);
for i = 1:n-3
    [v,Dv] = vcuatro(cosphi(i),c1,c2,c3);
     V4 =  V4 +  v;
    DV4 = DV4 + Dv*reshape(Dcosphi(i,:,:),[dim n]);
end

 VLJ = 0;
DVLJ = zeros([dim n]);
if lj
%   for i = 1:n-4
%       for j = i+3:n
%           if i < n-1 && j < n-1
%               sigma = sig22;saves
%               epsilon = eps22;
%           elseif i == n && j == n
%               sigma = sig33;
%               epsilon = eps33;
%           else
%               sigma = sig32;
%               epsilon = eps32;
%           end
%           ij = index(i,j);
%           [v,Dv] = vlj(d(ij),sigma,epsilon);
%            VLJ =  VLJ +  v;
%           DVLJ = DVLJ + Dv*squeeze(Dd(ij,:,:));
%       end
%   end
%   % hack: reduces checkings from the previous loops
%   %       works along "diagonals" according to (i,j)
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
            ij = index(i,j+(i-1));
            [v,Dv] = vlj(d(ij),sigma,epsilon);
             VLJ =  VLJ +  v;
            DVLJ = DVLJ + Dv*reshape(Dd(ij,:,:),[dim n]);
        end
    end
end

 V =  V2 +  V3 +  V4 +  VLJ;
DV = DV2 + DV3 + DV4 + DVLJ;
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
