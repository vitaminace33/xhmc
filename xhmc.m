function [q, accepted, N] = xHMC(fun, q0, options)

% default options
%      N = 1000; % Markov chain length
%   burn = 100;  % MCMC burn-in length
%   beta = 1;    % Boltzmann inverse temperature
%    psi = pi/3; % Horowitz's angle, psi = pi/2 gives HMC
%  extra = 3;    % Extra chances, extra = 0 gives (G)HMC
% integrator = @verlet;
%      h = .01;  % time step integration
%  steps = 10;   % integration steps
%  shift = .1;   % ±shift% for t.s.i. above (see `h' below)
% MaxInt = inf;  % Maximum number of integrations

eo = isempty(options);
if  eo || ~isfield(options,'N') || isempty(options.N)
    options.N = 1000;
end
if  eo || ~isfield(options,'burn') || isempty(options.burn)
    options.burn = 100;
end
if  eo || ~isfield(options,'beta') || isempty(options.beta)
    options.beta = 1;
end
if  eo || ~isfield(options,'psi') || isempty(options.psi)
    options.psi = pi/4;
end
if  eo || ~isfield(options,'extra') || isempty(options.extra)
    options.extra = 3;
end
if  eo || ~isfield(options,'integrator') || isempty(options.integrator)
    options.integrator = @verlet;
end
if  eo || ~isfield(options,'h') || isempty(options.h)
    options.h = .01;
end
if  eo || ~isfield(options,'steps') || isempty(options.steps)
    options.steps = 10;
end
if  eo || ~isfield(options,'shift') || isempty(options.shift)
    options.shift = .1;
end
if  eo || ~isfield(options,'MaxInt') || isempty(options.MaxInt)
    options.MaxInt = inf;
end

integrator = options.integrator;
 MaxInt = options.MaxInt;
   beta = options.beta;
    psi = options.psi;
      N = floor(abs(options.N));
   burn = floor(abs(options.burn));
      h = options.h;
  steps = options.steps;
  shift = options.shift;
  tries = options.extra+1;

cospsi = cos(psi);
sinpsi = sin(psi);

% dimension pre-processing
 dim = size(q0);
 Dim = prod(dim);
dims = length(dim);
isrow = false;
iscol = false;
if dims == 2
    if dim(1) == 1
        isrow = true;
    end
    if dim(2) == 1;
        iscol = true;
    end
end

% potential energy and forces definition (according to system dimensions)
% NOTE: f includes a leading `-', i.e. f=-DV
if iscol % || isscalar
    redim = dim;
    V = @(q)funV(fun,q);
    f = @(q)-funDV(fun,q);
    reshaped = false;
else
    redim = Dim;
    if isrow
        V = @(q)funV(fun,q');
        f = @(q)-funDV(fun,q')';
    else
        V = @(q)funV(fun,reshape(q,dim));
        f = @(q)reshape(-funDV(fun,reshape(q,dim)),[redim 1]);
    end
    reshaped = true;
end

% Silly code in case one wants to use a non-trivial mass matrix that, if it
% happens to be scalar, will simply use it as a scalar.
iM = eye(redim); % inverse mass matrix
lambda = trace(iM)/redim;
isscalar = false; %#ok<NASGU>
if all(all( iM == lambda*eye(redim) ))
        iM = lambda; % SPEEDS UP THINGS, CAREFUL!!!
        isscalar = true; %#ok<NASGU>
end

% initial burn-in phase
if burn > 0
    options.burn = 0;
    options.N = burn;
    options.MaxInt = Inf;
    q = xHMC(fun,q0,options); %#ok<NASGU>
    colons = ':';
    for i = 2:dims
        colons = [colons,',:']; %#ok<AGROW>[Dim 1]
    end
    eval(['q0 = q(',colons,',burn);'])
    clear q
end

irbeta = 1/sqrt(beta);

q0 = reshape(q0,[redim 1]);
p0 = zeros([redim 1]);
V0 = V(q0);

q = zeros([redim N]);

q1 = q0;
p1 = p0;
V1 = V0;
accepted = zeros(1,tries);
i = 0;
NumInt = 0;
while i < N && NumInt < MaxInt
    i = i+1;
    q0 = q1;
    p0 = p1;
    V0 = V1;
    p0 = cospsi*p0 + sinpsi*irbeta*randn([redim 1]);
    E0 = 1/2*(p0'*iM*p0) + V0;
    u = rand;
    trynum = 0;
    acceptance = 0;
    hplus = h*(1-shift+2*shift*rand);
    q1 = q0;
    p1 = p0;
    while acceptance < u && trynum < tries
        trynum = trynum + 1;
        [q1, p1] = integrator(iM,f,q1,p1,hplus,steps);
        V1 = V(q1);
        E1 = 1/2*(p1'*iM*p1) + V1;
        acceptance = min(1,max(acceptance,exp(-beta*(E1-E0))));
        E0 = E1;
    end
    if acceptance >= u
        accepted(trynum) = accepted(trynum) +1;
    else
        q1 =  q0;
        p1 = -p0;
        V1 =  V0;
    end
    q(:,i) = q1;
    NumInt = NumInt + trynum;
end
N = i;

if reshaped
    if isrow
        q = q';
    else
        q = reshape(q,[dim N]);
    end
end
end

function f = funV(fun,q)
[f,~] = fun(q);
end

function f = funDV(fun,q)
[~,f] = fun(q);
end

function [q1,p1] = verlet(iM,fun,q0,p0,h,steps)
q1 = q0;
pm = p0 + h/2*fun(q0); % fun(q0) could reuse fun(q1) from past integration
for i = 1:steps-1
    q1 = q1 + h*iM*pm;
    pm = pm + h*fun(q1);
end
q1 = q1 + h*iM*pm;
p1 = pm + h/2*fun(q1); % fun(q1) could be salvaged at next integration
end
