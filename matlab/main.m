clear all

% needs functions normalize(), regress_NOp(), IA_fit(), IA_fit_G2()


m = 100; % number of genetic markers
n = 1e4; % sample size
a0 = sqrt(.05); % linear effect of GRS on y
b1 = sqrt(.3); % linear effect of E on y
c1 = sqrt(.05); % interaction effect
d1 = 0; % correlation between E and GRS
skwE = 0; % skewness of E
krtE = 3; % kurtosis of E
skwN = 0; % skewness of the noise
krtN = 3; % kurtosis of the noise
pow = 1; % transformation power


sim_num = 1e2; % number of bootstrap / fake GRS
%% simulate data here - y (trait), GRS, E is masked

maf = rand(1,m);
G = binornd(2,ones(n,1)*maf,n,m);
eff = randn(m,1);
eff = eff/sqrt(sum(eff.^2));
GRS = normalize(G*eff,1);

E = d1*GRS+sqrt(1-d1^2)*pearsrnd(0,1,skwE,krtE,n,1);
noi = pearsrnd(0,1,skwN,krtN,n,1);
sig = sqrt(1-a0^2-b1^2-c1^2*(1+d1^2)-2*a0*b1*d1);

z = normalize(a0*GRS+b1*E+c1*(GRS.*E)+sig*noi);

if isnan(noi(1)) || isnan(E(1))
    error('Skew and kurtosis values not compatible (kurtosis > skew^2+1 not satisfied');
end
%%%%% transform trait
F = @(s,p1,p2) ((s-p1).^p2-1)/p2;

if pow~=0
    y = normalize(F(z,min(z)-1e-5,pow));
else
    y = normalize(log(z-min(z)+1e-5));
end

sel = find(abs(y)>10);
while ~isempty(sel)
    y(sel) = NaN;
    y = normalize(y);
    sel = find(abs(y)>10);
end


%% estimate interaction effect for GRS

[ao_het] = regress_NOp(y,[ones(length(y),1),GRS,GRS.^2]);

Xopt = zeros(4,sim_num);
% h = waitbar(0);
for j=1:sim_num
    j
%     waitbar(j/sim_num,h);

    ix = ceil(length(y)*rand(length(y),1)); % bootstrap sample index
    tmp = fminsearch(@(abc) IA_fit(abc,y(ix),GRS(ix)),[ao_het(2);0.1;0]); % first fit without quadratic GRS term
    Xopt(:,j) = fminsearch(@(abc) IA_fit_G2(abc,y(ix),GRS(ix)),[tmp(1);0;tmp(2:3)]);
end
% close(h);

%% fake GRS

fGRS = simulate_fGRS(y,GRS,sim_num);

Xopt0 = zeros(4,sim_num);
for j=1:sim_num
    j
    tmp = fminsearch(@(abc) IA_fit(abc,y,fGRS(:,j)),[ao_het(2);0.1;0]);
    tmp = fminsearch(@(abcd) IA_fit_G2(abcd,y,fGRS(:,j)),[tmp(1);0;tmp(2:3)]);
    Xopt0(:,j) = tmp;
end

%%
xopt0 = mean(Xopt0,2);
xopt = mean(Xopt,2);
SExopt0 = std(Xopt0,[],2);
SExopt = std(Xopt,[],2);
Pxopt0 = 2*normcdf(-abs(xopt0./SExopt0));
Pxopt = 2*normcdf(-abs(xopt./SExopt));

tdiff = (xopt0-xopt)./sqrt(SExopt.^2+SExopt0.^2);

[xopt xopt0]
[SExopt SExopt0]
[Pxopt Pxopt0]

