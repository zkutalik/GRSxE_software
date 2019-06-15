function [x,se,Cov] = regress_NOp(y,X,df)
% Model : y = X*beta
%%
nonsel = find(sum(isnan([y X]),2)==0);
X = X(nonsel,:);
y = y(nonsel,:);

NaN_field = find(max(X,[],1)==min(X,[],1));
noNaN_field = find(max(X,[],1)~=min(X,[],1));

if ~isempty(NaN_field)
    noNaN_field = [NaN_field(1) noNaN_field];
else
    y = y-mean(y);
end

x = NaN(size(X,2),1);
se = x;
X = X(:,noNaN_field);

if min(y)~=max(y)
    Z = inv(X'*X);
    tmp_x = Z*(X'*y);
    x(noNaN_field) = tmp_x;
else
    return
end

%%
if nargin<3
    df = length(tmp_x);
end

if nargout>1
    err = y - X*tmp_x;
    sig = sqrt(sum(err.^2,1)/(length(y)-df));
    tmp_se = sig.*sqrt(diag(Z));
    se(noNaN_field) = tmp_se;
    Cov = (sig.^2)*Z;
%     t = x./se;
%     p = 2*tcdfz(-abs(t),length(y)-df);
end

sel = find(imag(se)~=0);
x(sel) = NaN;
se(sel) = NaN;

