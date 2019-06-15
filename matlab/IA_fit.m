function v = IA_fit(abc,y,g)

a = abc(1,:);
b = abc(2,:);
c = abc(3,:);
s2 = (1-a.^2-b.^2-c.^2);
d1 = ones(length(y),1);
d2 = ones(1,length(a));

v0 = (1/2)*sum((y*d2-g*a).^2./((d1*b+g*c).^2 + d1*(s2))+log((d1*b+g*c).^2 + d1*(s2)),1); % minus log-likelihood
v = (1/2)*sum((y*d2).^2,1);
sel = find(s2>0);

% sel = find(min((d1*b+g*c).^2 + d1*(s2),[],1)>0);

v(sel) = v0(sel);

%%%%

