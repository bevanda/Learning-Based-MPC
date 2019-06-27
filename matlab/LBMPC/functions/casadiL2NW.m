function g = casadiL2NW(x,u,data)
X=data(1:3,:); Y=data(4:7,:); v=data(8,:); ksi=[x(1:2);u];
if nargin <2
    X=zeros(3,1);
    Y=zeros(4,1);
end 

bandwidth = 0.5;  
lambda = 0.001; 
s = size(Y,1);
n = size(X,2);
kval = cell(1,n);

for i=1:n
    kval{i} = exp(-((norm(X(:,i) - ksi).^2)/bandwidth^2));
end

skval = 0;
for j=1:n
    skval = skval + kval{j}*v(j);
end

y=zeros(s,1);
for k=1:n
    out = Y(:,k)*(kval{k}/(lambda + skval));
    y = y + out;
end
g = y;
end

