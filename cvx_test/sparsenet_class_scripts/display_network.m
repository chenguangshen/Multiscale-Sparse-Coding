function display_network(A,S_var)
%
%  display_network -- displays the state of the network (weights and 
%                     output variances)
%
%  Usage:
%
%    display_network(A,S_var);
%
%    A = basis function matrix
%    S_var = vector of coefficient variances
%

figure(1)

[L M]=size(A);

sz=sqrt(L);

buf=1;

if floor(sqrt(M))^2 ~= M
  n=sqrt(M/2);
  m=M/n;
else
  m=sqrt(M);
  n=m;
end

array=-ones(buf+n*(sz+buf),buf+m*(sz+buf));

k=1;

for j=1:m
  for i=1:n
    clim=max(abs(A(:,k)));
    array(buf+(i-1)*(sz+buf)+[1:sz],buf+(j-1)*(sz+buf)+[1:sz])=...
	reshape(A(:,k),sz,sz)/clim;
    k=k+1;
  end
end

subplot(212)
imagesc(array,[-1 1]), axis image off
title('basis functions')

if exist('S_var','var')
  subplot(211)
  bar(S_var), axis([0 M+1 0 max(S_var)])
  title('coeff variance')
end

drawnow
