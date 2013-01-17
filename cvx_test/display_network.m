function display_network(A,A2,S_var)
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

[L2 M2]=size(A2);
sz2=sqrt(L2);

buf=1;

if floor(sqrt(M))^2 ~= M
  n=sqrt(M/2);
  m=M/n;
else
  m=sqrt(M);
  n=m;
end

array=-ones(buf+n*(sz+buf),buf+m*(sz+buf));
array2=-ones(buf+n*(sz2+buf),buf+m*(sz2+buf));

k=1;

for j=1:m
  for i=1:n
    clim=max(abs(A(:,k)));
    array(buf+(i-1)*(sz+buf)+[1:sz],buf+(j-1)*(sz+buf)+[1:sz])=...
	reshape(A(:,k),sz,sz)/clim;
    k=k+1;
  end
end

k=1;

for j=1:m
  for i=1:n
    clim=max(abs(A2(:,k)));
    array2(buf+(i-1)*(sz2+buf)+[1:sz2],buf+(j-1)*(sz2+buf)+[1:sz2])=...
  reshape(A2(:,k),sz2,sz2)/clim;
    k=k+1;
  end
end

if exist('S_var','var')
  subplot(211)
  bar(S_var), axis([0 M+1 0 max(S_var)])
  title('Coefficient Variance')
  subplot(212)
end

subplot(121)
imagesc(array,[-1 1]), axis image off

title('Basis Functions L1')

subplot(122)
imagesc(array2,[-1 1]), axis image off

title('Basis Functions L2')

drawnow
