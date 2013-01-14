% sparsify.m - sparsifies analog coefficients via LCA
function a = sparsify(I,Phi,lambda,thresh_type,display_p)

sz = 32;
szb = 16;
OC = 2;
M = OC*(sz-(szb-1))^2;

if ~exist('display_p','var')
    display_p=0;
end

if ~exist('thresh_type','var')
    thresh_type='soft';
end

if size(Phi',1) ~= size(I,1)
    Phihat = reshape(Phi,szb,szb,M);

    Phi = zeros(sz,sz,M);

    for i=1:M/OC
        %current row
        j = ceil(i/(sz-(szb-1)));

        %current column
        k = mod(i,sz-(szb-1));
        if k == 0
            k = sz-(szb-1);
        end 

        for z=1:OC
            Phi(j:j+szb-1,k:k+szb-1,i+(z-1)*M/OC) = Phihat(:,:,i+(z-1)*M/OC);
        end
    end

    Phi = reshape(Phi,sz*sz,M);
end

batch_size=size(I,2);

[N M]=size(Phi);
sz=sqrt(N);

b=Phi'*I;
G=Phi'*Phi-eye(M);

num_iterations=75;
eta=0.1;

u = zeros(M,batch_size);

l=0.5*max(abs(b));
a=g(u,l);

if display_p
    figure(2)
    subplot(411)
    ha=bar(a(:,1)); axis([0 M+1 -2 2])
    subplot(412)
    hu=bar(u(:,1)); axis([0 M+1 -2 2])
    hold on
    hlu=plot([0 M+1],[l l]);
    hll=plot([0 M+1],-[l l]);
    hold off
    subplot(413)
    hIh=imagesc(zeros(sz),'EraseMode','none',[-1 1]); axis image
    subplot(414)
    imagesc(reshape(I(:,1),sz,sz),[-1 1]), axis image
    drawnow
end

for t=1:num_iterations
    
    u = eta*(b-G*a) + (1-eta)*u;
    a = g(u,l,thresh_type);

    if display_p
        set(ha,'YData',a(:,1))
        set(hu,'YData',u)
        set(hlu,'YData',[l l])
        set(hll,'YData',-[l l])
        Ihat=Phi*a(:,1);
        set(hIh,'CData',reshape(Ihat,sz,sz))
        drawnow
    end
    
    l=0.95*l;
    l(find(l<lambda))=lambda;
    
end
