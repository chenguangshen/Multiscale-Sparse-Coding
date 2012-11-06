% sparsenet.m - Olshausen & Field sparse coding algorithm
% 
% Before running you must first load the training data array IMAGES

num_trials=10000;
batch_size=100;

[imsize imsize num_images]=size(IMAGES);
BUFF=4;

% number of outputs
M=144;

% number of inputs
N=144;
sz=sqrt(N);

% initialize basis functions (comment out these lines if you wish to 
% pick up where you left off)
Phi=randn(N,M);
Phi=Phi*diag(1./sqrt(sum(Phi.*Phi)));

% learning rate (start out large, then lower as solution converges)
eta = 3.0*(1/batch_size);

% lambda
lambda = 0.1;

a_var=ones(M,1);
var_eta=.1;

% Increase size to have four times the non-padded basis function. TO DO: Generalize to full image size
I=zeros(1*N,batch_size);

display_every=100;
display_network(Phi,a_var);

for t=1:num_trials
    
    % choose an image for this batch

    imi=ceil(num_images*rand);

    % extract subimages at random from this image to make data array I

    for i=1:batch_size
        % learns image patches that are twice the size
        r=BUFF+ceil((imsize-(2*sz)-2*BUFF)*rand);
        c=BUFF+ceil((imsize-(2*sz)-2*BUFF)*rand);
        I(:,i)=reshape(IMAGES(r:r+(1*sz)-1,c:c+(1*sz)-1,imi),4*N,1);
    end

    % calculate coefficients for these data via LCA

    % need to zero-pad Phi for the rest of the size

    %for i=1:N
    %    Phi(:,i) = reshape(padarray(reshape(Phi(:,i),sz,sz),[sz/2, sz/2]), 4*N, 1);
    %end

    ahat = sparsify(I,Phi,lambda);

    % calculate residual error

    R=I-Phi*ahat;

    % update bases
    dPhi = 0;
    
    for b=1:batch_size
        dPhi = dPhi + R(:,b) * ahat(:,b)'; %  learning rule here
    end

    %dPhi = dPhi/batch_size;

    Phi = Phi + eta * dPhi;

    Phi=Phi*diag(1./sqrt(sum(Phi.*Phi))); % normalize bases

    % accumulate activity statistics
    
    a_var=(1-var_eta)*a_var + var_eta*mean(ahat.^2,2);

    % display

    if (mod(t,display_every)==0)
        display_network(Phi,a_var);
    end

    % To-Do: Use simulated annealing
    eta = .9999 * eta;

end
