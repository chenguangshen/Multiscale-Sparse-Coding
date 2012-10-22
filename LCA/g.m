% g.m - thresholding non-linearity
%
% function a = g(u,theta,type)

function a = g(u,theta,thresh_type)

if ~exist('thresh_type','var')
    thresh_type='soft';
end
    
M=size(u,1);

switch thresh_type
    
    case 'soft'
        a=abs(u)-repmat(theta,M,1);
        a(find(a<0))=0;
        a=sign(u).*a;
    case 'hard'
        a=u;
        a(find(abs(a)<repmat(theta,M,1)))=0;
    case 'hard+'
        a=u;
        a(find(a<repmat(theta,M,1)))=0;
        
end
