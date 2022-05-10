% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for finding best number of Gaussian node for quadrature integral 
%%for given x grid, y grid, gaussian quadrature weight and density
%%distributions

function node=find_nodes(XX1,YY1,coef,rho)
    %input:
    %   XX1=gridded x data
    %   YY1=gridded y data
    %   coef=topographic surface data into Gaussian quadrature transformed
    %        grid
    
    %output
    %   node=required number of Gaussian quadrature node given density and
    %        topographic surface
    
            vv(:,:,1)=rand(size(XX1));
            for n=1:100
                ss1=0;
                [y,w]=lgwt(n,-1,1);
                for jn=1:n
                  ss1=ss1+w(jn).*rho(XX1,YY1,coef.*y(jn)+coef);
                end
                ss1=ss1.*coef;
                vv(:,:,n+1)=ss1;
                err=norm(vv(:,:,n+1)-vv(:,:,n));
                if err<=10^-5
                    break
                end
            end
            node=n;
            if node<4
                node=4;
            end
            %node=floor(0.5*node);
end
    
