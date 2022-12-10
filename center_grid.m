% ***************************************************************
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Solid Earth Research Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************

%%Matlab function for creating center grid of any data using mean of adjacent grids
function [x_g,y_g,data_g]=center_grid(x,y,data)
    %input:
    %   x=all x data values having dimension (m+1 X 1)
    %   y=all y data values having dimension (n+1 X 1)
    %   data= 2D gridded data having dimension (m+1 X n+1)
    
    %output
    %   x=all x data center grid values having dimension (m X 1)
    %   y=all y data center grid values having dimension (n X 1)
    %   data= 2D center gridded data having dimension (m X n)
    
    %grid spacing
    dx=x(2)-x(1);
    dy=y(2)-y(1);
    %length of x and y dimensions
    m=length(x); n=length(y);
    %size of the data
    sz=size(data);
    %checking for dimensions of data and grid points
    vv=(m-sz(2))+(n-sz(1));
    if vv==0
        %taking mean of adjacent grids
        A1=data(1:n-1,1:m-1);
        A2=data(2:n,1:m-1);
        A3=data(1:n-1,2:m);
        A4=data(2:n,2:m);
        
        data_g=(A1+A2+A3+A4)/4;
        %new center grids x and data
        %x_g=x(1)+(dx/2):dx:x(end)-(dx/2);
        %y_g=y(1)+(dy/2):dy:y(end)-(dx/2);
        x_g=linspace(x(1)+(dx/2),x(end)-(dx/2),m-1);
        y_g=linspace(y(1)+(dy/2),y(end)-(dy/2),n-1);
    else
        msg = 'dimension mismatch of x, y and corresponding data';
        error(msg)
    end
end
