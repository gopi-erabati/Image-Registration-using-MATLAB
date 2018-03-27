function [e]=affine_registration_function(par,scale,Imoving,Ifixed,mtype,ttype)
% This function affine_registration_image, uses affine transfomation of the
% 3D input volume and calculates the registration error after transformation.
%
% I=affine_registration_image(parameters,scale,I1,I2,type);
%
% input,
%   parameters (in 2D) : Rigid vector of length 3 -> [translateX translateY rotate]
%                        or Affine vector of length 7 -> [translateX translateY
%                                           rotate resizeX resizeY shearXY shearYX]
%
%   parameters (in 3D) : Rigid vector of length 6 : [translateX translateY translateZ
%                                           rotateX rotateY rotateZ]
%                       or Affine vector of length 15 : [translateX translateY translateZ,
%                             rotateX rotateY rotateZ resizeX resizeY resizeZ,
%                             shearXY, shearXZ, shearYX, shearYZ, shearZX, shearZY]
%
%   scale: Vector with Scaling of the input parameters with the same lenght
%               as the parameter vector.
%   I1: The 2D/3D image which is affine transformed
%   I2: The second 2D/3D image which is used to calculate the
%       registration error
%   mtype: Metric type: sd: sum of squared differences m: mutual information.
%
% outputs,
%   I: An volume image with the registration error between I1 and I2
%
% example,
%
% Function is written by D.Kroon University of Twente (July 2008)
x=par.*scale;

switch ttype
    case 'r'
        M=[ cos(x(3)) sin(x(3)) x(1);
            -sin(x(3)) cos(x(3)) x(2);
            0 0 1];
    case 'a'
        M=[ x(3) x(4) x(1);
            x(5) x(6) x(2);
            0 0 1];
end

I3=affine_transform_2d_double(double(Imoving),double(M),3); % 3 stands for cubic interpolation

% metric computation
switch mtype
    case 'sd' %squared differences
        e=sum((I3(:)-Ifixed(:)).^2)/numel(I3);
    case 'm' %mutual information metric implementation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [COUNTS,X] = imhist(I3);
        prob_moving = COUNTS./numel(I3);
        %         prob_mov = [];
        %         for i=1:256
        %             if( prob_moving(i) == 0 )
        %                 %do nothing
        %             else
        %                 prob_mov = [prob_mov prob_moving(i)] ;
        %             end
        %         end
        
        [COUNTS,X] = imhist(Ifixed);
        prob_fixed = COUNTS./numel(Ifixed);
        %         prob_fix = [];
        %         for i=1:256
        %             if( prob_fixed(i)==0 )
        %                 %do nothing
        %             else
        %                 prob_fix = [prob_fix prob_fixed(i)] ;
        %             end
        %         end
        
        prob_mov_fix = zeros(255);
        I3 = im2uint8(I3);
        Ifixed = im2uint8(Ifixed);
        for i = 1:size(I3,1)
            for j = 1:size(I3,2)
                m = I3(i,j) + 1;
                n = Ifixed(i,j) + 1;
                prob_mov_fix(m,n) = prob_mov_fix(m,n) + 1;
            end
        end
        
        [r,c] = size(prob_mov_fix);
        prob_mov_fix = prob_mov_fix./(r*c);
        
        e = [];
        for i = 1:255
            for j = 1:255
                temp = prob_mov_fix(i,j)/(prob_moving(i)*prob_fixed(j));
                if (temp == 0 || temp == Inf || isnan(temp))
                    %do nothing
                else
                    
                    e = [e prob_mov_fix(i,j)*log10(temp)];
                end
            end
        end
        e = sum(e);
        
        % % calculating joint histogram for two images
        %         prob_mov_fix = zeros(256);
        %         I3 = im2uint8(I3);
        %         Ifixed = im2uint8(Ifixed);
        %         for i = 1:size(I3,1)
        %             for j = 1:size(I3,2)
        %                 m = I3(i,j) + 1;
        %                 n = Ifixed(i,j) + 1;
        %                 prob_mov_fix(m,n) = prob_mov_fix(m,n) + 1;
        %             end
        %         end
        % % prob_mov_fix
        %
        % [r,c] = size(prob_mov_fix);
        % jh_norm = prob_mov_fix./(r*c); % normalized joint histogram
        % y_marg = sum(jh_norm); %sum of the rows of normalized joint histogram
        % x_marg = sum(jh_norm,2);%sum of columns of normalized joint histogran
        %
        % Hy=0;
        % for i=1:c;    %  col
        %       if( y_marg(i)==0 )
        %          %do nothing
        %       else
        %          Hy = Hy + -(y_marg(i)*(log2(y_marg(i)))); %marginal entropy for image 1
        %       end
        % end
        % Hx=0;
        % for i=1:r;    %rows
        %    if( x_marg(i)==0 )
        %          %do nothing
        %       else
        %          Hx = Hx + -(x_marg(i)*(log2(x_marg(i)))); %marginal entropy for image 2
        %    end
        % end
        % h_xy = -sum(sum(jh_norm.*(log2(jh_norm+(jh_norm==0))))); % joint entropy
        % e = Hx + Hy - h_xy;% Mutual information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        error('Unknown metric type');
end;


