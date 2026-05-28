function gradientSym = gradient_sym(V,X,coordinate_system)
%% gradient_sym
% Copyright (c) 2018, mohamed fekry
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are 
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%POSSIBILITY OF SUCH DAMAGE.
%
%This Function calculates the gradient of 3D scalar function in Cartesian, Cylindrical, and Spherical coordinate system.
%function gradientSym = gradient_sym(V,X,coordinate_system)
%V is the 3D scalar function
%X is the parameter which the gradient will calculate with respect to.
%coordinate_system is the kind of coordinate system at which the vector field is specified.
%the gradient is calculated according to (Engineering Electromagnetics Sixth Edition William H. Hayt, Jr. . John A. Buck)
%Example (1):
%V=24*cos(pi*y/3)*sin(2*pi*z/3)
%gradient_sym(V,[x,y,z],'Cartesian')
%Example (2):
%V=Vo*exp(-2*r)*sin(3*phi)
%gradient_sym(V,[r,phi,z],'Cylindrical')
%Example (3):
%V=Vo*a/R*cos(2*th)
%gradient_sym(V,[R,th,phi],'Spherical')

%%
switch coordinate_system
    case {'cartesian','Cartesian'}
        gradientSym=gradient(V,X);
    case {'cylindrical','Cylindrical'}
        gradientSym=[diff(V,X(1)),1/X(1)*diff(V,X(2)),diff(V,X(3))];
    case {'spherical','Spherical'}
        gradientSym=[diff(V,X(1)),1/X(1)*diff(V,X(2)),1/(X(1)*sin(X(2)))*diff(V,X(3))];
end
end
