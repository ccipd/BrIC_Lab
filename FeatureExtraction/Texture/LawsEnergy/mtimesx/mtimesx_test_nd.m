% Test routine for mtimesx, multi-dimensional speed and equality to MATLAB
%******************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    mtimesx_test_nd
%  Filename:    mtimesx_test_nd.m
%  Programmer:  James Tursa
%  Version:     1.11
%  Date:        January 6, 2010
%  Copyright:   (c) 2009,2010 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
%  Syntax:
%
%    mtimesx_test_nd    % default n=3 is used
%    mtimesx_test_nd(n)
%
%    where n = number of repetitions (should be 3 or greater)
%
%  Output:
%
%    Prints out speed and equality test results.
%
%--------------------------------------------------------------------------

function mtimesx_test_nd(n)
mtimesx; % load the mex routine into memory
if( nargin == 0 )
    n = 3;
else
    n = floor(n);
    if( ~(n >= 3 && n <= 100) )
        n = 3;
    end
end
cn = sprintf('%g',n);

disp(' ');
disp(['mtimesx multi-dimensional test routine using ' cn ' repetitions']);

disp(' ');
disp('3x7x3x4x3x2x5x8 example');
A = rand(3,5,1,4,3,2,1,8);
B = rand(5,7,3,1,3,2,5);
% mtimes
tm = zeros(1,n);
for k=1:n
clear Cm Cx
A(1) = 2*A(1);
B(1) = 2*B(1);
tic
Cm = zeros(3,7,3,4,3,2,5,8);
for k1=1:3
    for k2=1:4
        for k3=1:3
            for k4=1:2
                for k5=1:5
                    for k6=1:8
                        Cm(:,:,k1,k2,k3,k4,k5,k6) = A(:,:,1,k2,k3,k4,1,k6) * B(:,:,k1,1,k3,k4,k5);
                    end
                end
            end
        end
    end
end
tm(k) = toc;
end
% mtimesx
tx = zeros(1,n);
for k=1:n
tic
Cx = mtimesx(A,B);
tx(k) = toc;
end
% results
tm = median(tm);
tx = median(tx);
faster = sprintf('%5.1f',100*(tm-tx)/tx-100);
disp(' ');
disp(['mtimesx is ' faster '% faster than mtimes for 3x7x3x4x3x2x5x8 example'])
if( isequal(Cx,Cm) )
    disp('mtimesx result matches mtimes for 3x7x3x4x3x2x5x8 example')
else
    disp('mtimesx result does not match mtimes for 3x7x3x4x3x2x5x8 example')
end

disp(' ');
disp('3x3x1000000 example');
A = rand(3,3,1000000);
B = rand(3,3,1000000);
% mtimes
tm = zeros(1,n);
for k=1:n
clear Cm Cx
A(1) = 2*A(1);
B(1) = 2*B(1);
tic
Cm = zeros(3,3,1000000);
for k1=1:1000000
    Cm(:,:,k1) = A(:,:,k1) * B(:,:,k1);
end
tm(k) = toc;
end
% mtimesx
tx = zeros(1,n);
for k=1:n
tic
Cx = mtimesx(A,B);
tx(k) = toc;
end
% results
tm = median(tm);
tx = median(tx);
faster = sprintf('%5.1f',100*(tm-tx)/tx-100);
disp(' ');
disp(['mtimesx is ' faster '% faster than mtimes for 3x3x1000000 example'])
if( isequal(Cx,Cm) )
    disp('mtimesx result matches mtimes for 3x3x1000000 example')
else
    disp('mtimesx result does not match mtimes for 3x3x1000000 example')
end

disp(' ');
disp('Done');
disp(' ');

end
