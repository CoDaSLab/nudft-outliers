function [data, outliers] = nudft_filter(tbl,smplR,p_cutoff)
%% Identifies outliers in non-uniformly spaced, but periodic data. Input p value for statistical cutoff, p is the expected fraction of data to exclude.
% smplR is expected sample frequency.
% tbl is the data read from a .csv while where the first column corresponds to the time, and the second column the measurements

% coded by: Michael Sorochan Armstrong (mdarmstr@ugr.es), and Jose Camacho
% Paez (jcamachop@ugr.es)
% last modification: 22/Dec/2022
%
% Copyright (C) 2023  University of Granada, Granada
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

T = readtable(tbl,"NumHeaderLines",1);
X = table2array(T(:,2));

L = size(X,1);
t = (0:L-1)*smplR; %length of time vector

outliers = false(L,1);
logical_vec = X~=-9999; %Ignore missing values

chi2cutoff = chi2inv(1-p_cutoff,1);

mu = mean(X(logical_vec));

Y = nufft(X(logical_vec)-mu,find(logical_vec),t/max(t));

ftfit = zeros(size(Y));

[~,indices] = maxk(abs(Y),6);
ftfit(indices) = Y(indices);

y = ifft(ftfit,"symmetric") + mu;

Z = (X(logical_vec) - y(logical_vec)) ./ std(X(logical_vec) - y(logical_vec));

outliers(logical_vec) = Z.^2 < chi2cutoff;
data = X;
