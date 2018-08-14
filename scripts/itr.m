function [ itr ] = itr(n, p, t)
% Calculate information transfer rate (ITR) for brain-computer interface 
% (BCI) [2]
% function [ itr ] = itr(n, p, t)
% 
% Input:
%   n   : # of targets
%   p   : Target identification accuracy (0 <= p <= 1) 
%   t   : Averaged time for a selection [s]
%
% Output:
%   itr : Information transfer rate [bits/min] 
%
% Reference:
%   [1] J. R. Wolpaw, N. Birbaumer, D. J. McFarland, G. Pfurtscheller, and
%       T. M. Vaughan,
%       "Brain¨Ccomputer interfaces for communication and control",
%       Clin. Neurophysiol. 113(6), 767¨C791, 2002.
% 
% Masaki Nakanishi & Yijun Wang, 19-June-2018
% E-mail: masaki@sccn.ucsd.edu, wangyj@semi.ac.cn

if nargin < 3
    error('stats:itr:LackOfInput', 'Not enough input arguments.'); end

if p < 0 || 1 < p
    error('stats:itr:BadInputValue',...
        'Accuracy need to be between 0 and 1.');
elseif p < 1/n
    warning('stats:itr:BadInputValue',...
        'The ITR might be incorrect because the accuracy < chance level.');
    itr = 0;
elseif p == 1
    itr = log2(n)*60/t;
else
    itr = (log2(n) + p*log2(p) + (1-p)*log2((1-p)/(n-1)))*60/t;
end