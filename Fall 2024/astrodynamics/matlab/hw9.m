%% Initialize
clear; clc;

%% Problem 5.8
%% Input times
years = [1914 1946 2010 2007 2024];
months = [8 4 9 10 11];
days = [14 18 1 16 12];
hours = [5 14 0 12 17];
mins = [30 0 0 0 0];
secs = [0 0 0 0 0];

%% calculate julian day w/ function from book
% convert to combined time
uts = hours + mins/60 + secs/3600;

% calc j0s
j0s = J0(years, months, days);

% then determine julian day
jds = j0s + uts/24;

%% output results
fprintf("Problem 5.8");
for i=1:length(years)
    fprintf("\nInput:");
    fprintf("\t%g:%g:%g UT on %g/%g/%g", hours(i), mins(i), secs(i), months(i), days(i), years(i));
    fprintf("\nJulian day:");
    fprintf("\t%11.3f\n", jds(i));
end

%% Problem 5.9
%% Input times (birthday first, then today)
years = [2002 2024];
months = [3 11];
days = [4 12];
hours = [12 12];
mins = [0 0];
secs = [0 0];

%% calculate julian day w/ function from book
% convert to combined time
uts = hours + mins/60 + secs/3600;

% calc j0s
j0s = J0(years, months, days);

% then determine julian day
jds = j0s + uts/24;

%% determine time between
jbday = jds(1);
jtday = jds(2);
days = jtday - jbday;
fprintf("\nProblem 5.9");
fprintf("\n%g days have passed between 12 UTC on my birthday and today\n", days);

%% Problem 5.10
%% Input locations and times
years = [2008 2007 2005 2006 2006 2024];
months = [1 12 7 2 3 11];
days = [1 21 4 15 21 12];
hours = [12 10 20 3 8 17];
mins = [0 0 0 0 0 0];
secs = [0 0 0 0 0 0];
elds = [18 144 -118+360 -43+360 131 -84+360]; % east long degs (convert west ones to east)
elms = [3 58 15 6 56 30];                     % east long mins
elss = [0 0 0 0 0 31.075];                    % east long mins

%% run local sidereal time function from book
% convert to decimals
ELs = elds + elms/60 + elss/3600;
uts = hours + mins/60 + secs/3600;
lsts = [0 0 0 0 0];
for i=1:length(years)
    lsts(i) = LST(years(i), months(i), days(i), uts(i), ELs(i));
end

%% output results
fprintf("\nProblem 5.10");
for i=1:length(years)
    fprintf("\nInput:");
    fprintf("\t%g:%g:%g UT on %g/%g/%g @ east long=%gdeg %g'", hours(i), mins(i), secs(i), months(i), days(i), years(i), elds(i), elms(i));
    fprintf("\nLocal sidereal time:");
    fprintf("\t%g deg\n", lsts(i));
end