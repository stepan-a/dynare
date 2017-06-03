function runDynareReport(dc_a, dc_q, db_a, db_q)
% Copyright (C) 2013-2015 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

tic
larange= dates('2007a'):dates('2014a');
trange = dates('2012q2'):dates('2014q4');
prange = dates('2007q1'):dates('2013q4');
forecast_date = dates('2012q2');
srange = forecast_date:prange(end);

startpoint = strings(prange(1));
shaded = strings(srange(1));
endpoint = strings(prange(end));

shortNames = {'US', 'EU', 'JA', 'EA6', 'LA6', 'RC6'};
longNames  = {'Coca Cola', 'Kinder Bueno', 'Pizza', ...
              'Vegetarianism Is Good', 'OS X', 'Dothraki'};

%% Begin Report
rep = report();


%% Page 1: GDP
rep = rep.addPage('title', 'Jan1 vs Jan2', ...
                  'titleFormat', '\large\bfseries');
rep = rep.addSection();
rep = rep.addVspace();

% Table 1
rep = rep.addTable('title', {'Real GDP Growth','subtitle 1', 'subtitle 2'}, ...
                   'range', larange, ...
                   'vlineAfter', dates('2011y'), ...
                   'highlightRows', {'gray!25','white','green!22'});
rep = AnnualTable(rep, db_a, dc_a, 'PCH_GROWTH4_', larange);
rep = rep.addVspace('number', 2);

% Table 2
rep = rep.addTable('title', 'Potential GDP Growth', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
rep = AnnualTable(rep, db_a, dc_a, 'PCH_GROWTH4_BAR_', larange);


%% Page 2: Headline & Core Inflation
rep = rep.addPage('title', 'Jan1 vs Jan2', ...
                  'titleFormat', '\large\bfseries');
rep = rep.addSection();
rep = rep.addVspace();

% Table 1
rep = rep.addTable('title', 'Headline CPI Inflation', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
rep = AnnualTable(rep, db_a, dc_a, 'PCH_PIE4_', larange);
rep = rep.addVspace('number', 2);

% Table 2
rep = rep.addTable('title', 'Core CPI Inflation', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
rep = AnnualTable(rep, db_a, dc_a, 'PCH_PIEX4_', larange);


%% Page 3: Gas & Food Inflation
rep = rep.addPage('title', 'Jan1 vs Jan2', ...
                  'titleFormat', '\large\bfseries');
rep = rep.addSection();
rep = rep.addVspace();

% Table 1
rep = rep.addTable('title', 'Gas Inflation', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
rep = AnnualTable(rep, db_a, dc_a, 'PCH_PIE4_GAS_', larange);
rep = rep.addVspace('number', 2);

% Table 2
rep = rep.addTable('title', 'Food Inflation', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
rep = AnnualTable(rep, db_a, dc_a, 'PCH_PIE4_CONSFOOD_', larange);


%% Page 4: i & Output Gap
rep = rep.addPage('title', 'Jan1 vs Jan2', ...
                  'titleFormat', '\large\bfseries');
rep = rep.addSection();
rep = rep.addVspace();

% Table 1
rep = rep.addTable('title', 'Nominal Interest Rate', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
rep = AnnualTable(rep, db_a, dc_a, 'RS_', larange);
rep = rep.addVspace('number', 2);

% Table 2
rep = rep.addTable('title', 'Output Gap', 'range', larange, ...
                   'vlineAfter', dates('2011y'));
db_a = db_a.tex_rename('Y_WORLD', 'World');
rep = rep.addSeries('data', db_a{'Y_WORLD'});
delta = db_a{'Y_WORLD'}-dc_a{'Y_WORLD'};
delta = delta.tex_rename('$\Delta$');
rep = rep.addSeries('data', delta, ...
                    'tableShowMarkers', true, ...
                    'tableAlignRight', true);
rep = AnnualTable(rep, db_a, dc_a, 'Y_', larange);

%% Country Pages
for i=1:length(shortNames)
    rep = rep.addPage('title', {'Jan1 vs Jan2', longNames{i}}, ...
                      'titleFormat', {'\large\bfseries', '\large'});
    rep = rep.addSection('cols', 5);
    rep = CountryGraphPage(rep, shortNames{i}, db_q, dc_q, prange, srange);

    rep = rep.addPage('title', 'Jan1 vs Jan2', ...
                      'titleFormat', '\large\bfseries');
    rep = rep.addSection();
    rep = CountryTablePage(rep, shortNames{i}, longNames{i}, db_q, dc_q, ...
                           db_a, dc_a, trange, dates('2012q2'));
end

%% Residual Reports
% Countries
for i=1:length(shortNames)
    rep = rep.addPage('title', 'Residual Report Jan1 vs Jan2', ...
                      'titleFormat', '\large\bfseries');
    rep = rep.addSection();
    rep = ResidTablePage(rep, shortNames{i}, longNames{i}, db_q, dc_q, trange, dates('2012q2'));
end

% Commodities
rep = rep.addPage('title', 'Residual Report Jan1 vs Jan2', ...
                  'titleFormat', '\large\bfseries');
rep = rep.addSection();
rep = CommResidTablePage(rep, db_q, dc_q, trange, dates('2012q2'));

%% Commodities Graphs
%Page 24
rep = rep.addPage('title', 'Jan1 vs Jan2', ...
                  'titleFormat', '\large\bfseries');
rep = rep.addSection();

rep = rep.addGraph('title', {'World Real Oil Price Index','SUBTITLE'}, ...
                   'xrange', prange, ...
                   'shade', srange, ...
                   'xTicks', [1,5,10,15,find(srange(1)==prange),length(prange)], ...
                   'xTickLabels',{startpoint{:},'2008Q1','2009Q2','2010Q3',shaded{:}, endpoint{:}},...
                   'xTickLabelRotation', 0);
rep = rep.addSeries('data', db_q{'LRPFOOD_BAR_WORLD'}, ...
                    'graphBar', true, ...
                    'graphBarColor', 'red', ...
                    'graphBarFillColor', 'gray', ...
                    'graphBarWidth', 1);
db_q = db_q.tex_rename('LRPOIL_WORLD', 'Oil Price');
rep = rep.addSeries('data', db_q{'LRPOIL_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineWidth', 1.5, ...
                    'graphMarker', 'triangle*', ...
                    'graphMarkerEdgeColor','black', ...
                    'graphMarkerSize',4);
db_q = db_q.tex_rename('LRPOIL_BAR_WORLD', 'Equilibrium Oil Price');
rep = rep.addSeries('data', db_q{'LRPOIL_BAR_WORLD'}, ...
                    'graphLineColor', 'green', ...
                    'graphLineStyle', 'solid', ...
                    'graphLineWidth', 1.5);


rep = rep.addGraph('title', 'World Real Food Price Index', ...
                   'xrange', prange, ...
                   'shade', srange, ...
                   'xTicks', [1,5,10,15,find(srange(1)==prange),length(prange)], ...
                   'xTickLabels',{startpoint{:},'2008Q1','2009Q2','2010Q3',shaded{:}, endpoint{:}},...
                   'xTickLabelRotation', 0, ...
                   'showLegend', true, ...
                   'legendAt', [.5,.5]);
rep = rep.addSeries('data', db_q{'LRPFOOD_BAR_WORLD'}, ...
                    'graphBar', true, ...
                    'graphBarColor', 'green', ...
                    'graphBarFillColor', 'yellow', ...
                    'graphBarWidth', 1);
db_q = db_q.tex_rename('LRPFOOD_WORLD', 'Food Price');
rep = rep.addSeries('data', db_q{'LRPFOOD_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineWidth', 1.5);
db_q = db_q.tex_rename('LRPFOOD_BAR_WORLD', 'Equilibrium Food Price');
rep = rep.addSeries('graphVline', dates('2009q2'), ...
                    'graphLineColor', 'red', ...
                    'graphLineWidth', 1.5);

% Page 25
rep = rep.addPage('title', {'Jan1 vs Jan2', 'World Oil and Food Prices'}, ...
                  'titleFormat', {'\large\bfseries', '\large'});
rep = rep.addSection('cols', 1);
rep = rep.addParagraph('text', 'Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.', ...
                       'cols', 2, ...
                       'heading', '\textbf{My First Paragraph Has Two Columns}');

rep = rep.addSection('cols', 1);
rep = rep.addParagraph('text', 'Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\newline', ...
                       'heading', '\textbf{My Next Paragraphs Only Have One}', ...
                       'indent', false);
rep = rep.addParagraph('text', 'Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.\newline');

rep = rep.addSection('cols', 2);

rep = rep.addGraph('title', 'World Real Oil Price', ...
                   'xrange', prange, ...
                   'shade', srange, ...
                   'xTicks', [1,5,10,15,find(srange(1)==prange),length(prange)], ...
                   'xTickLabels',{startpoint{:},'2008Q1','2009Q2','2010Q3',shaded{:}, endpoint{:}},...
                   'xTickLabelRotation', 0);
rep = rep.addSeries('data', db_q{'LRPOIL_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineWidth', 1.5);
rep = rep.addSeries('data', dc_q{'LRPOIL_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineStyle', 'dashed', ...
                    'graphLineWidth', 1.5);

srange1 = prange(1):forecast_date;
rep = rep.addGraph('title', 'Equilibrium World Real Oil Price', ...
                   'xrange', prange, ...
                   'shade', srange1);
rep = rep.addSeries('data', db_q{'LRPOIL_BAR_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineWidth', 1.5);
rep = rep.addSeries('data', dc_q{'LRPOIL_BAR_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineStyle', 'dashed', ...
                    'graphLineWidth', 1.5);

rep = rep.addGraph('title', 'World Real Food Price', ...
                   'xrange', prange, ...
                   'shade', srange, ...
                   'xTickLabels','ALL',...
                   'xTickLabelRotation', 45,...
                   'xAxisTight',false,...
                   'yAxisTight',true);
rep = rep.addSeries('data', db_q{'LRPFOOD_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineWidth', 1.5);
rep = rep.addSeries('data', dc_q{'LRPFOOD_WORLD'}, ...
                    'graphLineColor', 'blue', ...
                    'graphLineStyle', 'dashed', ...
                    'graphLineWidth', 1.5);
rep = rep.addSeries('graphHline', 460, ...
                    'graphLineColor', 'red', ...
                    'graphLineWidth', 1.5);

a=dseries([1:200]', '1984q1');
b=a;
c=a;
d=a;
b(dates('2012q2'):dates('2015q2'))=b(dates('2012q2'):dates('2015q2'))+2;
c(dates('2012q2'):dates('2015q2'))=c(dates('2012q2'):dates('2015q2'))+4;
d(dates('2012q2'):dates('2015q2'))=d(dates('2012q2'):dates('2015q2'))+6;

rep = rep.addGraph('title', 'Equilibrium World Real Food Price', ...
                   'xrange', prange, ...
                   'shade', srange, ...
                   'showLegend', true, ...
                   'xTickLabelRotation', 0);
rep = rep.addSeries('data', a, ...
                    'graphLineColor', 'blue', ...
                    'graphLineWidth', 1.5, ...
                    'graphLegendName', 'baseline', ...
                    'graphMiscTikzAddPlotOptions', 'mark=halfcircle*,color=red');
rep = rep.addSeries('data', b, ...
                    'graphLineColor', 'blue', ...
                    'graphLineStyle', 'dashed', ...
                    'graphLineWidth', 1.5, ...
                    'graphLegendName', 'control', ...
                    'graphMiscTikzAddPlotOptions', 'mark=halfcircle*,mark options={rotate=90,scale=3}', ...
                    'graphFanShadeColor', 'red', 'graphFanShadeOpacity', 40);
rep = rep.addSeries('data', c, ...
                    'graphLineColor', 'blue', ...
                    'graphLineStyle', 'dashed', ...
                    'graphLineWidth', 1.5, ...
                    'graphLegendName', 'control', ...
                    'graphFanShadeColor', 'red', 'graphFanShadeOpacity', 30);
rep = rep.addSeries('data', d, ...
                    'graphLineColor', 'blue', ...
                    'graphLineStyle', 'dashed', ...
                    'graphLineWidth', 1.5, ...
                    'graphLegendName', 'control', ...
                    'graphFanShadeColor', 'red', 'graphFanShadeOpacity', 20);

%% Write & Compile Report
rep.write();
rep.compile();
toc
end
