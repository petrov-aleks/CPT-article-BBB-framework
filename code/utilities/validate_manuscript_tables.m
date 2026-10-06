function validate_manuscript_tables()
%VALIDATE_MANUSCRIPT_TABLES Compare displayed CSV values with manuscript TeX.
% This function is read-only: it never generates or modifies a TeX table.

codeDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(codeDir));
csvDir = fullfile(repoRoot, 'results', 'tables');
texDir = fullfile(csvDir, 'tex');

failures = strings(0,1);
failures = [failures; validateTable1(csvDir, texDir)];
failures = [failures; validateS1(csvDir, texDir)];
failures = [failures; validateS2(csvDir, texDir)];
failures = [failures; validateS3(csvDir, texDir)];
failures = [failures; validateS4(csvDir, texDir)];
failures = [failures; validateS5(csvDir, texDir)];

if ~isempty(failures)
    fprintf(2, '%s\n', join(failures, newline));
    error('Table validation failed with %d numerical mismatch(es).', numel(failures));
end
fprintf('PASS: all displayed numerical values in Tables 1 and S1-S5 match their CSV sources.\n');
end

function failures = validateTable1(csvDir, texDir)
t = readtable(fullfile(csvDir,'table_1.csv'),'TextType','string');
rows = dataRows(fullfile(texDir,'permeability_table.tex'), 6, 3, true);
expected = strings(height(t),3);
expected(:,1) = displayNumber(t.PS0_1e4_mL_per_g_per_s);
expected(:,2) = displayNumber(t.P0_1e6_cm_per_s);
expected(:,3) = displayNumber(t.PS_BBB_L_per_h);
failures = compareRows('Table 1', expected, rows(:,3:5));
end

function failures = validateS1(csvDir, texDir)
t = readtable(fullfile(csvDir,'table_s1.csv'),'TextType','string');
rows = dataRows(fullfile(texDir,'phys_params.tex'), 5, 3);
failures = compareRows('Table S1', string(t.DisplayValue), rows(:,3));
end

function failures = validateS2(csvDir, texDir)
t = readtable(fullfile(csvDir,'table_s2.csv'),'TextType','string');
rows = dataRows(fullfile(texDir,'drug_plasma_params.tex'), 8, 3);
expected = string(t{:,5:10});
failures = compareRows('Table S2', expected, rows(:,3:8));
end

function failures = validateS3(csvDir, texDir)
compounds = ["ethanol","diazepam","paracetamol","ibuprofen","indomethacin","mannitol"];
values = strings(0,numel(compounds));
for j = 1:numel(compounds)
    t = readtable(fullfile(csvDir,'table_s3_'+compounds(j)+'.csv'),'TextType','string');
    if j == 1
        values = strings(height(t),numel(compounds));
        parameters = string(t.Parameter);
    else
        assert(isequal(parameters,string(t.Parameter)), ...
            'Table S3 compound CSV files have inconsistent row order.');
    end
    values(:,j) = string(t.DisplayValue);
end
rows = dataRows(fullfile(texDir,'drug_CNS_params.tex'), 8, 3, true);
failures = compareRows('Table S3', values, rows(:,3:8));
end

function failures = validateS4(csvDir, texDir)
t = readtable(fullfile(csvDir,'table_s4.csv'),'TextType','string');
rows = dataRows(fullfile(texDir,'clinical_trials_summary.tex'), 8, 3);
expected = [string(t.DosingRegimen),string(t.NSubjects),string(t.AgeYears),string(t.MalePercent)];
failures = compareRows('Table S4', expected, rows(:,2:5));
end

function failures = validateS5(csvDir, texDir)
t = readtable(fullfile(csvDir,'table_s5.csv'),'TextType','string');
rows = dataRows(fullfile(texDir,'insilico_models.tex'), 7, 2);
expected = [compose('%.2f',t.R2),string(t.N), ...
    string(t.Minimum_logP)+" "+string(t.Maximum_logP), ...
    string(t.Minimum_MW_g_per_mol)+" "+string(t.Maximum_MW_g_per_mol)];
failures = compareRows('Table S5 statistics', expected, rows(:,2:3), rows(:,5:6));

coefficientCounts = [3;4;3;2];
roundedCoefficients = {[-1.24,0.887,-0.5],[-2.06,0.448,-0.366,100], ...
    [round(t.Intercept(3),2),round(t.LogP_Coefficient(3),2),-0.5],[-6.210,0.939]};
for i = 1:height(t)
    actual = numericTokens(rows(i,1));
    base = find(actual == 10,1,'first');
    actual = actual(base+1:min(base+coefficientCounts(i),numel(actual)));
    if ~isequaln(actual(:).', roundedCoefficients{i})
        failures(end+1,1) = "Table S5 equation row "+i+ ...
            ": CSV coefficients "+mat2str(roundedCoefficients{i})+ ...
            "; TeX coefficients "+mat2str(actual(:).'); %#ok<AGROW>
    end
end
end

function rows = dataRows(texFile, columnCount, firstNumericColumn, allowTextValues)
if nargin < 4, allowTextValues = false; end
parts = split(string(fileread(texFile)), "\\");
rows = strings(0,columnCount);
for i = 1:numel(parts)
    cells = split(parts(i), '&').';
    if numel(cells) ~= columnCount, continue; end
    cells = strip(cells);
    if contains(cells(1),'toprule') || contains(join(cells),'thead') || ...
            contains(join(cells),'Value [Ref.]') || contains(join(cells),'Descriptor(s)')
        continue
    end
    probe = cleanCell(cells(firstNumericColumn));
    if isempty(numericTokens(probe)) && ~(allowTextValues && ...
            any(contains(lower(probe),["neutral","acid","--"])))
        continue
    end
    rows(end+1,:) = cells; %#ok<AGROW>
end
end

function failures = compareRows(label, expected, varargin)
actual = [varargin{:}];
failures = strings(0,1);
if ~isequal(size(expected),size(actual))
    failures(end+1,1) = label+": expected "+sizeText(expected)+ ...
        ", found "+sizeText(actual)+" data cells.";
    return
end
for i = 1:size(expected,1)
    for j = 1:size(expected,2)
        e = numericTokens(expected(i,j));
        a = numericTokens(cleanCell(actual(i,j)));
        if ~isequaln(e,a)
            failures(end+1,1) = label+" row "+i+", value column "+j+ ...
                ": CSV '"+expected(i,j)+"'; TeX '"+strip(actual(i,j))+"'."; %#ok<AGROW>
        end
    end
end
end

function text = cleanCell(text)
text = regexprep(string(text),'\\cite\{[^}]*\}','');
text = regexprep(text,'\\textsuperscript\{[^}]*\}','');
text = regexprep(text,'\$\^\{\\text\{[^}]*\}\}\$','');
end

function values = numericTokens(text)
if ismissing(text)
    values = [];
    return
end
text = regexprep(string(text),'-\s+','-');
tokens = regexp(char(text),'(?<![A-Za-z])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?','match');
values = str2double(tokens);
end

function values = displayNumber(values)
values = string(values);
values(ismissing(values) | lower(values)=="nan") = "--";
end

function text = sizeText(value)
text = string(size(value,1))+"x"+string(size(value,2));
end
