%synaptome model
clear
clf

%uncomment the alternative you want to use
use_pattern=1; %theta burst
%use_pattern=2; %gamma
%use_pattern=3; %theta
%use_pattern=4; %high gamma pairs

%uncomment the alternative you want to use
%age=1;  %1w
age=3;  %3M
%age=18; %18M

%uncomment the alternative you want to use
save_figure=0; %do not save the figure(s)
%save_figure=1; %save the figure(s) as png and pdf

%uncomment the alternative you want to use
threeDplot=0; %do not show extra 3D hairpin plot
%threeDplot=1; %show extra 3D hairpin plot

%nothing below needs to be changed, but welcome to experiment yourself

%%%load synapse data
gradientData = load ('gradientData.txt'); %gradient data for the synapses
xyCoord = load ('xyCoord.txt'); %x-y coordinates for the synapses

%%%initialize some parameters
v=[0.25  130    5.4   13    0.45  860]; %see Zhu et al., 2018 for furhter info on these values which were set using 10 epochs theta burst data to scale EPSP amplitudes and time constants
epspThresh=2.5;
plotthresh=15;
%totTime = 1450; %stimulate with only 8 spikes
totTime = 3850; %for 20 spike stimulation
t = linspace(0, totTime, totTime);
maxsumpeaks=0;
minsumpeaks=10000; %probably big enough

%%%set up stimulation pattern
%20 spikes
st10=[25  50   75  100  225  250  275  300  425  450  475  500  625  650  675  700  825  850  875  900]; %theta burst
st20=[25   50   75  100  125  150  175  200  225  250  275  300  325  350  375  400  425  450  475  500]; %gamma
st30=[25   225   425   625   825  1025  1225  1425  1625  1825  2025  2225  2425  2625  2825  3025  3225  3425  3625  3825]; %theta
st40=[25    38   338   351   651   664   964   977  1277  1290  1590  1603  1903  1916  2216  2229  2529  2542  2842  2855]; %high gamma pairs

switch use_pattern
case 1
pattmat=[st10]; %theta burst
case 2
pattmat=[st20]; %gamma
case 3
pattmat=[st30]; %theta
case 4
pattmat=[st40]; %high gamma pairs
end

%%% modify gradient data to account for age differences
switch age
case 1 %1W
redfac1=2.1;
redfac2=5.4;
redfac3=3.0;
redfac4=6.6;
case 3 %3M
redfac1=1.6;
redfac2=2.5;
redfac3=2.1;
redfac4=3.7;
case 18 %18M
redfac1=1.0;
redfac2=1.0;
redfac3=1.0;
redfac4=1.0;
end

minv=min(gradientData(:,1));
maxv=max(gradientData(:,1));
diffshift=(maxv-minv)*(1-1/redfac1);
gradientData(:,1)=((gradientData(:,1)-minv)/redfac1)+minv+diffshift;

minv=min(gradientData(:,2));
maxv=max(gradientData(:,2));
diffshift=(maxv-minv)*(1-1/redfac2);
gradientData(:,2)=((gradientData(:,2)-minv)/redfac2)+minv+diffshift;

minv=min(gradientData(:,3));
maxv=max(gradientData(:,3));
diffshift=(maxv-minv)*(1-1/redfac3);
gradientData(:,3)=((gradientData(:,3)-minv)/redfac3)+minv+diffshift;

minv=min(gradientData(:,4));
maxv=max(gradientData(:,4));
diffshift=(maxv-minv)*(1-1/redfac4);
gradientData(:,4)=((gradientData(:,4)-minv)/redfac4)+minv+diffshift;

counts = zeros(length(pattmat(:,1)), length(xyCoord));
%%% iterate over all inputs and over all synapses
for patt=1:length(pattmat(:,1))
spikeTimes=pattmat(patt,:);
for all=1:length(xyCoord(:,1))

%compute EPSP shape
%depression
Ad=v(1)*gradientData(all,1);
tauD=v(2)*gradientData(all,2);
tD = linspace(0, 5*tauD, 5*tauD+1);
depr = Ad * exp(-tD/tauD);

%fast facilitation
fSaturation = 3.3;
Af=v(3)*gradientData(all,3);
tauF=v(4)*gradientData(all,4);
tF = linspace(0, 5*tauF, 5*tauF+1);
fac = Af * exp(-tF/tauF);

%slow facilitation
Aff=v(5);
tauFF=v(6);
tFF = linspace(0, 5*tauFF, 5*tauFF+1);
sfac = Aff * exp(-tFF/tauFF);

%the basic EPSP waveform
Ae = 1;
tau1 = 3.0;
tau2 = 0.4;
tE = linspace(0, tau1*10, tau1*10+1);
psp = Ae * (exp(-tE/tau1) - exp(-tE/tau2)); %tau1>tau2

totEpsp = zeros(1, totTime);
totE = zeros(1, totTime);

%first spike, no fac or depr (as no preceeding spikes)
if spikeTimes(1) < totTime
preFill = zeros(1, spikeTimes(1));
if length(preFill)+length(psp) < totTime
postFill = zeros(1, totTime-length(preFill)-length(psp));
epspTrace = [preFill psp postFill];
else
epspTrace = [preFill psp(1:totTime-length(preFill))];
end %if 
end %if
totEpsp=totEpsp+epspTrace;
totE=totE+epspTrace;

%second to last spike
for j=2:length(spikeTimes);
totDepr = ones(1, totTime);
totFac = ones(1, totTime);
totsFac = ones(1, totTime);
for i=1:j-1;
if spikeTimes(i) < totTime
preFill = zeros(1, spikeTimes(i));

%depression time profile
if length(preFill)+length(depr) < totTime
postFillDepr = zeros(1, totTime-length(preFill)-length(depr));
depression = [preFill depr postFillDepr];
else
depression = [preFill depr(1:totTime-length(preFill))];
end %if 

%facilitation time profile
if length(preFill)+length(fac) < totTime
postFillFac = zeros(1, totTime-length(preFill)-length(fac));
facilitation = [preFill fac postFillFac];
else
facilitation = [preFill fac(1:totTime-length(preFill))];
end %if 

%slow facilitaiton time profile
if length(preFill)+length(sfac) < totTime
postFillsFac = zeros(1, totTime-length(preFill)-length(sfac));
sfacilitation = [preFill sfac postFillsFac];
else
sfacilitation = [preFill sfac(1:totTime-length(preFill))];
end %if 

end %if spikeTimes(i) < totTime
%for each synapse, all input-induced traces are summed (superimposed)
totDepr=totDepr-depression;
totFac=totFac+facilitation+sfacilitation;  %should we add sfacilitation here?
end %for i=1:j-1;
totDepr = max(totDepr,0.0);
totFac = min(totFac,fSaturation);

Ae = totDepr(spikeTimes(j))*totFac(spikeTimes(j));
psp = Ae * (exp(-tE/tau1) - exp(-tE/tau2)); %tau1>tau2
if spikeTimes(j) < totTime
preFill = zeros(1, spikeTimes(j));
if length(preFill)+length(psp) < totTime
postFill = zeros(1, totTime-length(preFill)-length(psp));
epspTrace = [preFill psp postFill];
else
epspTrace = [preFill psp(1:totTime-length(preFill))];
end %if 
end %if
totEpsp=totEpsp+epspTrace;
totE=totE+(epspTrace);
end %for j=2:length(spikeTimes);

%normalize to amplitude of first EPSP
indm=find (t>20 & t<30);
peakE1=max(totEpsp(indm));
totEpsp=totEpsp/peakE1;

%%%find and sum up peak amplitudes of the EPSPs
searchinterv=[];
for i=1:length(spikeTimes)
searchinterv=[searchinterv; spikeTimes(i)-5 spikeTimes(i)+10];
end

peakVek=[];
for i=1:length(spikeTimes)
peakVek=[peakVek; max(totEpsp(searchinterv(i,1):searchinterv(i,2)))];
end

sumpeaks=sum(peakVek);
if sumpeaks > maxsumpeaks
maxsumpeaks = sumpeaks;
end
if sumpeaks < minsumpeaks
minsumpeaks = sumpeaks;
end

%set color scale, all simulations use the same color scale
eps=0.0001; %min and max below are truncated values, eps for col to stay strictly in [0,1]
  %run the script once for each case of age*pattern to find out the min and max, printing the minsumpeaks and maxsumpeaks below, use the smallest (largest) of all cases
minsumpeaksreg=8.1027-eps;
maxsumpeaksreg=62.1981+eps;
rangepeaks=maxsumpeaksreg-minsumpeaksreg;
counts(patt, all) = (sumpeaks-minsumpeaksreg)/rangepeaks;

end %for all %all spikes
end %for patt %all patterns

%minsumpeaks %uncomment to print value, update minsumpeaksreg if necessary
%maxsumpeaks %uncomment to print value, update maxsumpeaksreg if necessary

%%%plot the result
%define z values for plot
zMat=[];
for i=1:11
zMat = [zMat; counts(1,(i-1)*11+1:(i-1)*11+11)]; %for plot3, rows are y, cols are x 
end
figure(1)
imagesc(zMat);
colormap(fake_parula());
maxLim = 1.0;
minLim = 0.0;
hAx = gca;
hAx.CLim = [minLim maxLim];
if save_figure
print -dpng 2Dplot.png
print -dpdf 2Dplot.pdf
end

%3dplot
if threeDplot
%define x-values for plot
x = (1:1:11).';
xMat = repmat(x, 1, 11);
%define y values for plot
y = 1:1:11;
yMat = repmat(y, numel(x), 1);

figure(2)
for i=1:11
for j=1:11
clr=[counts(1,(i-1)*11+j) 0.0 1-counts(1,(i-1)*11+j)];
stem3(xMat(i,j), yMat(i,j), zMat(i,j), 'o', 'Color', clr), hold on
end
end
grid;
view(70, 30);
axis([0 12 0 12 0 1]);
% Add title and axis labels
zlabel('Normalized amplitude');
xlabel('Tangential direction');
ylabel('Radial direction');
if save_figure
print -dpng 3Dplot.png
print -dpdf 3Dplot.pdf
end
end %3d plot
