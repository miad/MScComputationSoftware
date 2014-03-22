%close all;
clf;
%arr = {'test_output'};
potDepths = [-5:0.05:-0.05 0.05:0.05:4.95];
arr=cell(size(potDepths, 2), 1);

nmax = [30];
cl = {'bx', 'go', 'm+'};

for npr=1:size(nmax, 2)
    for i=1:size(potDepths, 2)
        p = potDepths(i);
        chr = '';
        if abs(p) < 1E-6
            p = 0;
            chr='';
        end
        arr{i, 1} = sprintf('/net/data2/riklund/BUSCH/output_%s%3.2f', chr,p);
    end
    PTS=PlotEnergies(arr, potDepths,0, [cl{npr}]);
end
hold on;
MinusOneOverG = -5:0.1:4.9;


redo = 0;
mass=2.2;
hbar=0.303;
omega=400.4;
mu=mass/2;
aref=sqrt(hbar/(mu*omega));

if ~exist('RedPoints','var') || redo ~= 0
    NoRedPoints = 16;
    RedPoints = zeros(1,1:NoRedPoints*size(MinusOneOverG,2));
    gPoints = RedPoints;
    count = 1;

    for startvalue = 0:15
        en=zeros(size(MinusOneOverG));
        for i=1:size(en, 2)
           RedPoints(count) = 1/2+  fsolve(@(e) gamma(-e/2+1/4)/gamma(-e/2+3/4)-2*MinusOneOverG(i), (startvalue-0.9));
           %RedPoints(count) = 1 +  fsolve(@(e) gamma(-e/2)/gamma(-e/2+1/2)*1/2*mu/hbar^2-aref*MinusOneOverG(i), (startvalue-0.9));
           gPoints(count) = MinusOneOverG(i);
           count=count+1;
        end
    end
end
plot(gPoints, RedPoints,'r.');
ylabel('Energy /(hbar*omega)');
xlabel('Coupling coefficient /(aref*hbar*omega)');