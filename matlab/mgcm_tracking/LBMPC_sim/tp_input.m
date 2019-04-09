function [tp]=tp_input(n_tp,t_bereich,y_bereich,Ts)
%tp_input(n_tp,                     t_bereich,      y_bereich,          Ts)
%tp_input(Zahl der Trainingsdaten,  Zeitbereich,    Amplitudenbereich, Abtastzeit)

tp_save=zeros(n_tp,1);

amplituden=rand(n_tp,1)*(y_bereich(2)-y_bereich(1))+y_bereich(1);
zeiten=round(((rand(n_tp,1))*(t_bereich(2)-t_bereich(1))+t_bereich(1))/Ts);

if t_bereich(1)<0
    display('keine negativen Zeitwerte!');
end;

indexakt=1;

for i=1:n_tp
    amplakt=amplituden(i);
    zeitakt=zeiten(i);

    for j=1:1:zeitakt;

        if indexakt>n_tp
            tp=struct('time',[],'signals',[],'blockname','Trainingsdaten');
            tp.signals=struct('values',tp_save,'dimensions',1,'label','');            
            return;
        end;

        tp_save(indexakt)=amplakt;
        indexakt=indexakt+1;
    end;
end;