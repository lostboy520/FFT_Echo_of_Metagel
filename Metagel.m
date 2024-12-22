clear
addpath(genpath(pwd));
[filename,filepath]=uigetfile('*.bin','OPEN','MultiSelect', 'on');

    address=[filepath,filename];
    fid=fopen(address,'r');                     
    a=fread(fid,'ubit14',2);                   
    fclose(fid);                                
    sig=(10/16384)*double(a(51757:end))-5;     
    clear a; 
fs=2e7;
fpass=[1.3e6 3e6];
y=fir_butter(sig,fpass,fs,0);
plot(y(1:1000));
clear locrec
kk=1;
jj=1;
mm=150; 
nn=1;
oo=1;
while mm<length(y)
     if y(mm)>0.3
           locrec(nn)=mm;
            mm=mm+2e4-40;     
            nn=nn+1;
            oo=1;
         elseif oo>400
              mm=mm+2e4-400;
              oo=1;
         else
             mm=mm+1; 
             oo=oo+1;
     end
end
num_rec=length(locrec); 
fs=2*10^7;  
Ts=1/fs;
N0=2*10^5;
L0=100000;
ya=(zeros(1,L0))';
aa=1;
bb=1;
empty=[];
window=4;
d1=20;
d2=150;
freqstart1=0.8e4; freqstop1=4e4;
while window*(aa+2)<num_rec        
     siga=([]);
     for m0=window*aa-1:1:window*(aa+1)+1      
     sigm=y(locrec(m0)-d1:locrec(m0)+d2);
     sigm=[sigm;ya];    
     siga=[siga;sigm];  
     end 
    sigrec=siga';
    clear siga
N2=length(sigrec);
spectral_resolution=fs/N2;
freqstart=fix(1.6e6/spectral_resolution); freqstop=fix(3e6/spectral_resolution); 
    n2=0:N2-1;
    f2=n2*fs/N2;
    Xrec=fft(sigrec,N2);
    ampl2=abs(Xrec);
    frec=f2;
    amplrec=ampl2;
    [ecoup,ecodown] = envelope(amplrec(freqstart:freqstop),70,'peak');
    recup(aa,:)=ecoup;
    clear ecoup， Xrec
   [p,n]=findpeaks(recup(aa,freqstart1:freqstop1),'MinPeakHeight',1);   
   [~,num_p]=max(p); 
    l=n(num_p);
    TF=isempty(l);
    if TF==0
    dia_u1(bb,:)=l+freqstart+freqstart1;  
    freqstop1=freqstart1+l+1000;
    freqstart2=freqstart1+l-1000;
    freqstart1=freqstart2;
    else
    bb=bb-1;    
    empty=[empty,aa];
    end
   aa=aa+1;         
   bb=bb+1;
   clear l，freqstart2，sigrec，amplrec，ampl2
   clear m
   clear n
end
t0=4/1000;  
reflec1=dia_u1*fs/N2;  
REF=reflec1(1:end,1);
T=t0:t0:t0*length(REF);
figure();plot(T,REF,'r');
xlabel('Time'); ylabel('freq(Hz)');
