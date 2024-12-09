function output= true_strain_stress(strain,stress,indx_Rp2,indx_Rm)
%Calculate true stress, true strain and hardening exponent values
%Input the tensile tests data

%Calculate True Stress X True Strain
strain_1_Rm=strain(1:indx_Rm); %Strain resume Rm
stress_1_Rm=stress(1:indx_Rm); %Stress resume Rm
resume_Rm_length=length(strain_1_Rm); %Vector length to output values

true_strain=log(1+strain_1_Rm);%True strain until Rm
true_stress=stress_1_Rm; %True stress until Rm

%Calculate true stress
for count=1:length(strain_1_Rm)
true_stress(count)=stress_1_Rm(count)*(1+strain_1_Rm(count));%True stress
end

true_strain_Rp2_Rm=true_strain(indx_Rp2:end); %resume true strain to n calculus
true_stress_Rp2_Rm=true_stress(indx_Rp2:end); %resume true stress to n calculus
resume_Rp2_Rm_length=length(true_strain_Rp2_Rm); %vector length to output values

%Calculate hardening exponent (n) between 2% of plastic deformation and Ag
true_strain_log=log10(true_strain_Rp2_Rm);%log convert
true_stress_log=log10(true_stress_Rp2_Rm);%log convert

%Linear regression---------------------------------------------------
p1=polyfit(true_strain_log,true_stress_log,1); %Line equation
n=p1(1); %Hardening exponent

yfit = p1(1)*true_strain_log+p1(2);
yresid = true_stress_log - yfit;
SQresid = sum(yresid.^2);
SQtotal = (length(true_stress_log)-1)*var(true_stress_log);
R2 = 1 - SQresid/SQtotal;

output=zeros(1500,6);

for i=1:length(true_strain)
   output(i,1)=true_strain(i);
end

for i=1:length(true_stress)
   output(i,2)=true_stress(i);
end

for i=1:length(true_strain_log)
   output(i,3)=true_strain_log(i);
end

for i=1:length(true_stress_log)
   output(i,4)=true_stress_log(i);
end

for i=1:length(true_stress_log)
   output(i,4)=true_stress_log(i);
end

output(1,5)=n;

output(1,6)=p1(2);
output(1,7)=R2;
output(1,8)=resume_Rm_length;
output(1,9)=resume_Rp2_Rm_length;
end

