%Adriano Roberto Krevicz
%Last Update: 07/05/2022 - 17:00

% Matlab code to work the tensile tests data
%Entrar com os dados dos provetes em L0, W0, t0
%Entrar os valores de deformações plásticas de interesse para o gom ma
%matriz anisotropy_plastic_strain
%Entrar o número de extensômetros utilizados no Gom Correlate
%Após verificar os resultados no Gom Correlate, preencher os valores na
%função anisotropy, para que sejam calculados os valores de anisotropia
close all
clear all %#ok<*CLALL>

marks='*sxph+^v<>o*sxph+^v<>o';
color='grbcmkygrbcmky';
line_l={'-','-.',':','-','-.',':','-','-.',':','-','-.',':'};

figuresize=[0 0 19 9*19/16];% FIGURE SIZE IN CENTIMETERS ASPECT RATIO 4:3

disp('##### A ler Ficheiro de input. ######')
FileNames=dir('*.csv');% files with the ansys simulation results

%User set -----------------------------------------------------------------
L0=[50 50 50 50 50 50 50 50 50];
W0=[12.267 12.083 12.083 11.967 11.967 12.000 11.367 11.417 11.383]; %Specimens Width (mm)
t0=[1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20]; %Specimens Thickness (mm)
anisotropy_plastic_strain=[2 2.5 5 7.5 10 12.5 15 17.5 25]; %Values of plastic deformation to Gom Correlate
measured_deformation=length(anisotropy_plastic_strain);
number_extensometers=4; %User set
%--------------------------------------------------------------------------
number_tests=length(FileNames); %Amount of tests


%Calculate Young Modulus and plot Graph------------------------------------
for k=1:number_tests   
[folder, test_name, extension] = fileparts(FileNames(k).name);
names{k}=test_name;
InputData = csvread(FileNames(k).name);

ndata=InputData(1:end,1);
Disp=InputData(1:end,2);
Load=InputData(1:end,3);

disp('##### Setting geometric parameters. ######')
% Input values
A0(k)=W0(k).*t0(k); %mm*mm

strain=Disp./L0(k); %strain values in mm/mm
stress=1000*Load./A0(k); %Stress values in MPa

strain_length(k)=length(strain);
%Create output matrix
if k==1
    strain_out=zeros(1500,number_tests);
    stress_out=zeros(1500,number_tests);
end

for i=1:length(strain)
    strain_out(i,k)=strain(i);
    stress_out(i,k)=stress(i);
end
%% CALCULAR O MODULO DE ELASTICIDADE COM VÁRIOS SEGMENTOS
for j=1:1:5
NS=3+(j-1)*2 %NUM SEGMENTOS
[Fm,indx]=max(Load); %Tensão de rotura à tracção
strain_Fm=strain(indx); % deformação em percentagem existente à Tensão de rotura à tracção
stress_Fm=stress(indx);
Load_seg=[0:Fm/NS:Fm];
stress_seg=Load_seg/A0(k)*1000;


for i=1:NS+1
[value(i), ind(i)] =min(abs((Load(1:indx)-Load_seg(i))));
end

% calculo do modulo de elasticidade em cada segmento ---------------------
for i=1:NS
p=polyfit(strain(ind(i):ind(i+1)), stress(ind(i):ind(i+1)), 1)';  %p is polynomial coefficients in descending powers: p(x)=p(1)*x + p(2)
E_Modulus(i)=p(1)/1000; %Young Modulus in GPa
end
E_Modulus(find(E_Modulus<0))=0;

[val,idx]=min(abs(E_Modulus-200));
E_mod(j)=E_Modulus(idx);

end
[val,idx]=min(abs(E_mod-200));
E_modul(k)=E_mod(idx);   

%Plot figure T/D - E_Tens_def_k---------------------------------------
% plot(strain*100, stress,'LineStyle',line_l{k},'Color',color(k),'LineWidth',1,'DisplayName',test_name);
% line([0 50], [0 50*E_modul(k)*10],'Color','red')
% xlim([0 0.4*100]);
% ylim([0 500]);
% grid on
% box on
% nome=strcat('E_Tens_def_',num2str(k),'.png');
% print(figure(1), '-r1000', '-dpng', nome);
% close
p(2)=0;p(1)=E_modul(k)*1000;

% para determinar a intersecção da recta de proporcionalidade 0,2% com a curva tensão vs Deformação
DefConvi = 0.2;
DefConv=0.2/100;
stress_Proporc02=p(1)*(strain-DefConv)+p(2); %p1 é o módulo de elasticidade e p2 a ordenada na origem
[m(k),indx_Rp02(k)]=min(abs(stress-stress_Proporc02));
strain_Rp02=strain(indx_Rp02(k));  %valor da deformação na intersecção
Rp02=stress(indx_Rp02(k)); %valor de Rp [MPa] para a deformação convencional definida 
strain_Rp02_out=strain_Rp02*100; %para passar a valor percentual
disp(['Tensão Limite Convencional de Proporcionalidade, Rp02= ', num2str(Rp02), ' [MPa]']);
disp(['Extensão total da Tensão Limite Convencional = ', num2str(strain_Rp02_out), ' [%]']);
Rp02_out(k)=Rp02;

% para determinar a intersecção da recta de proporcionalidade 1% com a curva tensão vs Deformação
DefConvi = 1;
DefConv=1/100;
stress_Proporc1=p(1)*(strain-DefConv)+p(2); %p1 é o módulo de elasticidade e p2 a ordenada na origem
[m(k),indx_Rp1(k)]=min(abs(stress-stress_Proporc1));
strain_Rp1=strain(indx_Rp1(k));  %valor da deformação na intersecção
Rp1=stress(indx_Rp1(k)); %valor de Rp [MPa] para a deformação convencional definida 
strain_Rp1_out=strain_Rp1*100; %para passar a valor percentual
disp(['Tensão Limite Convencional de Proporcionalidade, Rp1= ', num2str(Rp1), ' [MPa]']);
disp(['Extensão total da Tensão Limite Convencional = ', num2str(strain_Rp1_out), ' [%]']);
Rp1_out(k)=Rp1;

% para determinar a intersecção da recta de proporcionalidade 2% com a curva tensão vs Deformação
DefConvi = 2;
DefConv=2/100;
stress_Proporc2=p(1)*(strain-DefConv)+p(2); %p1 é o módulo de elasticidade e p2 a ordenada na origem
[m,indx_Rp2(k)]=min(abs(stress-stress_Proporc2));
strain_Rp2=strain(indx_Rp2(k));  %valor da deformação na intersecção
Rp2=stress(indx_Rp2(k)); %valor de Rp [MPa] para a deformação convencional definida 
strain_Rp2_out(k)=strain_Rp2*100; %para passar a valor percentual
disp(['Tensão Limite Convencional de Proporcionalidade, Rp1= ', num2str(Rp2), ' [MPa]']);
disp(['Extensão total da Tensão Limite Convencional = ', num2str(strain_Rp2_out), ' [%]']);
Rp2_out(k)=Rp2;

%%%% Cálculo da Tensão de rotura à tracção ------------------------------
disp('##### Processo de Cálculo da Tensão de rotura à tracção. ######')

[Rm,indx_Rm(k)]=max(stress); %Tensão de rotura à tracção
strain_Rm=strain(indx_Rm(k))*100; % deformação em percentagem existente à Tensão de rotura à tracção
       disp(['Tensão de rotura à tracção, Rm= ', num2str(Rm), '[MPa]']);
       disp(['Deformação total = ', num2str(strain_Rm), '[%]']);
Ag=((strain_Rm/100)-(Rm/(1000*E_modul(k))))*100;
Rm_out(k)=Rm;
Agt_out(k)=strain_Rm;
Ag_out(k)=Ag;

% At Extensão total na rotura [%] -----------------------------------
disp('##### Processo de Cálculo Extensão Total na Rotura. ######')
[At,indx_at(k)]=max(strain);
while stress(indx_at(k))<0
    indx_at(k)=indx_at(k)-1;
end

At=strain(indx_at(k))*100;
disp(['Extensão total na rotura, At= ', num2str(At), '[%]']);
Rf=stress(indx_at(k));
A=((At/100)-(Rf/(1000*E_modul(k))))*100;
At_out(k)=At;
A_out(k)=A;
Rf_out(k)=Rf;

%True Strain x True Stress Calculation------------------------------------
if k==1 %create zeros matrix to write the vectors
true_strain_out=zeros(1500,number_tests);
true_stress_out=zeros(1500,number_tests);
true_strain_log_out=zeros(1500,number_tests);
true_stress_log_out=zeros(1500,number_tests);
end

true_values=true_strain_stress(strain,stress,indx_Rp2(k),indx_Rm(k));

true_strain_out(:,k)=true_values(:,1);
true_stress_out(:,k)=true_values(:,2);
true_strain_log_out(:,k)=true_values(:,3);
true_stress_log_out(:,k)=true_values(:,4);
n_out(k)=true_values(1,5);
P1_out(k)=true_values(1,6);
R2_out(k)=true_values(1,7);
resume_Rm_length(k)=true_values(1,8);
resume_Rp2_Rm_length(k)=true_values(1,9);

linear_equation_out{k}=strcat(num2str(n_out(k),4),'x+',num2str(P1_out(k),4),'(R2=',num2str(R2_out(k),4),')');

%Calculate relation between Plastic Strain and Total Strain (to search values on Gom Correlate)
plastic_deformation(k,:)=plastic_deformation_value(anisotropy_plastic_strain,strain,stress,Ag,E_modul(k)*1000);

close all
hold all
%Plot_Tens_def_-----------------------------------------------------------
plot(strain*100, stress,'LineStyle',line_l{k},'Color',color(k),'LineWidth',1,'DisplayName',test_name);
%Plot Rp02%
plot_Rpvalues(strain*100,stress_Proporc02,strain_Rp02*100,Rp02,0.02)

%Plot Rp1%
plot_Rpvalues(strain*100,stress_Proporc1,strain_Rp1*100,Rp1,1)

%Plot Rp2%
plot_Rpvalues(strain*100,stress_Proporc2,strain_Rp2*100,Rp2,2)

%Plot Rm%
plot_Rpvalues(strain_Rm,Rm,strain_Rm,Rm,3)

xlabel('\epsilon [%]','FontSize',12);
ylabel('\sigma [MPa]','FontSize',12);
xlim([0 1.1*max(strain)*100]);
ylim([0 1.1*max(stress)]);
grid on
box on

nome=strcat('Plot_Tens_def_',num2str(k),'.png');
print(figure(1), '-r1000', '-dpng', nome);
close

end

for i=1:number_tests
names_out(i)=convertCharsToStrings(names{i});
end
%Plot all Engineering stress x strain--------------------------------------
plot_graphics(number_tests,strain_length,strain_out,stress_out,line_l,color,names_out,linear_equation_out,1);

% %Plot All True stress X true strain--------------------------------------
plot_graphics(number_tests,indx_Rm,true_strain_out,true_stress_out,line_l,color,names_out,linear_equation_out,2);

%Plot Line hardening coeficient and Resume Log True Stress x Log Strain----
plot_graphics(number_tests,resume_Rp2_Rm_length,true_strain_log_out,true_stress_log_out,line_l,color,names_out,linear_equation_out,3);

% %Anisotropy evaluate-------------------------------------------------------
% r_out_c=anisotropy(anisotropy_plastic_strain,measured_deformation,number_extensometers,L0,W0,Ag_out,number_tests,names_out);
% close all
% for k=1:number_tests %convert 3D matrix o 2D matrix
%     for i=1:measured_deformation
%         r_out(k,i)=r_out_c(i,1,k);
%     end
% end
% 
% for k=1:number_tests
% r_out_45(1,k)=r_out_c(7,1,k);
% end
% 
%Adjust to ignore Anisotropy calc
%r_out_45=[0 0 0 0 0 0 0 0 0];
% %Mean values---------------------------------------------------------------
%mean_value=mean_value(names_out, L0, t0, W0, A0, E_modul, Rp02_out, Rp1_out, Rp2_out, Rm_out, Rf_out, Ag_out, Agt_out, A_out, At_out,n_out,r_out_45);
mean_value=mean_value(names_out, L0, t0, W0, A0, E_modul, Rp02_out, Rp1_out, Rp2_out, Rm_out, Rf_out, Ag_out, Agt_out, A_out, At_out,n_out);
mean_value_variable={'00','45','90'}';

% Write Equivalent Deformation Excel---------------------------------------
disp('##### Excel data writing. ######')
variable_def={'Test','2%','2,5%','5%','7,5%','10%','12,5%','15%','17.5%','Ag%'};
values_def=[names_out',plastic_deformation*100];
xlswrite('0912_Equivalent-Deformation',variable_def,'out','A1');
xlswrite('0912_Equivalent-Deformation',values_def,'out','A2');

% Write data Excel---------------------------------------------------------
file_name=extractBetween(names_out(1),1,4);
disp('##### Excel data writing. ######')
%variable={'Teste','L0(mm)','t0(mm)','W0(mm)','A0(mm^2)','E_Modulus(GPa)','Rp02(MPa)','Rp1(MPa)','Rp2(MPa)','Rm(MPa)','Rf(MPa)','Ag(%)','Agt(%)','A(%)','At(%)','n(-)','r(xx/15)'};
variable={'Teste','L0(mm)','t0(mm)','W0(mm)','A0(mm^2)','E_Modulus(GPa)','Rp02(MPa)','Rp1(MPa)','Rp2(MPa)','Rm(MPa)','Rf(MPa)','Ag(%)','Agt(%)','A(%)','At(%)','n(-)'};
%values=[names_out', L0', t0', W0', A0', E_modul', Rp02_out', Rp1_out', Rp2_out', Rm_out', Rf_out', Ag_out', Agt_out', A_out', At_out',n_out',r_out_45'];
values=[names_out', L0', t0', W0', A0', E_modul', Rp02_out', Rp1_out', Rp2_out', Rm_out', Rf_out', Ag_out', Agt_out', A_out', At_out',n_out'];
xlswrite(file_name,variable,'out','A1');
xlswrite(file_name,values,'out','A2');
xlswrite(file_name,variable,'media','A1')
xlswrite(file_name,mean_value,'media','A2')
xlswrite(file_name,mean_value_variable,'media','A2')

%xlswrite(file_name,variable_def,'anisotropia','A1')
% xlswrite(file_name,r_out,'anisotropia','B2')
% xlswrite(file_name,mean_value_variable,'anisotropia','A2')