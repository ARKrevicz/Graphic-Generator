function [r_out] = anisotropy(anisotropy_plastic_strain,measured_deformation,number_extensometers,L0,width,Ag_out,number_tests,names_out)
% Anisotropy Calculation
% Entry the detailed matrice for the deformation values from each specimen
% test
%The values are defined by gom correlate analysis
%Ext 1 is the length extensometer

%Anisotropy Calculation--------------------------------------------------
marks='*sxph+^v<>o';
color='grbcmky';
line_l={'-','-.',':','-','-.',':'};
anisotropia=zeros(measured_deformation,number_extensometers,number_tests);
% Matrix:  %Elongation x Result extensometers from Gom Correlate
%Test 1-------------------------------------------------------------------
%Ext1  
delta_L_1=[1.076,1.351,2.598,3.860,5.119,6.357,7.604,8.904,9.815]';
                   
%           Ext2   Ext3   Ext4
delta_w_1=[-0.097 -0.131 -0.129
           -0.129 -0.167 -0.164
           -0.269 -0.302 -0.294
           -0.416 -0.442 -0.427
           -0.542 -0.582 -0.581
           -0.671 -0.705 -0.698
           -0.772 -0.848 -0.818
           -0.875 -0.974 -0.932
           -0.937 -1.071 -1.030];
%Test 2-------------------------------------------------------------------
%Ext1  
delta_L_2=[1.084,1.354,2.596,3.852,5.107,6.366,7.607,8.895,10.130]';
                   
%           Ext2   Ext3   Ext4
delta_w_2=[-0.114 -0.123 -0.134
           -0.142 -0.159 -0.166
           -0.279 -0.305 -0.303
           -0.402 -0.450 -0.455
           -0.512 -0.579 -0.585
           -0.632 -0.721 -0.710
           -0.731 -0.834 -0.842
           -0.829 -0.955 -0.963
           -0.906 -1.085 -1.077];
       
%Test 3-------------------------------------------------------------------
%Ext1  
delta_L_3=[1.093,1.334,2.605,3.842,5.099,6.368,7.585,8.847,10.004]';
                   
%           Ext2   Ext3   Ext4
delta_w_3=[-0.129 -0.089 -0.129
           -0.157 -0.116 -0.151
           -0.313 -0.251 -0.308
           -0.450 -0.390 -0.439
           -0.584 -0.518 -0.567
           -0.728 -0.628 -0.702
           -0.846 -0.752 -0.816
           -0.966 -0.856 -0.928
           -1.078 -0.948 -1.030];
%-------------------------------------------------------------------------
    delta_L_global(:,:,1)=delta_L_1;
    delta_L_global(:,:,2)=delta_L_2;
    delta_L_global(:,:,3)=delta_L_3;

    delta_w_global(:,:,1)=delta_w_1;
    delta_w_global(:,:,2)=delta_w_2;
    delta_w_global(:,:,3)=delta_w_3;
%--------------------------------------------------------------------------

for k=1:number_tests
    for i=1:measured_deformation
        for j=1:(number_extensometers-1)
            w_final_global(i,j,k)=width(k)+delta_w_global(i,j,k);
        end
    end
    wm_global=mean(w_final_global,2);

    %Calculte Final Dimensions
    for i=1:measured_deformation
        l_final_global(i,1,k)=L0(k)+delta_L_global(i,1,k);
    end          

    for i=1:measured_deformation
        anisotropy_global(i,1,k)=(log(wm_global(i,1,k)/width(k)))/(log((width(k)*L0(k))/(wm_global(i,1,k)*l_final_global(i,1,k))));
        strain_global(i,1,k)=((l_final_global(i,1,k)-L0(k))/L0(k))*100;
    end
    %anisotropy_out=anisotropy_global;
end

for k=1:number_tests
    
    anisotropy_plastic_strain(end)=Ag_out(k);
    test_name=names_out(k);
    plot(anisotropy_plastic_strain,anisotropy_global(:,1,k),'LineStyle',line_l{k},'Color',color(k),'LineWidth',1,'DisplayName',strcat(test_name,'-Anisotropy'))
    hold on
    %Data for 15%strain
    %r_out(k)=anisotropy_global(7,1,k);
    r_out=anisotropy_global;
end
xlabel('\epsilon [%]','FontSize',12);
ylabel('r','FontSize',12);
xlim([0 25]);
ylim([0 1.2]);
legend('Location','southeast')
grid on
box on
nome=strcat('Anisotropy_all','.png');
print(figure(1), '-r1000', '-dpng', nome);

close all
end

