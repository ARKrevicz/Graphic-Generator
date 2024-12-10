function graphic_Rvalues = plot_Rpvalues(strain,stress_Proporc,strain_Rp,Rp,selection)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if selection==0.02
    Rp_id='Rp.02=';
    ptext=strcat(Rp_id,num2str(Rp,5), ' MPa');
    points=[' \leftarrow ' ptext];
elseif selection==1
    Rp_id='Rp1.0=';
    ptext=strcat(Rp_id,num2str(Rp,5), ' MPa');
    points=[' \leftarrow ' ptext];    
elseif selection==2
    Rp_id='Rp2.0=';
    ptext=strcat(Rp_id,num2str(Rp,5), ' MPa');
    points=[' \leftarrow ' ptext];   
else
    Rp_id='Rm=';
    clear points
    ptext=strcat(Rp_id,num2str(Rp,5), ' MPa');
    points={ptext};
end

plot(strain, stress_Proporc);
plot(strain_Rp,Rp,'r.','MarkerSize',20)
%ptext=strcat(Rp_id,num2str(Rp,5), ' MPa');
%text(strain_Rp,Rp,[' \leftarrow ' ptext],'FontSize',12)
if selection==3
    text(strain_Rp,1.05*Rp,points,'FontSize',12) 
else
    text(strain_Rp,Rp,points,'FontSize',12) 
end

end

