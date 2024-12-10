function equivalent_deformation = plastic_deformation_value(anisotropy_plastic_strain,strain,stress,Ag,Emodul)
%Calculation to equivalence values of plastic deformation on total deformation
%   Detailed explanation goes here

anisotropy_plastic_strain(length(anisotropy_plastic_strain))=Ag; %Insert Value of Ag
for a=1:length(anisotropy_plastic_strain)
    plastic_strain=strain-(stress/(Emodul));
    NumSearch_plastic_strain=(anisotropy_plastic_strain(a)/100);
    indx_plastic_strain=dsearchn(plastic_strain, NumSearch_plastic_strain);
    equivalent_deformation(1,a)=strain(indx_plastic_strain);
end
end

