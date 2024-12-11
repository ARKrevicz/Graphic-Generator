function all_graphics=plot_graphics(number_tests,vector_length,x_values,y_values,line_l,color,names_out,linear_equation_out,selection)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch selection
    case 1
        strain_length=vector_length;
        strain_out=x_values;
        stress_out=y_values;
        
    for k=1:number_tests
       for i=1:strain_length(k)
            strain_plot(i)=strain_out(i,k);
            stress_plot(i)=stress_out(i,k);
       end
    hold on
    plot(strain_plot*100, stress_plot,'LineStyle',line_l{k},'Color',color(k),'LineWidth',1,'DisplayName',names_out(k));
    clear strain_plot
    clear stress_plot
    end
    xlabel('e [%]','FontSize',12);
    ylabel('\sigma [MPa]','FontSize',12);
    xlim([0 40]);
    ylim([0 600]);
    legend('Location','southwest')
    grid on
    box on
    nome=strcat('Tens_def_all',num2str(k),'.png');
    print(figure(1), '-r1000', '-dpng', nome);

    xlim([0 0.05*100]);
    nome=strcat('Tens_def_Zoom_',num2str(k),'.png');
    print(figure(1), '-r1000', '-dpng', nome);
    close  
    
    case 2
         %Plot All True stress X true strain-------------------------------------------
         indx_Rm=vector_length;
         true_strain_out=x_values;
         true_stress_out=y_values;
                  
        for k=1:number_tests
            true_strain_plot=true_strain_out(1:indx_Rm(k),k);
            true_stress_plot=true_stress_out(1:indx_Rm(k),k);
            %Plot True Stress-Strain-----------------------------------------------
            plot(true_strain_plot*100, true_stress_plot,'LineStyle',line_l{k},'Color',color(k),'LineWidth',1,'DisplayName',names_out(k));
            hold on
        end
        xlabel('\epsilon [%]','FontSize',12);
        ylabel('\sigma _v[MPa]','FontSize',12);
        xlim([0 40]);
        ylim([0 650]);
        legend('Location','southwest')
        grid on
        box on
        nome=strcat('Plot_True_stress_strain_all',num2str(k),'.png');
        print(figure(1), '-r1000', '-dpng', nome);
        close
        
    case 3
   %Plot Line hardening coeficient and Resume Log True Stress x Log Strain----
    true_strain_log_out=x_values;
    true_stress_log_out=y_values;
    resume_Rp2_Rm_length=vector_length;
    
        for n_plot=1:number_tests %#ok<*UNRCH>
           for i=1:resume_Rp2_Rm_length(n_plot)
                true_strain_log_plot(i,1)=true_strain_log_out(i,n_plot);
                true_stress_log_plot(i,1)=true_stress_log_out(i,n_plot);
            end
        %test_name=convertCharsToStrings(names{n_plot});
        plot(true_strain_log_plot,true_stress_log_plot,'LineStyle',line_l{n_plot},'Color',color(n_plot),'LineWidth',1,'DisplayName',names_out(n_plot));
        hold on
        %plot(true_strain_resume_log_plot,polyval(p1,true_strain_resume_log_plot),'LineStyle',':','Color',color(n_plot),'LineWidth',1,'DisplayName',strcat(test_name(n_plot),'-Linear'))
        linear_equation = convertCharsToStrings(linear_equation_out);

        if n_plot==1
            legend_x=1.4*max(true_strain_log_plot); %Text position x
            legend_y=0.5*(max(true_stress_log_plot)-min(true_stress_log_plot))+min(true_stress_log_plot); %Text position y 
        end
        ptext=linear_equation(n_plot);
        text(legend_x,legend_y,{ptext},'Color',color(n_plot),'FontSize',8) 
        hold on
        clear true_strain_log_plot
        clear true_stress_log_plot
        %legend_x=legend_x+0.3;
        legend_y=legend_y-0.01;
        end 
        xlabel('Log \epsilon [mm/mm]','FontSize',12);
        ylabel('Log \sigma_v [MPa]','FontSize',12);
        legend('Location','northwest')
        grid on
        box on

        nome=strcat('Plot_Hardening_Exponent_n_All','.png');
        print(figure(1), '-r1000', '-dpng', nome);
        close all 
        
end

end

