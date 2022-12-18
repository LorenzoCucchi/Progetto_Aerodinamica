function [found_af, diff_found] = find_airfoil(naca, alpha, obj)

error=1; % percentage [%]
means=nan(length(alpha),1);
diff=nan(length(alpha),length(naca));
for j=1:length(naca)
    af=naca(j)
    Re=1e6;
    Mach=0;
    [pol,foil] = xfoil(convertStringsToChars(strcat("NACA ",af)),alpha,Re,Mach);
    naca_fit=polyfit(pol.alpha,pol.CL,5);
    obj_fit=polyfit(obj.cl_alfa(:,1),obj.cl_alfa(:,2),5);
    
    
    for i=1:length(alpha)
        if abs(polyval(obj_fit,alpha(i)))>0.1
            diff(i,j)=abs(polyval(naca_fit,alpha(i)) - polyval(obj_fit,alpha(i)))./ abs(polyval(obj_fit,alpha(i))) *100;
        else 
            diff(i,j)=0;
        end
    end
    means(j)=mean(diff(:,j));
    if (all(diff<error))
        found_af=af;
        return 
    end


end

warning('airfoil not found, returned the better one')
[M,index]=min(means);

found_af=naca(index);
diff_found=diff(:,index);
end