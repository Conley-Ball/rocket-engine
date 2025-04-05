function [gamma, cp,dlogV_dlogP_T,dlogV_dlogT_P] = get_thermo(gas)

    dP = pressure(gas)*1e-5;

    n1 = 1/meanMolecularWeight(gas);
    P1 = pressure(gas);
    P2 = pressure(gas) + dP;
    gas2 = Solution('sandiego20161214.yaml');
    set(gas2,'T',temperature(gas),'P',P2,'X',moleFractions(gas));
    equilibrate(gas2,'TP');
    n2 = 1/meanMolecularWeight(gas2);
    dlogn = log(n2)-log(n1);
    dlogP = log(P2)-log(P1);
    
    dlogn_dlogP_T = dlogn/dlogP;
    dlogV_dlogP_T = -1 + dlogn_dlogP_T;
    
    dT = temperature(gas)*1e-5;

    n1 = 1/meanMolecularWeight(gas);
    T1 = temperature(gas);
    T2 = temperature(gas) + dT;
    gas3 = Solution('sandiego20161214.yaml');
    set(gas3,'T',T2,'P',pressure(gas),'X',moleFractions(gas));
    equilibrate(gas3,'TP');
    n2 = 1/meanMolecularWeight(gas3);
    dlogn = log(n2)-log(n1);
    dlogT = log(T2)-log(T1);
    
    dlogn_dlogT_P = dlogn/dlogT;
    dlogV_dlogT_P = 1 + dlogn_dlogT_P;
    
    
    dT = temperature(gas)*1e-5;
    h1 = enthalpy_mass(gas);
    
    cpgas = Solution('sandiego20161214.yaml');
    set(cpgas,'T',temperature(gas)+dT,'P',pressure(gas),'X',moleFractions(gas));
    equilibrate(cpgas,'TP');
    h2 = enthalpy_mass(cpgas);
    cp = (h2-h1)/dT;
    
    cv = cp + 8314.4598/meanMolecularWeight(gas) * dlogV_dlogT_P^2 / dlogV_dlogP_T;
    
    gamma = cp/cv;
    
    gamma = -gamma/dlogV_dlogP_T;

end