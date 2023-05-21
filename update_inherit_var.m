function dummy = update_inherit_var(export, dummy, stat, prop)

    if export
        dummy.poro = prop(stat).cell_poro;
        dummy.So = prop(stat).cell_So;
        dummy.Sg = prop(stat).cell_Sg;
        dummy.Sw = prop(stat).cell_Sw;
        dummy.Xi  = prop(stat).cell_Xi;
        dummy.Yi  = prop(stat).cell_Yi;
        dummy.Ki  = prop(stat).cell_Ki;      % For new PVT subroutine  
        dummy.Qai = prop(stat).cell_Qai;     % For adsorption subroutine 
        dummy.GKapp = prop(stat).cell_GKapp; % For gas slippage subroutine
        dummy.GFD = prop(stat).cell_GFD;     % For gas general Fickian diffusion subroutine
        dummy.No = prop(stat).cell_No;
        dummy.Ng = prop(stat).cell_Ng;
        dummy.Vfluid = prop(stat).cell_Vfluid;
        dummy.fv = prop(stat).cell_fv;
        dummy.MDen_o = prop(stat).cell_MDen_o;
        dummy.MDen_g = prop(stat).cell_MDen_g;
        dummy.MDen_w = prop(stat).cell_MDen_w;
        dummy.MolarDensHC = prop(stat).cell_MolarDensHC;
        dummy.MolarDensLiq = prop(stat).cell_MolarDensLiq;
        dummy.MolarDensVap = prop(stat).cell_MolarDensVap;
        dummy.MolarDensWat = prop(stat).cell_MolarDensWat;
        dummy.MolarVolLiq = prop(stat).cell_MolarVolLiq;
        dummy.MolarVolVap = prop(stat).cell_MolarVolVap;
        dummy.kro  = prop(stat).cell_kro;
        dummy.krg  = prop(stat).cell_krg;
        dummy.krw  = prop(stat).cell_krw;
        dummy.Vsic_o = prop(stat).cell_Vsic_o;
        dummy.Vsic_g = prop(stat).cell_Vsic_g;
        dummy.Vsic_w = prop(stat).cell_Vsic_w;
        dummy.pcow  = prop(stat).cell_pcow;
        dummy.pcgo  = prop(stat).cell_pcgo;
        dummy.PV = prop(stat).cell_PV;
        dummy.HCPV  = prop(stat).cell_HCPV;
    else
        return 
    end

end