function [dummy] = alloc_inherit_var()

    dummy.poro = [];
    dummy.So = [];
    dummy.Sg = [];
    dummy.Sw = [];
    dummy.Xi  = [];
    dummy.Yi  = [];
	dummy.Ki  = [];  % For new PVT SUBROUTINE
    dummy.Qai = [];  % For adsorption SUBROUTINE 
    dummy.GFD = [];  % For general Fickian SUBROUTINE 
    dummy.No = [];
    dummy.Ng = [];
    dummy.Vfluid = [];
    dummy.fv = [];
    dummy.MDen_o = [];
    dummy.MDen_g = [];
    dummy.MDen_w = [];
    dummy.MolarDensHC = [];
    dummy.MolarDensLiq = [];
    dummy.MolarDensVap = [];
    dummy.MolarDensWat = [];
    dummy.MolarVolLiq = [];
    dummy.MolarVolVap = [];
    dummy.kro  = [];
    dummy.krg  = [];
    dummy.krw  = [];
    dummy.Vsic_o = [];
    dummy.Vsic_g = [];
    dummy.Vsic_w = [];
    dummy.pcow  = [];
    dummy.pcgo  = [];
    dummy.PV = [];
    dummy.HCPV  = [];
    dummy.GKapp = []; % For gas phase slippage flux

end