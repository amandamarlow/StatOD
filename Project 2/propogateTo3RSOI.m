function [t, r3SOI_N, v3SOI_N, P_3SOI] = propogateTo3RSOI(t_DCO, X_DCO, P_DCO, constants)
    n = size(X_DCO,1);
    [t23RSOI, Xvec_to3RSOI, STM_fromDCO] = integrateTo3RSOI([t_DCO, t_DCO*100], X_DCO, eye(n), constants);
    r3SOI_N = Xvec_to3RSOI(1:3,end);
    v3SOI_N = Xvec_to3RSOI(4:6,end);
    P_3SOI_temp = STM_fromDCO(:,:,end)*P_DCO*STM_fromDCO(:,:,end)';
    P_3SOI = P_3SOI_temp(1:6,1:6);
    t = t23RSOI(end);
end

