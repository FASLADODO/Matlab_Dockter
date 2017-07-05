%[s_cr, P_cr, A] = ColumnBuckling(1,207*10^9,350*10^6,'PinPin','Circle',[0.07],2.5)

%[s_cr,S_allow, P_cr, P_allow, A] = ColumnBuckling(1,203*10^9,689*10^6,'PinPin','Circle',[0.0378],2.5)

[s_cr,S_allow, P_cr, P_allow, A] = ColumnBuckling(20,10.4*10^6,25*10^3,'PinPin','Rectangle',[1,2],4)