%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:        CNTChargeNeut.m 
% enforces initial charge neutrality

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Em_cnt_old  = -kBT*anti_dummy(Nd'./N1D,delta,Nband);
Em_old      = zeros(bound(13),1);

for ii_row=1:bound(12)
    for ii_col=1:bound(11)
        ii_node=(ii_row-1)*bound(11)+ii_col;
        Em_old(ii_node)=Em_cnt_old(ii_col);
    end
end