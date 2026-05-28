function reflect = test_reflect(tx_alt, rx_alt)
%% Test whether we expect a ray to reflect or not
if tx_alt < 100 && rx_alt < 100
    reflect = true; 
elseif tx_alt > 400 && rx_alt > 400
    reflect = true;
else
    reflect = false;
end