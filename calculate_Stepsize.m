function alpha = calculate_Stepsize(Df, Df_past, x, x_past)

    
    alpha = ( (x - x_past)' * (Df - Df_past)' ) ./ ((Df-Df_past)*(Df-Df_past)');
    

end 