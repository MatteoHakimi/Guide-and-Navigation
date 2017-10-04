function [ X_dot ] = Cessna_182_t_los_weight_hexagon( t , X )


global W_P tau_V tau_psi_1 tau_gamma_1 Cs beta P_max Cl_E S g eta_prop Cd_0 AR e i_count 





% Definizione  dello spazio di stato

V     = X ( 1 ) ;
psi   = X ( 2 ) ;
gamma = X ( 3 ) ; 
x     = X ( 4 ) ;
y     = X ( 5 ) ;
z     = X ( 6 ) ;
W     = X ( 7 ) ;

% Definizione di funzioni utili e di velocità richiesta

[ ~ , ~ , ~ , rho ] = atmosisa( z )                           ;
P = @(z, delta_t) P_max * exp( -beta * z ) * delta_t          ;
V_star = sqrt ( 2 * W / ( rho * S * Cl_E ) )                  ;

% Definizione del vettore l e degli angoli di visuale gamma e psi

l       = W_P ( i_count , : ) - [ x , y , z ]  ;
psi_r   = atan2 ( l ( 2 ) , l ( 1 ) ) ;
gamma_r = asin( l ( 3 ) / norm( l ) ) ;



% Sistema di equazioni differenziali

V_dot     = ( V_star - V ) / ( tau_V )            ;
psi_dot   = ( psi_r - psi ) / ( tau_psi_1 )       ;
gamma_dot = ( gamma_r - gamma ) / ( tau_gamma_1 ) ;
    
x_dot     = V * cos( gamma ) * cos( psi ) ;
y_dot     = V * cos( gamma ) * sin( psi ) ;
z_dot     = V * sin( gamma )              ;





n  = 1 / g * sqrt ( ( V * psi_dot * cos ( gamma ) ) ^2 + ( V_dot * gamma_dot + g * cos ( gamma ) )^2 ) ;



Cl = 2 * W * n / ( S * rho * V ^2 )            ;
Cd = Cd_0 + Cl ^ 2 / (pi * AR *e)              ;
D  = 1 / 2 * rho * Cd * V ^2 * S               ;
T  = D + W * ( sin ( gamma ) + 1 / g * V_dot ) ;  
                                                                   


delta_t = T * V / ( eta_prop * P ( z , 1 ) )                           ;
W_dot     = - Cs * P ( z ,delta_t )                                    ;

% Condizione di soglia per switch del way point fissata a 200 m

if norm ( l ) < 200 && i_count < size ( W_P , 1 ) 
    i_count  = i_count + 1 ;
end


X_dot=[ V_dot , psi_dot , gamma_dot , x_dot , y_dot , z_dot , W_dot ]' ;


end