function [ X_dot ] = Cessna_182_t_los_weight_mobile_target( t , X )


global target_los tau_V tau_psi_1 tau_gamma_1 Cs beta P_max Cl_E S g eta_prop Cd_0 AR e V_t





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

% Definizione della componente x del target mobile

x_t = target_los ( 1 ) + V_t * t      ;
target_los ( 1 ) = x_t                ;

% Definizione del vettore l e degli angoli di visuale gamma e psi 

l       = target_los - [ x , y , z ]  ;
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




X_dot=[ V_dot , psi_dot , gamma_dot , x_dot , y_dot , z_dot , W_dot ]' ;


end