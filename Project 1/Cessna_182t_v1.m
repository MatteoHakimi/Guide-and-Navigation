function [ X_dot ] = Cessna_182t_v1 ( t ,X )


global P_max eta_prop delta_t delta_h Cl Cd S g phi C_s


% Definizione  dello spazio di stato

V     = X ( 1 ) ;
psi   = X ( 2 ) ;
gamma = X ( 3 ) ; 
% x     = X ( 4 ) ;
% y     = X ( 5 ) ;
z     = X ( 6 ) ;
W     = X ( 7 ) ;


% Calcolo di termini noti

[~, ~, ~, rho] = atmosisa( z )       ;  
P       = P_max * delta_t * delta_h  ;
T       = eta_prop * P / V           ;
L       = 1/2 * rho * V ^ 2 * S * Cl ;
D       = 1/2 * rho * V ^ 2 * S * Cd ;
m       = W / g                      ;


% Sistema di equazioni differenziali

V_dot      = 1 / m * ( T - D - W * sin ( gamma ) )                 ;
psi_dot    = 1 /( m * V * cos ( gamma ) ) * L * sin ( phi )        ;
gamma_dot  = 1 /( m * V ) * ( L * cos( phi ) -  W * cos( gamma ) ) ;

x_dot      = V * cos( gamma ) * cos( psi ) ;
y_dot      = V * cos( gamma ) * sin( psi ) ;
z_dot      = V * sin( gamma )              ; 
W_dot      = - C_s * P                     ;


X_dot=[ V_dot , psi_dot , gamma_dot , x_dot , y_dot , z_dot , W_dot ]' ;

end