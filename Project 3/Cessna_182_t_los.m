function [ X_dot ] = Cessna_182_t_los( t , X )


global target_los V_star tau_V tau_psi tau_gamma


% Definizione  dello spazio di stato

V     = X ( 1 ) ;
psi   = X ( 2 ) ;
gamma = X ( 3 ) ; 
x     = X ( 4 ) ;
y     = X ( 5 ) ;
z     = X ( 6 ) ;

% Definizione del vettore l e degli angoli di visuale gamma e psi

l       = target_los - [ x , y , z ]  ;
psi_r   = atan2 ( l ( 2 ) , l ( 1 ) ) ;
gamma_r = asin( l ( 3 ) / norm( l ) ) ;


% Sistema di equazioni differenziali

V_dot     = ( V_star - V ) / ( tau_V )          ;
psi_dot   = ( psi_r - psi ) / ( tau_psi )       ;
gamma_dot = ( gamma_r - gamma ) / ( tau_gamma ) ;
    
x_dot     = V * cos( gamma ) * cos( psi ) ;
y_dot     = V * cos( gamma ) * sin( psi ) ;
z_dot     = V * sin( gamma )              ; 


X_dot=[ V_dot , psi_dot , gamma_dot , x_dot , y_dot , z_dot ]' ;


end

