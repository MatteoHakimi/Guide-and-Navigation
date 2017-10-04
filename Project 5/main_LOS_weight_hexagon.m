tic
close all 
clc, clear


global W_P tau_V tau_psi_1 tau_gamma_1 Cs beta P_max Cl_E S g eta_prop Cd_0 AR e i_count 

% Parametri velivolo

b        = 11                          ;
S        = 16.2                        ;
AR       = b ^ 2 / S                   ;
P_max    = 172e3                       ;
Cd_0     = 0.022                       ;
e        = 0.8                         ;
eta_prop = 0.82                        ;
beta     = 1e-4                        ;
g        = 9.81                        ;
Cs       = 5e-7                        ;
Cl_E     = sqrt ( Cd_0 * pi * AR * e ) ;

% Definizioni funzioni utili

d2r = @(x) x * pi / 180                              ;
r2d = @(x) x / pi * 180                              ;
P = @(z, delta_t) P_max * exp( - beta * z) * delta_t ;


% Velocità di stallo e peso del velivolo

W                = 14100    ;
V_s              = 91 / 3.6 ;


% Condizioni iniziali del velivolo

V_0      = 1.2 * V_s    ;
psi_0    = d2r ( - 45 ) ;
gamma_0  = 0            ;
x_0      = 0            ;
y_0      = 0            ;
z_0      = 0            ;   
W_0      = W            ;


% Posizione target

R = 2e4 ; 

W_P =zeros ( 6 ,3 ) ;

W_P ( 1 , : ) = [ R * cos( d2r ( 30 ) ) , R * sin( d2r ( 30 ) ) , 1.5e3 ]     ;
W_P ( 2 , : ) = [ 0 , R , 1.5e3 ]                                             ;
W_P ( 3 , : ) = [ - R * cos( d2r ( 30 ) ) ,  R * sin( d2r ( 30 ) ) , 1.5e3 ]  ;
W_P ( 4 , : ) = [ - R * cos( d2r ( 30 ) ) , - R * sin( d2r ( 30 ) ) , 1.5e3 ] ;
W_P ( 5 , : ) = [ 0  , -R , 1.5e3 ]                                           ;
W_P ( 6 , : ) = [ R * cos( d2r ( 30 ) ) , - R * sin( d2r ( 30 ) ) , 1.5e3 ]   ;


% Tempi tau


tau_V     = 30 ;
tau_psi   = 20 ;
tau_gamma = 20 ;

% Lunghezza vettori

n_length = length ( tau_psi )      ;
m_length = length ( tau_gamma )    ;

% Integrazione delle equazioni con Runge Kutta 4° ordine
% Definizione del tempo di integrazione e delle condizioni iniziali

dt               = 0.1                                               ;
t_span           = 0 : dt : 2340                                     ;
X_0              = [ V_0 , psi_0 , gamma_0 , x_0 , y_0 , z_0 , W_0]  ;

% Definizione cell

V_c       = cell ( n_length , m_length ) ;
psi_c     = cell ( n_length , m_length ) ;
gamma_c   = cell ( n_length , m_length ) ;
x_c       = cell ( n_length , m_length ) ;
y_c       = cell ( n_length , m_length ) ;
z_c       = cell ( n_length , m_length ) ;
W_c       = cell ( n_length , m_length ) ;
phi_c     = cell ( n_length , m_length ) ;
Cl_c      = cell ( n_length , m_length ) ;
delta_t_c = cell ( n_length , m_length ) ;
n_c       = cell ( n_length , m_length ) ;
t_c       = cell ( n_length , m_length ) ;

i_count = 1     ; 

for i = 1 : n_length     ; 
    
    for j = 1 : m_length ;


tau_psi_1 = tau_psi ( i )     ;
tau_gamma_1 = tau_gamma ( j ) ;


[t,X]            = ode45( 'Cessna_182_t_los_weight_hexagon', t_span , X_0 ) ;

V     = X ( : , 1 ) ;
psi   = X ( : , 2 ) ;
gamma = X ( : , 3 ) ; 
x     = X ( : , 4 ) ;
y     = X ( : , 5 ) ;
z     = X ( : , 6 ) ;


V_dot     = gradient ( V , t )     ;
psi_dot   = gradient ( psi , t )   ;
gamma_dot = gradient ( gamma , t ) ; 
x_dot     = gradient ( x , t )     ; 
y_dot     = gradient ( y , t )     ;
z_dot     = gradient ( z , t )     ;

gamma  = asin ( z_dot ./ V )                                                                                    ;
psi    = atan2 ( y_dot , x_dot )                                                                                ;
psi    = r2d ( psi )                                                                                            ;
phi    = atan2 ( V .* cos ( gamma ) .* psi_dot , V .* gamma_dot + g * cos ( gamma ) )                           ; 
phi    = r2d ( phi )                                                                                            ;
n      = 1 / g * sqrt ( ( V .* psi_dot .* cos ( gamma ) ) .^2 + ( V_dot .* gamma_dot + g * cos ( gamma ) ).^2 ) ;

[ ~ , ~ , ~ , rho ] = atmosisa( z ) ;
rho                 = rho'          ;

Cl = 2 * W * n ./ ( S * rho .* V .^2 )         ;
Cd = Cd_0 + Cl .^ 2 / (pi * AR *e)             ;
D  = 1 / 2 * rho .* Cd .* V .^2 * S            ;
T  = D + W * ( sin ( gamma ) + 1 / g * V_dot ) ;  

gamma    = r2d ( gamma )                       ;                                                                     

delta_h = exp( -beta * z )                     ;
delta_t = T .* V ./ ( eta_prop * P ( z , 1 ) ) ;

V_c     { i , j }     = X ( : , 1 )                  ;
psi_c   { i , j }     = psi                          ;
gamma_c { i , j }     = gamma                        ; 
x_c     { i , j }     = X ( : , 4 )                  ;
y_c     { i , j }     = X ( : , 5 )                  ;
z_c     { i , j }     = X ( : , 6 )                  ;
W_c     { i , j }     = X ( : , 7 ) / W_0 *100 - 100 ;

phi_c     { i , j } = phi                                    ;
delta_t_c { i , j } = T .* V ./ ( eta_prop * P ( z , 1 ) )   ;
Cl_c      { i , j } = Cl                                     ;
n_c       { i , j } = n                                      ;
t_c       { i , j } = t                                      ;


    end
end


for i = 1 : n_length ;
    for j = 1 : m_length ;
        
        if max ( delta_t_c { i , j } )  <= 1 && min ( delta_t_c { i , j } )  > 0
            
        fprintf( 'La manovra con tau_psi = %4.1f e tau_gamma = %4.1f e'' ammissibile con Delta_Max W / W_0 = %6.2f \n ' , tau_psi ( i ) , tau_gamma ( j ) , max ( abs ( W_c { i , j } ) ) ) ;
      
        end
        
    end
end



%% Plot

legend_name = cell ( 1 , m_length * n_length ) ;

d = - 1 ;
c =   2 ;

% TROVARE SOLUZIONE ALTERNATIVA A LISTARE LA CELL
for i = 1 : n_length ;
    for j = 1 : m_length ;
        
        k = i + j + d ; 
        legend_name { 1 , k } = [' \tau_{ \psi }=', num2str( tau_psi ( i ) ) ,' \tau_{ \gamma }=', num2str( tau_gamma ( j ) ) ];

    end
        d = d + c ;
end



% Traiettoria sul piano xy

figure ( 1 ) 
hold on
cellfun ( @plot , y_c , x_c )
xlabel ( ' y  [ m ] ' )
ylabel ( ' x  [ m ] ' )
title ( ' Traiettoria sul piano xy ' )
legend ( legend_name , 'Location' , 'northwest'  )

                
% Angolo di rotta

figure ( 2 )
hold on
cellfun ( @plot , t_c , psi_c )
xlabel (' t  [ s ] ' )
ylabel ( ' \psi [ deg ] ' )
title ( ' Andamento nel tempo angolo di rotta \psi' )
legend ( legend_name , 'Location' , 'southeast'  )

% Angolo di rampa

figure ( 3 )
hold on
cellfun ( @plot , t_c , gamma_c )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \gamma  [ deg ] ' )
title ( ' Andamento nel tempo angolo di rampa \gamma ' )
legend ( legend_name , 'Location' , 'southeast'  )

% Fattore di carico

figure ( 4 )
hold on
cellfun ( @plot , t_c , n_c )
% plot ( t1 , n_min * ones ( 1 , length ( t1 ) ) , 'r' )
% plot ( t1 , n_max * ones ( 1 , length ( t1 ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' n_{z}  ' )
title ( ' Andamento nel tempo del fattore di carico n_{z}' )
legend ( legend_name )

% Coefficiente di portanza

figure ( 5 )
hold on
cellfun ( @plot , t_c , Cl_c )
% plot ( t , Cl_min * ones ( 1 , length ( t ) ) , 'r' )
% plot ( t , Cl_max * ones ( 1 , length ( t ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' C_{L}  ' )
title ( ' Andamento nel tempo del coefficiente di portanza C_{L} ' )
legend ( legend_name )

% Grado di ammissione

figure ( 6 )
hold on
cellfun ( @plot , t_c , delta_t_c )
plot ( t , ones ( 1 , length ( t ) ) , 'r' )
plot ( t , zeros ( 1 , length ( t ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \delta_{T}  ' )
title ( ' Andamento nel tempo del grado di ammissione \delta_{T} ' )
axis( [0 t( end ) -1 2 ] )
legend ( legend_name , 'Location' , 'southeast'  )

% Velocità

figure ( 7 )
hold on
cellfun ( @plot , t_c , V_c )
xlabel ( ' t  [ s ] ' )
ylabel ( ' V  [ m/s ] ' )
title ( ' Andamento nel tempo della velocità V ' )
legend ( legend_name , 'Location' , 'southeast'  )

% Varizione del peso in termini percentuali

figure ( 8 )
hold on
cellfun ( @plot , t_c , W_c )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \Delta W / W_{0}   [ % ] ' )
title ( ' Andamento nel tempo della variazione di peso del velivolo \Delta W / W_{0}  ' )
legend ( legend_name )

toc