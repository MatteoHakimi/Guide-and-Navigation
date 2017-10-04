tic
close all 
clc, clear

global target_los V_star tau_V tau_psi tau_gamma



% Parametri velivolo

b        = 11        ;
S        = 16.2      ;
AR       = b ^ 2 / S ;
P_max    = 172e3     ;
Cd_0     = 0.022     ;
e        = 0.8       ;
eta_prop = 0.82      ;
beta     = 1e-4      ;
g        = 9.81      ;


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


% Posizione target

target_los = 1e3*[ 10 , 10 , 1.5 ] ;


% Velocità richiesta e tempi tau

V_star = 45 ;

tau_V     = 30 ;
tau_psi   = 30 ;
tau_gamma = 30 ;


% Integrazione delle equazioni con Runge Kutta 4° ordine

dt               = 0.01                                              ;
t_span           = 0 : dt : 339                                      ;
X_0              = [ V_0 , psi_0 , gamma_0 , x_0 , y_0 , z_0 ]       ;
[t,X]            = ode45( 'Cessna_182_t_los', t_span , X_0 )         ;


% Calcolo delle prestazioni

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


%% Plot



% Traiettoria sul piano xy

figure ( 1 ) 
hold on
plot ( y , x , 'b' )
xlabel ( ' y  [ m ] ' )
ylabel ( ' x  [ m ] ' )
title ( ' Traiettoria sul piano xy ' )


% Traiettoria componente verticale

figure ( 2 )
hold on
plot ( sqrt( y .^2 + x .^2 ) , z , 'b' )
xlabel ( ' r  [ m ] ' )
ylabel ( ' h  [ m ] ' )
title ( ' Traiettoria componente verticale ' )


% Angolo di sbandamento 

figure ( 3 )
plot ( t , phi , 'b' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \phi  [ deg ] ' )
title ( ' Andamento nel tempo angolo di sbandamento \phi ' )


% Angolo di rotta

figure ( 4 )
hold on
plot ( t , psi , 'b' )
xlabel (' t  [ s ] ' )
ylabel ( ' \psi [ deg ] ' )
title ( ' Andamento nel tempo angolo di rotta \psi' )


% Angolo di rampa

figure ( 5 )
hold on
plot ( t , gamma , 'b' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \gamma  [ deg ] ' )
title ( ' Andamento nel tempo angolo di rampa \gamma ' )


% Fattore di carico

figure ( 6 )
hold on
plot ( t , n , 'b' )
% plot ( t1 , n_min * ones ( 1 , length ( t1 ) ) , 'r' )
% plot ( t1 , n_max * ones ( 1 , length ( t1 ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' n_{z}  ' )
title ( ' Andamento nel tempo del fattore di carico n_{z}' )


% Coefficiente di portanza

figure ( 7 )
hold on
plot ( t , Cl , 'b' )
% plot ( t , Cl_min * ones ( 1 , length ( t ) ) , 'r' )
% plot ( t , Cl_max * ones ( 1 , length ( t ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' C_{L}  ' )
title ( ' Andamento nel tempo del coefficiente di portanza C_{L} ' )


% Grado di ammissione

figure ( 8 )
hold on
plot ( t , delta_t , 'b' )
plot ( t , ones ( 1 , length ( t ) ) , 'r' )
plot ( t , zeros ( 1 , length ( t ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \delta_{T}  ' )
title ( ' Andamento nel tempo del grado di ammissione \delta_{T} ' )
axis( [0 t( end ) -1 2 ] )


% Velocità

figure ( 9 )
hold on
plot ( t , V , 'b' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' V  [ m/s ] ' )
title ( ' Andamento nel tempo della velocità V ' )



toc