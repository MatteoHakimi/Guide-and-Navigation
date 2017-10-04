close all 
clc, clear

global P_max eta_prop S g  C_s 

tic

% Definizione variabili e costanti

b        = 11         ;
beta     = 1e-4       ;
g        = 9.81       ;
W        = 14100      ;
S        = 16.2       ;
AR       = b ^ 2 / S  ;
P_max    = 172e3      ;
Cd_0     = 0.022      ;
e        = 0.8        ;
eta_prop = 0.82       ;
C_s      = 0          ;

% Tempo volo

T_v = 200             ;
eps = 1e-3            ;                            % 210 s per 1e-4 meglio non andare oltre POTREBBE SATURARE!
t1  = 0 : eps  : T_v  ;

% Coordinate r_f

x_e = 1e4   ;
y_e = 4e3   ;
h_e = 1.5e3 ;

% Definizione della traiettoria

x = @(t) ( x_e / T_v ) * t ;
y = @(t) y_e * ( sin ( ( pi / 2 ) * ( t / T_v ) ) ) .^ 2 ;
h = @(t) h_e * ( sin ( ( pi / 2 ) * ( t / T_v ) ) ) .^ 2 ;



% Definizione derivate prime traiettoria

x_dot = ( x ( t1 + eps ) - x ( t1 - eps ) ) / ( 2 * eps ) ;
y_dot = ( y ( t1 + eps ) - y ( t1 - eps ) ) / ( 2 * eps ) ;
h_dot = ( h ( t1 + eps ) - h ( t1 - eps ) ) / ( 2 * eps ) ;


% Definizione derivate seconde traiettoria

x_ddot = ( x ( t1 + eps ) + x ( t1 - eps ) - 2 * x ( t1 ) ) / eps ^2 ;
y_ddot = ( y ( t1 + eps ) + y ( t1 - eps ) - 2 * y ( t1 ) ) / eps ^2 ;
h_ddot = ( h ( t1 + eps ) + h ( t1 - eps ) - 2 * h ( t1 ) ) / eps ^2 ;


% Definizione velocità richiesta

V     = sqrt ( x_dot .^2 + y_dot .^2 + h_dot .^2)                    ; 
V_dot = ( x_dot .* x_ddot + y_dot .* y_ddot + h_dot .* h_ddot ) ./ V ;


% Definizione ratei richiesti
gamma_dot    = ( h_ddot -  h_dot .* ( V_dot ./ V ) ) ./ sqrt ( x_dot .^2 + y_dot .^2 ) ;
psi_dot      = ( y_ddot .* x_dot - x_ddot .* y_dot ) ./ ( x_dot .^2  + y_dot .^2 )     ;


% Calcolo delle prestazioni

gamma  = asin ( h_dot ./ V )                                                                                    ;
psi    = atan2 ( y_dot , x_dot )                                                                                ;
phi    = atan2 ( V .* cos ( gamma ) .* psi_dot , V .* gamma_dot + g * cos ( gamma ) )                           ; 
n      = 1 / g * sqrt ( ( V .* psi_dot .* cos ( gamma ) ) .^2 + ( V_dot .* gamma_dot + g * cos ( gamma ) ).^2 ) ;

[ ~ , ~ , ~ , rho ] = atmosisa( h ( t1 ) ) ; 

Cl = 2 * W * n ./ ( S * rho .* V .^2 )         ;
Cd = Cd_0 + Cl .^ 2 / (pi * AR *e)             ;
D  = 1 / 2 * rho .* Cd .* V .^2 * S            ;
T  = D + W * ( sin ( gamma ) + 1 / g * V_dot ) ;  


delta_h = exp( -beta * h ( t1 ) )                  ;
delta_t = T .* V ./ ( eta_prop * P_max * delta_h ) ;

% Condizioni iniziali

V_0      = V ( 1 )        ;
psi_0    = psi ( 1 )      ;
gamma_0  = gamma ( 1 )    ;
x_0      = x ( t1 ( 1 ) ) ;
y_0      = y ( t1 ( 1 ) ) ;
z_0      = h ( t1 ( 1 ) ) ;          
W_0      = W              ;


% Integrazione delle equazioni con Runge Kutta 4° ordine

t_span       = t1                                               ;
X_0          = [ V_0 , psi_0 , gamma_0 , x_0 , y_0 , z_0 , W_0 ];
[ t_v , X ]  = ode45( @Cessna_182t_v2 , t_span , X_0 , [] , t1, delta_t , delta_h , Cl , Cd , phi )                 ;



%% Plot


% Traiettoria sul piano xy

figure ( 1 ) 
hold on
plot ( y ( t1 ) , x ( t1 ) , 'b' )
plot ( X ( : , 5 ) , X ( : , 4 ), 'r' )
xlabel ( ' y  [ m ] ' )
ylabel ( ' x  [ m ] ' )
title ( ' Traiettoria sul piano xy ' )
legend ( ' Traiettoria problema inverso ' , 'Traiettoria Ode45 ' , 'Location' , 'northwest' )

% Traiettoria componente verticale

figure ( 2 )
hold on
plot ( sqrt( y ( t1 ) .^2 + x ( t1 ) .^2 ) , h ( t1 ) , 'b' )
plot ( sqrt( X ( : , 5 ) .^2 + X ( : , 4 ) .^2 ) , X ( : , 6 ) ,'r' )
xlabel ( ' r  [ m ] ' )
ylabel ( ' h  [ m ] ' )
title ( ' Traiettoria componente verticale ' )
legend ( ' Traiettoria problema inverso ' , 'Traiettoria Ode45 ' , 'Location' , 'northwest' )

% Angolo di sbandamento 

figure ( 3 )
plot ( t1 , phi , 'b' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \phi  [ rad ] ' )
title ( ' Andamento nel tempo angolo di sbandamento \phi ' )

% Angolo di rotta

figure ( 4 )
hold on
plot ( t1 , psi , 'b' )
plot ( t1 , X ( : , 2 ) , 'r' )
xlabel (' t  [ s ] ' )
ylabel ( ' \psi [ rad ] ' )
title ( ' Andamento nel tempo angolo di rotta \psi' )
legend ( ' \phi_{P.I.} ' , ' \phi_{Ode45}' , 'Location' , 'northwest' )

% Angolo di rampa

figure ( 5 )
hold on
plot ( t1 , gamma , 'b' )
plot ( t1 , X ( : , 3 ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \gamma  [ rad ] ' )
title ( ' Andamento nel tempo angolo di rampa \gamma ' )
legend ( ' \gamma_{P.I.} ' , ' \gamma_{Ode45}' , 'Location' , 'northwest' )

% Fattore di carico

figure ( 6 )
hold on
plot ( t1 , n , 'b' )
% plot ( t1 , n_min * ones ( 1 , length ( t1 ) ) , 'r' )
% plot ( t1 , n_max * ones ( 1 , length ( t1 ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' n_{z}  ' )
title ( ' Andamento nel tempo del fattore di carico n_{z}' )

% Coefficiente di portanza

figure ( 7 )
hold on
plot ( t1 , Cl , 'b' )
% plot ( t1 , Cl_min * ones ( 1 , length ( t1 ) ) , 'r' )
% plot ( t1 , Cl_max * ones ( 1 , length ( t1 ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' C_{L}  ' )
title ( ' Andamento nel tempo del coefficiente di portanza C_{L} ' )

% Grado di ammissione

figure ( 8 )
hold on
plot ( t1 , delta_t , 'b' )
plot ( t1 , ones ( 1 , length ( t1 ) ) , 'r' )
plot ( t1 , zeros ( 1 , length ( t1 ) ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' \delta_{T}  ' )
title ( ' Andamento nel tempo del grado di ammissione \delta_{T} ' )
axis( [0 T_v -1 2 ] )

% Velocità

figure ( 9 )
hold on
plot ( t1 , V , 'b' )
plot ( t1 , X( : , 1 ) , 'r' )
xlabel ( ' t  [ s ] ' )
ylabel ( ' V  [ m/s ] ' )
title ( ' Andamento nel tempo della velocità V ' )
legend ( ' V_{P.I.} ' , ' V_{Ode45}' , 'Location' , 'northwest' )



toc